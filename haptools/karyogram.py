import os
import argparse
import pylab
import brewer2mpl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.collections as mcol
import numpy as np

def GetChrom(chrom):
    """
    Extract a numerical chromosome

    Parameters
    ----------
    chrom : str
       Chromosome string

    Returns
    -------
    chrom : int
       Integer-value for the chromosome
       X gets set to 23
    """
    if chrom.startswith("chr"):
        if "X" not in chrom:
            return int(chrom[3:])
        else: return 23
    else: return int(chrom)

def GetHaplotypeBlocks(bp_file, sample_name):
    """
    Extract haplotype blocks for the desired sample
    from the bp file

    Parameters
    ----------
    bp_file : str
       Path to .bp file with breakpoints
    sample_name : str
       Sample ID to extract
    
    Returns
    -------
    sample_blocks : list of [hap_blocks]
       each hap_block is a dictionary with keys
       'pop', 'chrom', 'start', 'end'
    """
    sample_blocks = [] # blocks for the two copies

    parsing_sample = False # keep track of if we're in the middle of parsing a sample
    blocks = [] # keep track of current blocks

    with open(bp_file, 'r') as f:
        for line in f:
            line = line.strip().split()
            if len(line) == 1:
                assert line[0].endswith("_1") or line[0].endswith("_2")
                # Check if we're done parsing a previous sample
                # If so, add the blocks and reset
                if parsing_sample:
                    sample_blocks.append(blocks.copy())

                # Check if we're done
                if len(sample_blocks) == 2:
                    parsing_sample = False
                    break

                # Check if we should start processing the next sample
                if sample_name == "_".join(line[0].split("_")[:-1]):
                    blocks = []
                    parsing_sample = True
                    continue
                else:
                    parsing_sample = False

            # If we're in the middle of parsing a sample, add the block
            if parsing_sample:
                if len(blocks) == 0 or (blocks[-1]['chrom'] != GetChrom(line[1])):
                    start = 0.0001
                else:
                    start = blocks[-1]['end'] + 0.0001
                    assert(float(line[-1]) > start) # Check the file is in sorted order
                hap_block = {'pop': line[0], 'chrom': GetChrom(line[1]), 
                             'start': start, 'end': float(line[-1])}
                blocks.append(hap_block)
    # Check if we still need to add the last block
    # Happens if the sample is at the end of the file
    if parsing_sample:
        sample_blocks.append(blocks)

    return sample_blocks

def GetCmRange(sample_blocks):
    """
    Get the min and max cM coordinates from the sample_blocks
    Parameters
    ----------
    sample_blocks : list of [hap_blocks]
       each hap_block is a dictionary with keys
       'pop', 'chrom', 'start', 'end'

    Returns
    -------
    min_val, max_val : float, float
       min_val is the minimum coordinate
       max_val is the maximum coordinate
    """
    min_val = np.inf
    max_val = -1*np.inf
    for i in range(len(sample_blocks)):
        for sb in sample_blocks[i]:
            if sb['start'] < min_val: min_val = sb['start']
            if sb['end'] > max_val: max_val = sb['end']
    return min_val, max_val

def GetPopList(sample_blocks):
    """
    Get a list of populations in the sample_blocks
    ----------
    sample_blocks : list of [hap_blocks]
       each hap_block is a dictionary with keys
       'pop', 'chrom', 'start', 'end'

    Returns
    -------
    poplist : list of str
       list of populations represented in the blocks
    """
    poplist = set()
    for i in range(len(sample_blocks)):
        for sb in sample_blocks[i]:
            poplist.add(sb['pop'])
    return list(poplist)

def GetChromOrder(sample_blocks):
    """
    Get a list of chroms in sorted order
    Parameters
    ----------
    sample_blocks : list of [hap_blocks]
       each hap_block is a dictionary with keys
       'pop', 'chrom', 'start', 'end'

    Returns
    -------
    chroms : list of int
       list of chromsomes in sorted order
    """
    chroms = set()
    for i in range(len(sample_blocks)):
        for sb in sample_blocks[i]:
            chroms.add(sb['chrom'])
    chroms = list(chroms)
    chroms.sort()
    return chroms

# TODO figure out centromeres
def PlotKaryogram(bp_file, sample_name, out_file,
        centromeres_file=None, title=None, colors=None):
    """
    Plot a karyogram based on breakpoints output by haptools simgenotypes

    Parameters
    ----------
    bp_file : str
       Path to .bp file with breakpoints
    sample_name : str
       Sample ID to plot
    out_file : str
       Name of output file
    centromeres_file : str, optional
       Path to bed file with centromere coordinates.
       If None, no centromere locations are shown
       (NOTE: centromeres not yet implemented)
    title : str, optional
       Plot title. If None, no title is annotated
    colors : dict of str:str, optional
       Dictionary of colors to use for each population
       If not set, reasonable defaults are used.
       In addition to strings, you can specify RGB or RGBA tuples.
    """
    # Parse haplotype blocks from the bp file for the 
    # specified sample
    sample_blocks = GetHaplotypeBlocks(bp_file, sample_name)

    # Extract metadata about the blocks
    min_cm, max_cm = GetCmRange(sample_blocks)
    chrom_order = GetChromOrder(sample_blocks)
    pop_list = GetPopList(sample_blocks)

    # Set up the figure for plotting
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(min_cm-5, max_cm+5)
    ax.set_ylim(len(chrom_order)+3, -3)
    ax.set_xlabel('Genetic position (cM)')
    ax.set_ylabel('Chromosome')
    if title is not None: ax.set_title(title)
    ax.set_yticks(range(len(chrom_order)))
    ax.set_yticklabels(chrom_order)    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    # Set up colors
    if colors is None:
        num_colors = len(pop_list)
        if num_colors < 3: num_colors = 3
        if num_colors > 9: num_colors = 9
        bmap = brewer2mpl.get_map('Set1', 'qualitative', num_colors)
        colors = dict(zip(pop_list, bmap.mpl_colors))

    # Plot the actual haplotype blocks
    for i in range(2):
        for info in sample_blocks[i]:
            PlotHaplotypeBlock(info, i, chrom_order, colors, ax)

    #write a legend
    p = []
    for i in range(len(pop_list)):
        p.append(plt.Rectangle((0, 0), 1, 1, color=colors[pop_list[i]]))
    p.append(plt.Rectangle((0, 0), 1, 1, color='k'))
    labs = pop_list
    leg = ax.legend(p, labs, loc=4, fancybox=True)
    leg.get_frame().set_alpha(0)
    fig.savefig(out_file)

def PlotHaplotypeBlock(block, hapnum, chrom_order, colors, ax):
    """
    Plot a haplotype block on the axis

    Parameters
    ----------
    block : dictionary with keys
       'pop', 'chrom', 'start', 'end'
    hapnum : int
       0 or 1 for the two haplotypes
    chrom_order : list of int
       chromosomes in sorted order
    colors : dict of str:str, optional
       Dictionary of colors to use for each population
       If not set, reasonable defaults are used.
       In addition to strings, you can specify RGB or RGBA tuples.
    ax : matplotlib axis to use for plotting
    """
    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]
    chrom_coord = chrom_order.index(block['chrom'])
    start = float(block['start'])
    stop = float(block['end'])

    padding = 0.4
    verts = [
        (start, chrom_coord - hapnum*padding),
        (start, chrom_coord + (1-hapnum)*padding),
        (stop, chrom_coord + (1-hapnum)*padding),
        (stop, chrom_coord - hapnum*padding),
        (0,0)
    ]
    clip_path = Path(verts, codes)
    col = mcol.PathCollection([clip_path], facecolor=colors[block['pop']], linewidths=0)
    ax.add_collection(col) 