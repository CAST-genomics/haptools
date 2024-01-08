"""
This script is inspired by Alicia Martin's karyogram code
originally published here: 
https://github.com/armartin/ancestry_pipeline/blob/master/plot_karyogram.py
"""

import os
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.collections as mcol
import numpy as np
import os
import sys

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

PADDING = 0.4  # PADDING between chromosomes


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
    if "X" in chrom:
        return 23
    if "Y" in chrom:
        return 24
    if chrom.startswith("chr"):
        return int(chrom[3:])
    else:
        return int(chrom)


def GetHaplotypeBlocks(bp_file, sample_name, centromeres_file=None):
    """
    Extract haplotype blocks for the desired sample
    from the bp file

    Parameters
    ----------
    bp_file : str
       Path to .bp file with breakpoints
    sample_name : str
       Sample ID to extract
    centromeres_file : str, optional
        If not None then use the chromosome ends listed to extend
        chromosomes to proper end coordinates

    Returns
    -------
    sample_blocks : list[list[hap_blocks]]
       each hap_block is a dictionary with keys
       'pop', 'chrom', 'start', 'end'
    """
    sample_blocks = []  # blocks for the two copies
    parsing_sample = False  # keep track of if we're in the middle of parsing a sample
    blocks = []  # keep track of current blocks

    if not os.path.exists(bp_file):
        sys.stderr.write("ERROR: Breakpoints file %s not found.\n" % bp_file)
        sys.exit(1)

    with open(bp_file, "r") as f:
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
                if len(blocks) == 0 or (blocks[-1]["chrom"] != GetChrom(line[1])):
                    start = 0.0001
                else:
                    start = blocks[-1]["end"] + 0.0001
                hap_block = {
                    "pop": line[0],
                    "chrom": GetChrom(line[1]),
                    "start": start,
                    "end": float(line[-1]),
                }
                blocks.append(hap_block)
    # Check if we still need to add the last block
    # Happens if the sample is at the end of the file
    if parsing_sample:
        sample_blocks.append(blocks)

    # If the centromeres file is given update the blocks so the last
    # block of the chromosome is extended to its actual end in cM
    if centromeres_file:
        # collect all actual chromosome ends
        chrom_ends = {}
        with open(centromeres_file, "r") as cfile:
            for line in cfile:
                chrom_data = line.strip().split()
                chrom_ends[GetChrom(chrom_data[0])] = float(chrom_data[-1])

        # update current haplotype tracts to end at actual chrom ends
        for hap, block in enumerate(sample_blocks):
            prev_chrom = block[0]["chrom"]
            for tind, tract in enumerate(block):
                cur_chrom = tract["chrom"]
                if cur_chrom != prev_chrom:
                    sample_blocks[hap][tind - 1]["end"] = chrom_ends[prev_chrom]
                prev_chrom = cur_chrom
            # Update last chromosome since chromosome won't change
            sample_blocks[hap][tind - 1]["end"] = chrom_ends[prev_chrom]

    return sample_blocks


def GetCmRange(sample_blocks):
    """
    Get the min and max cM coordinates from the sample_blocks

    Parameters
    ----------
    sample_blocks : list[list[hap_blocks]]
       each hap_block is a dictionary with keys
       'pop', 'chrom', 'start', 'end'

    Returns
    -------
    min_val, max_val : float, float
       min_val is the minimum coordinate
       max_val is the maximum coordinate
    """
    min_val = np.inf
    max_val = -1 * np.inf
    for i in range(len(sample_blocks)):
        for sb in sample_blocks[i]:
            if sb["start"] < min_val:
                min_val = sb["start"]
            if sb["end"] > max_val:
                max_val = sb["end"]
    return min_val, max_val


def GetPopList(sample_blocks):
    """
    Get a list of populations in the sample_blocks

    Parameters
    ----------
    sample_blocks : list[list[hap_blocks]]
       each hap_block is a dictionary with keys
       'pop', 'chrom', 'start', 'end'

    Returns
    -------
    poplist : list[str]
       list of populations represented in the blocks
    """
    poplist = set()
    for i in range(len(sample_blocks)):
        for sb in sample_blocks[i]:
            poplist.add(sb["pop"])
    return list(poplist)


def GetChromOrder(sample_blocks):
    """
    Get a list of chroms in sorted order

    Parameters
    ----------
    sample_blocks : list[list[hap_blocks]]
       each hap_block is a dictionary with keys
       'pop', 'chrom', 'start', 'end'

    Returns
    -------
    chroms : list[int]
       list of chromsomes in sorted order
    """
    chroms = set()
    for i in range(len(sample_blocks)):
        for sb in sample_blocks[i]:
            chroms.add(sb["chrom"])
    chroms = list(chroms)
    chroms.sort()
    return chroms


def PlotKaryogram(
    bp_file, sample_name, out_file, log, centromeres_file=None, title=None, colors=None
):
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
    log: log object
        Outputs messages to the appropriate channel.
    centromeres_file : str, optional
       Path to file with centromere coordinates.
       Format: chrom, chromstart_cm, centromere_cm, chromend_cm
       If None, no centromere and telomere locations are shown
    title : str, optional
       Plot title. If None, no title is annotated
    colors : dict(str->str), optional
       Dictionary of colors to use for each population
       If not set, reasonable defaults are used.
       In addition to strings, you can specify RGB or RGBA tuples.
    """
    # Parse haplotype blocks from the bp file for the
    # specified sample
    log.info("Collecting Haplotype Blocks...")
    sample_blocks = GetHaplotypeBlocks(bp_file, sample_name, centromeres_file)
    if len(sample_blocks) == 0:
        sys.stderr.write("ERROR: no haplotype blocks identified for %s. " % sample_name)
        sys.stderr.write(
            "Make sure %s_1 and %s_2 are in the .bp file.\n"
            % (sample_name, sample_name)
        )
        sys.exit(1)

    # Extract metadata about the blocks
    min_cm, max_cm = GetCmRange(sample_blocks)
    chrom_order = GetChromOrder(sample_blocks)
    pop_list = GetPopList(sample_blocks)

    # Set up the figure for plotting
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel("Genetic position (cM)")
    ax.set_ylabel("Chromosome")
    if title is not None:
        ax.set_title(title)
    ax.set_yticks(range(len(chrom_order)))
    ax.set_yticklabels(chrom_order)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")

    # Set up colors
    if colors is None:
        log.info("Colors not given. Setting up colors...")
        num_pops = len(pop_list)
        cmap = plt.cm.get_cmap("rainbow", num_pops)
        colors = dict(zip(pop_list, cmap(list(range(num_pops)))))

    # Optionally, plot centromeres/telomeres
    clipmask_perchrom = None
    if centromeres_file is not None:
        log.info("Centromeres present, adding into figure...")
        clipmask_perchrom = GetCentromereClipMask(centromeres_file, chrom_order)

    # Plot the actual haplotype blocks
    for i in range(2):
        for info in sample_blocks[i]:
            PlotHaplotypeBlock(
                info, i, chrom_order, colors, ax, clipmask_perchrom=clipmask_perchrom
            )

    # Write a legend
    legend_elements = []
    for i in range(len(pop_list)):
        legend_elements.append(plt.Rectangle((0, 0), 1, 1, color=colors[pop_list[i]]))
    legend_elements.append(plt.Rectangle((0, 0), 1, 1, color="k"))
    leg = ax.legend(legend_elements, pop_list, loc=4, fancybox=True)
    leg.get_frame().set_alpha(0)

    # Make sure axis limits are wide enough
    ax.set_xlim(min_cm - 5, max_cm + 5)
    ax.set_ylim(len(chrom_order) + 3, -3)

    fig.savefig(out_file)
    log.info(f"Karyogram Complete! Saved to {out_file}")


def GetCentromereClipMask(centromeres_file, chrom_order):
    """
    Get clipping mask for the centromeres and telomeres

    Parameters
    ----------
    centromeres_file : str, optional
       Path to file with centromere coordinates.
       Format: chrom, chromstart_cm, centromere_cm, chromend_cm
       If None, no centromere and telomere locations are shown
    chrom_order : list[int]
       chromosomes in sorted order

    Returns
    -------
    clipmask_perchrom : dict[str, matplotlib.Path]
       Clip region for telomeres/centromeres for each chromosome
    """
    clipmask_perchrom = {}
    if not os.path.exists(centromeres_file):
        sys.stderr.write("ERROR: centromeres file %s not found.\n" % centromeres_file)
        sys.exit(1)
    with open(centromeres_file, "r") as f:
        for line in f:
            items = line.strip().split()
            chrom = GetChrom(items[0])
            # If centromeres file has a chomosome not in blocks file pass on it
            if chrom not in chrom_order:
                continue
            chrom_ind = chrom_order.index(chrom)
            centro_coords = [float(item) for item in items[1:]]
            if len(centro_coords) == 2:  # acrocentric
                mask = [
                    (
                        centro_coords[0] + 2,
                        chrom_ind - PADDING,
                    ),  # add +/- 2 at the end of either end
                    (centro_coords[1] - 2, chrom_ind - PADDING),
                    (centro_coords[1] + 2, chrom_ind),
                    (centro_coords[1] - 2, chrom_ind + PADDING),
                    (centro_coords[0] + 2, chrom_ind + PADDING),
                    (centro_coords[0] - 2, chrom_ind),
                    (centro_coords[0] + 2, chrom_ind - PADDING),
                ]
                mask_codes = [
                    Path.MOVETO,
                    Path.LINETO,
                    Path.CURVE3,
                    Path.LINETO,
                    Path.LINETO,
                    Path.CURVE3,
                    Path.LINETO,
                ]
                clip_mask = Path(vertices=mask, codes=mask_codes)
            else:
                mask = [
                    (
                        centro_coords[0] + 2,
                        chrom_ind - PADDING,
                    ),  # add +/- 2 at the end of either end
                    (centro_coords[1] - 2, chrom_ind - PADDING),
                    (centro_coords[1] + 2, chrom_ind + PADDING),
                    (centro_coords[2] - 2, chrom_ind + PADDING),
                    (centro_coords[2] + 2, chrom_ind),
                    (centro_coords[2] - 2, chrom_ind - PADDING),
                    (centro_coords[1] + 2, chrom_ind - PADDING),
                    (centro_coords[1] - 2, chrom_ind + PADDING),
                    (centro_coords[0] + 2, chrom_ind + PADDING),
                    (centro_coords[0] - 2, chrom_ind),
                    (centro_coords[0] + 2, chrom_ind - PADDING),
                ]

                mask_codes = [
                    Path.MOVETO,
                    Path.LINETO,
                    Path.LINETO,
                    Path.LINETO,
                    Path.CURVE3,
                    Path.LINETO,
                    Path.LINETO,
                    Path.LINETO,
                    Path.LINETO,
                    Path.CURVE3,
                    Path.LINETO,
                ]
                clip_mask = Path(vertices=mask, codes=mask_codes)
            clipmask_perchrom[chrom] = clip_mask
    return clipmask_perchrom


def PlotHaplotypeBlock(block, hapnum, chrom_order, colors, ax, clipmask_perchrom=None):
    """
    Plot a haplotype block on the axis

    Parameters
    ----------
    block : dict
       dictionary with keys
       'pop', 'chrom', 'start', 'end'
    hapnum : int
       0 or 1 for the two haplotypes
    chrom_order : list[int]
       chromosomes in sorted order
    colors : dict[str, str], optional
       Dictionary of colors to use for each population
       If not set, reasonable defaults are used.
       In addition to strings, you can specify RGB or RGBA tuples.
    ax : matplotlib axis to use for plotting
    clipmask_perchrom : dict[str, matplotlib.Path], optional
       Clip region for telomeres/centromeres for each chromosome
    """
    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]
    chrom_coord = chrom_order.index(block["chrom"])
    start = float(block["start"])
    stop = float(block["end"])

    verts = [
        (start, chrom_coord - hapnum * PADDING),
        (start, chrom_coord + (1 - hapnum) * PADDING),
        (stop, chrom_coord + (1 - hapnum) * PADDING),
        (stop, chrom_coord - hapnum * PADDING),
        (0, 0),
    ]
    clip_path = Path(verts, codes)
    col = mcol.PathCollection([clip_path], facecolor=colors[block["pop"]], linewidths=0)

    # Optionally, deal with centromeres
    if clipmask_perchrom:
        col.set_clip_path(clipmask_perchrom[block["chrom"]], ax.transData)

    ax.add_collection(col)
