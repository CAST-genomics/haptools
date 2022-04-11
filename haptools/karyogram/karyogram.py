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


def splitstr(option, opt, value, parser):
  return(setattr(parser.values, option.dest, value.split(',')))

# TODO determine whether its better to plot bp or cM and about updating output so it outputs both
# TODO update so we plot in cM because plots are much cleaner
# TODO figure out centromere
def plot_karyogram(sample_file, title, centromeres, out, sample_name="Sample_1", chrX=False, colors=None):
    """
    Arguments
        sample_file - contains all samples, chrom, pop, and hap block location
        title - plot title
        centromeres - centromeres location file
        out - output file to save figure
        sample_name - sample to plot on karyogram
        chrX - include chromosome X? 
        colors - colors for respective populations

    """
    #read in bed files and get individual name
    samples = []
    sample = []
    pop_order = []
    plot_sample = False

    with open(sample_file,'r') as sample_file:
        for line in sample_file:
            line = line.strip().split('\t')

            # check if header
            if len(line) == 1:
                # new header write out haplotype
                if sample:
                    samples.append(sample)

                if sample_name in line[0]:
                    plot_sample = True
                    sample = []
                    continue
                else:
                    plot_sample = False
            
                # Already collected our two haplotypes
                if not plot_sample and len(samples) == 2:
                    break
            
            if plot_sample:
                if not line[0] in pop_order:
                    pop_order.append(line[0])

                if not sample:
                    start = 0.0001
                else:
                    start = sample[-1]['end'] + 0.0001
                hap_block = {'pop': line[0], 'chrom': int(line[1]), 
                             'start': start, 'end': float(line[-1])}
                sample.append(hap_block)

    #define plotting space
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(-5,300)
    if chrX:
      ax.set_ylim(24,0)
    else:
      ax.set_ylim(23,0)
    plt.xlabel('Genetic position (cM)')
    plt.ylabel('Chromosome')
    plt.title(title)
    if chrX:
      plt.yticks(range(1,24))
      yticks = range(1,23)
      yticks.append('X')
      ax.set_yticklabels(yticks)
    else:
      plt.yticks(range(1,23))

    #define colors
    def hex_to_rgb(value):
        value = value.lstrip('#')
        lv = len(value)
        return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

    bmap = brewer2mpl.get_map('Set1', 'qualitative', 4)
    if not colors:
      colors=bmap.mpl_colors
      colors.append((0,0,0))

    #define centromeres
    # TODO figure out what to do with centromeres
    """
    centro = open(args.centromeres)
    centromeres = {}
    for line in centro:
        line = line.strip().split()
        if chrX and line[0] == 'X':
          line[0] = '23'
        centromeres[line[0]] = line
    """

    #plot rectangles
    for info in samples[0]:
        try:
          plot_rects(info['pop'], info['chrom'], info['start'], 
                     info['end'], 'A', pop_order, colors, ax)
        except ValueError: #flexibility for chrX
          plot_rects(info['pop'], 23, info['start'], 
                     info['end'], 'A', pop_order, colors, ax)
    for info in samples[1]:
        try:
          plot_rects(info['pop'], info['chrom'], info['start'], 
                     info['end'], 'B', pop_order, colors, ax)
        except ValueError: #flexibility for chrX
          plot_rects(info['pop'], 23, info['start'], 
                     info['end'], 'B', pop_order, colors, ax)

    #write a legend
    p = []
    for i in range(len(pop_order)):
        p.append(plt.Rectangle((0, 0), 1, 1, color=colors[i]))
    p.append(plt.Rectangle((0, 0), 1, 1, color='k'))
    labs = list(pop_order)
    labs.append('UNK')
    leg = ax.legend(p, labs, loc=4, fancybox=True)
    leg.get_frame().set_alpha(0)

    #get rid of annoying plot features
    spines_to_remove = ['top', 'right']
    for spine in spines_to_remove:
        ax.spines[spine].set_visible(False)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')

    fig.savefig(out)

def plot_rects(anc, chrom, start, stop, hap, pop_order, colors, ax):
    # TODO update so we can work with centromeres but rn dont use
    """  
    centro_coords = map(float, centromeres[str(chrom)])
    if len(centro_coords) == 3: #acrocentric chromosome
        mask = [
        (centro_coords[1]+2,chrom-0.4), #add +/- 2 at the end of either end
        (centro_coords[2]-2,chrom-0.4),
        (centro_coords[2]+2,chrom),
        (centro_coords[2]-2,chrom+0.4),
        (centro_coords[1]+2,chrom+0.4),
        (centro_coords[1]-2,chrom),
        (centro_coords[1]+2,chrom-0.4)
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
    
    else: #need to write more complicated clipping mask with centromere masked out
        mask = [
        (centro_coords[1]+2,chrom-0.4), #add +/- 2 at the end of either end
        (centro_coords[2]-2,chrom-0.4),
        (centro_coords[2]+2,chrom+0.4),
        (centro_coords[3]-2,chrom+0.4),
        (centro_coords[3]+2,chrom),
        (centro_coords[3]-2,chrom-0.4),
        (centro_coords[2]+2,chrom-0.4),
        (centro_coords[2]-2,chrom+0.4),
        (centro_coords[1]+2,chrom+0.4),
        (centro_coords[1]-2,chrom),
        (centro_coords[1]+2,chrom-0.4)
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
    """   
    if hap == 'A': #bed_a ancestry goes on top
        verts = [
            (float(start), chrom), #left, bottom
            (float(start), chrom + 0.4), #left, top
            (float(stop), chrom + 0.4), #right, top
            (float(stop), chrom), #right, bottom
            (0, 0), #ignored
        ]
    else: #bed_b ancestry goes on bottom
        verts = [
            (float(start), chrom - 0.4), #left, bottom
            (float(start), chrom), #left, top
            (float(stop), chrom), #right, top
            (float(stop), chrom - 0.4), #right, bottom
            (0, 0), #ignored
        ]

    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]
    
    clip_path = Path(verts, codes)
    if anc in pop_order:
        col=mcol.PathCollection([clip_path],facecolor=colors[pop_order.index(anc)], linewidths=0)
    else:
        col=mcol.PathCollection([clip_path],facecolor=colors[-1], linewidths=0)
    #if 'clip_mask' in locals():
    #    col.set_clip_path(clip_mask, ax.transData)
    ax.add_collection(col)


