import re
import glob
import time

import numpy as np

from .admix_storage import GeneticMarker, HaplotypeSegment


# TODO at a certain point we are going to need to ensure populations in model file are also in invcf files
#      This is only required for outputting the haplotypes in a vcf 

def simulate_gt(model_file, coords_dir, seed=None):
    """
    Simulate admixed genotypes based on the parameters of model_file. 
    Arguments
        model
            File with the following structure. (Must be tab delimited)
            Header = # samples, Admixed, {all pop labels}
            Below  = generation#, frac, frac
            ex: 40    Admixed    CEU   YRI
                1       0        0.05  0.95
                2       0.20     0.05  0.75
        coords_dir
            Directory containing files ending in .map with genetic map coords 
                in cM used for recombination points
        seed
            Seed used for randomization.
    Return

    """
    # initialize seed used for breakpoints
    if seed:
        np.random.seed(seed)
        print(f"Using seed {seed}")

    # load population samples and labels to be simulated 
    mfile = open(model_file, 'r')
    num_samples, *pops = mfile.readline().strip().split()
    num_samples = int(num_samples)

    # coord file structure chr variant cMcoord bpcoord
    # NOTE coord files in directory should have chr{1-22, X} in the name
    coords = []

    def numeric_alpha(x):
        chrom = re.search(r'(?<=chr)X|\d+', x).group()
        if chrom == 'X':
            return 23
        else:
            return int(chrom)

    # sort coordinate files to ensure coords read are in sorted order
    all_coord_files = glob.glob(f'{coords_dir}/*.map')
    all_coord_files.sort(key=numeric_alpha)

    # coords list has form chroms x coords
    for coords_file in all_coord_files:
        file_coords = []
        with open(coords_file, 'r') as cfile:
            prev_coord = None
            for line in cfile:
                # create marker from each line and append to coords
                data = line.strip().split()
                if data[0] == 'X':
                    chrom = 23
                else:
                    chrom = int(data[0])
                gen_mark = GeneticMarker(chrom, float(data[2]), 
                                         int(data[3]), prev_coord)
                prev_coord = gen_mark
                file_coords.append(gen_mark)
        coords.append(file_coords)

    # store end coords 
    end_coords = [chrom_coord[-1] for chrom_coord in coords]

    # convert coords to numpy array for easy masking
    max_coords = max([len(chrom_coord) for chrom_coord in coords])
    np_coords = np.zeros((len(coords), max_coords)).astype(object)
    
    # precalculate recombination probabilities (given map pos in cM and we want M)
    recomb_probs = -1*np.ones((len(coords), max_coords))
    for chrom, chrom_coords in enumerate(coords):
        prev_map_pos = chrom_coords[0].get_map_pos()
        np_coords[chrom,:len(chrom_coords)] = chrom_coords[:]
        for cind, coord in enumerate(chrom_coords):
            # get current position
            cur_map_pos = coord.get_map_pos()
            dist = cur_map_pos - prev_map_pos
            recomb_probs[chrom,cind] = 1-np.exp(-dist/100)
            prev_map_pos = cur_map_pos
    coords = np_coords

    # number of haplotypes simulated per generation
    haps_per_gen = max(10000, 20 * num_samples)

    # starting generation is 0
    prev_gen = 0
    next_gen_samples = []

    # Time code
    start = time.time()

    # iterate over generations in model file
    for gen in mfile:
        # setup population proportions and generations to simulate
        cur_gen, *pop_fracs = gen.strip().split()
        cur_gen = int(cur_gen)
        pop_fracs = np.array(pop_fracs).astype(np.float) 
        sim_gens = cur_gen - prev_gen
        
        assert sim_gens > 0
        assert np.sum(pop_fracs) == 1

        # sim generation
        print(f"Simulating generation {prev_gen+1}")
        next_gen_samples = _simulate(haps_per_gen, pops, pop_fracs, prev_gen, 
                                     coords, end_coords, recomb_probs, next_gen_samples)

        # simulate remaining generations
        for i in range(1, sim_gens):
            print(f"Simulating generation {prev_gen+i+1}")
            # simulate next generations using previous generations to sample from for admixture
            next_gen_samples = _simulate(haps_per_gen, pops, pop_fracs, prev_gen+i, 
                                         coords, end_coords, recomb_probs, next_gen_samples)

        prev_gen = cur_gen 

    end = time.time()
    print(f"Time elapsed for simulation: {end - start}")

    mfile.close()
    return num_samples, next_gen_samples

def write_breakpoints(samples, breakpoints, out):
    breakpt_file = out + '.bp'
    print(f"Outputting breakpoint file {breakpt_file}")

    # randomly sample breakpoints to get the correct amount of samples to output
    breakpoints = np.array(breakpoints, dtype=object)
    breakpoints = np.random.choice(breakpoints, size=2*samples, replace=False)

    with open(breakpt_file, 'w') as output:
        for ind, sample in enumerate(breakpoints):
            # Get sample number and haplotype number
            haplotype = (ind)%2 + 1
            sample_num = ind//2 + 1

            # write header for sample
            output.write(f"Sample_{sample_num}_{haplotype}\n")

            # write all segments for sample
            for segment in sample:
                # write all segments for current sample
                pop = segment.get_pop()
                chrom = segment.get_chrom()
                end_coord = segment.get_end_coord()
                end_pos = segment.get_end_pos()
                output.write(f"{pop}\t{chrom}\t{end_coord}\t{end_pos}\n")
    return

def _simulate(samples, pops, pop_fracs, pop_gen, coords, end_coords, recomb_probs, prev_gen_samples=None):
    # generate all samples
    hap_samples = []
    
    # pre compute haplotypes and parent population 
    # if there is no previous generation randomly choose population based on frac]
    parent_pop = np.zeros(samples)
    if not prev_gen_samples:
        parent_pop = np.random.choice(pops, size=samples, p=pop_fracs)

    # choose a haplotype to copy
    # if previous generation of samples randomize two haplotypes
    # otherwise only choose a single haplotype
    haplotypes = -1*np.ones(2*samples)
    if not parent_pop[0]:
        # choose two of previous samples as haplotypes to choose from
        haplotypes = np.random.randint(samples, size=2*samples)

        # ensure no two haplotypes are the same
        for hap in range(samples):
            while haplotypes[2*hap] == haplotypes[2*hap+1]:
                haplotypes[2*hap+1] = np.random.randint(samples)

    # generate all samples
    for sample in range(samples):
        segments = []
        p_pop = parent_pop[sample]
        homolog = np.random.randint(2)
        haps = haplotypes[2*sample:2*sample+2]

        # store all probabilities to compare for recombination
        prob_vals = np.random.rand(recomb_probs.shape[0], 
                                   recomb_probs.shape[1])
        recomb_events = prob_vals < recomb_probs
        true_coords = coords[recomb_events]

        # sort true coords by chrom and pos
        def coord_sort(x):
            return (x.get_chrom(), x.get_map_pos())
        true_coords = sorted(true_coords, key=coord_sort)

        # TODO fix so we only fill chroms present in our coords files
        prev_chrom = 1
        for coord in true_coords:
            # information to generate segments
            prev_coord = coord.get_prev_coord()
            cur_chrom = coord.get_chrom()

            # check if we've yet to complete a chrom
            if segments and segments[-1].get_chrom() == prev_chrom:
                start_bp = segments[-1].get_end_coord()+1
            else:
                start_bp = 0

            # swapping chroms so store segments for each chrom we miss in between swap
            if cur_chrom != prev_chrom:
                # check if we've output the end of the prev chrom
                if segments and segments[-1].get_end_coord() == end_coords[prev_chrom-1].get_bp_pos():
                    prev_chrom = prev_chrom + 1
                    start_bp = 0

                # for each chromosome in between recombination events
                for i in range(cur_chrom - prev_chrom):
                    # end_bp = end of chromosome
                    end_bp = end_coords[prev_chrom+i-1].get_bp_pos()
                    prev_map_pos = end_coords[prev_chrom+i-1].get_map_pos()

                    # output segments of prev_chrom+i
                    segments.extend(get_segment(p_pop, haps[homolog], prev_chrom+i,
                                    start_bp, end_bp, prev_map_pos,
                                    prev_gen_samples))
                    
                    # change homolog
                    homolog = np.random.randint(2)
                    start_bp = 0

                prev_chrom = cur_chrom
            
            # get end bp coord and prev map pos since it updates inside swapping chrom
            end_bp = prev_coord.get_bp_pos()
            prev_map_pos = prev_coord.get_map_pos()

            # Store haplotype segments switching homologs
            segments.extend(get_segment(p_pop, haps[homolog], cur_chrom,
                                        start_bp, end_bp, prev_map_pos,
                                        prev_gen_samples))
            homolog = 1-homolog
            prev_chrom = cur_chrom

        # Check if we've output all chromosomes and if not output them
        prev_chrom = cur_chrom
        if segments[-1].get_end_coord() == end_coords[prev_chrom-1].get_bp_pos():
            prev_chrom += 1
        else:
            start_bp = segments[-1].get_end_coord()+1

        for i in range(23 - (prev_chrom-1)):
            # end_bp = end of chromosome
            end_bp = end_coords[prev_chrom+i-1].get_bp_pos()
            prev_map_pos = end_coords[prev_chrom+i-1].get_map_pos()

            # output segments of prev_chrom+i
            segments.extend(get_segment(p_pop, haps[homolog], prev_chrom+i,
                            start_bp, end_bp, prev_map_pos,
                            prev_gen_samples))
            
            # change homolog
            homolog = np.random.randint(2)
            start_bp = 0

        hap_samples.append(segments)
    return hap_samples

def get_segment(pop, haplotype, chrom, start_coord, end_coord, end_pos, prev_gen_samples):
    """
    Create a segment or segments for an individual of the current generation
    using either a population label (>0) or the previous generation's samples if 
    the admix pop type (0) is used. 
    Arguments
        pop - population to generate haplotype segments
        haplotype - index of range [0, len(prev_gen_samples)] to identify
                    the parent haplotype to copy segments from
        chrom - chromosome the haplotype segment lies on
        start_coord - starting coordinate from where to begin taking segments
                      from previous generation samples
        end_coord - ending coordinate of haplotype segment
        end_pos - ending coordinate in centimorgans
        prev_gen_samples - the previous generation simulated used as the parents
                           for the current generation
    Returns
        A list of HaplotypeSegments storing the population type and end coordinate
    """
    # Take from population data not admixed data
    if pop:
        return [HaplotypeSegment(pop, chrom, end_coord, end_pos)]
    
    # Take from previous admixed data
    else:
        segments = []
        prev_gen_segments = prev_gen_samples[haplotype]

        # iterate over haplotype segments to collect relevant ones
        # use binary search to find starting segment to collect information
        start_seg = start_segment(start_coord, chrom, prev_gen_segments)
        for prev_segment in prev_gen_segments[start_seg:]:
            if prev_segment.get_end_coord() >= end_coord or prev_segment.get_chrom() > chrom:
                break

            segments.append(prev_segment)
            
        # If there is not a segment within the boundary use the segment
        #    that spans the boundary respective population
        if not segments:
            out_pop = prev_segment.get_pop()
        else:
            out_pop = segments[-1].get_pop()

        # Append last segment using previous segments population
        segments.append(HaplotypeSegment(out_pop, chrom, end_coord, end_pos))
        return segments

def start_segment(start, chrom, segments):
    """
    Find first segment that is on chrom and its end coordinate is > start via binary search.
    """
    low = 0
    high = len(segments)-1
    mid = 0

    # first segment > start implies segment prior end coord < start
    while low <= high:
        mid = (high+low) // 2
       
        # collect coordinate and chrom information
        cur_coord = segments[mid].get_end_coord()
        cur_chrom = segments[mid].get_chrom()
        if mid == 0:
            prev_coord = -1
            prev_chrom = -1
        else:
            prev_coord = segments[mid-1].get_end_coord()
            prev_chrom = segments[mid-1].get_chrom()
    
        # check if chromosomes match otherwise update
        if chrom == cur_chrom:
            # check for current coords loc
            if cur_coord < start:
                low = mid + 1

            elif cur_coord >= start:
                if prev_chrom < cur_chrom:
                    return mid

                if prev_chrom == cur_chrom and prev_coord < start:
                    return mid
                else:
                    high = mid - 1
                    
            else:
                return len(segments)

        elif chrom < segments[mid].get_chrom():
            high = mid - 1

        else:
            low = mid + 1

    return len(segments)

