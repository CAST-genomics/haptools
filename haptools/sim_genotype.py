# fmt: off
from __future__ import annotations

import os
import re
import glob
import numpy as np
from cyvcf2 import VCF
from collections import defaultdict
from .admix_storage import GeneticMarker, HaplotypeSegment
from .data import GenotypesVCF, GenotypesPLINK
from .transform import GenotypesAncestry


def output_vcf(
        breakpoints, 
        chroms, 
        model_file, 
        variant_file, 
        sampleinfo_file, 
        region, 
        pop_field, 
        sample_field, 
        no_replacement,
        out,
        log,
        chunk_size = None,
    ):
    """
    Takes in simulated breakpoints and uses reference files, vcf and sampleinfo, 
    to create simulated variants output in file: out + .vcf

    Parameters
    ----------
    breakpoints: list[list[HaplotypeSegment]]
        the simulated breakpoints
    chroms: list[str]
        List of chromosomes that were used to simulate
    model_file: str
        file with the following structure. (Must be tab delimited)

        * Header: number of samples, Admixed, {all pop labels}
        * Below: generation number, frac, frac, frac

        For example,

        .. code-block::

                40    Admixed    CEU   YRI
                1       0        0.05  0.95
                2       0.20     0.05  0.75
                
    variant_file: str
        file path that contains samples and respective variants. Can be in
        VCF, BCF, VCF.GZ, or PGEN format.
    sampleinfo_file: str
        file path that contains mapping from sample name in vcf to population
    region: dict(str->str/int/int)
        Dictionary with the keys "chr", "start", and "end" holding chromosome (str),
        start position (int) and end position (int) allowing the simulation process to 
        only allow variants within that region.
    pop_field: boolean
        Flag to determine whether to have the population field in the VCF file output
    sample_field: boolean
        Flag to determine whether to have the sample field in the VCF file output
    no_replacement: boolean
        Flag to determine whether we sample from the reference VCF with or without
        replacement. When True there will be no replacement.
    out: str
        output prefix
    log: log object
        Outputs messages to the appropriate channel.
    chunk_size: int, optional
        The max number of variants to write to a PGEN file together
    """

    log.info(f"Outputting file {out}")

    # details to know
    # vcf file: how to handle samples and which sample is which haplotype block randomly choose out of current population types
    # need to go line by line of the vcf when creating the new vcf

    # read on populations from model file
    mfile = open(model_file, 'r')
    num_samples, *pops = mfile.readline().strip().split()
    pop_dict = {}
    for pop_ind, pop in enumerate(pops):
        pop_dict[pop_ind] = pop
    log.debug(f"Created pop_dict {pop_dict}.")

    # filter sampleinfo so only populations from model file are there
    samplefile = open(sampleinfo_file, 'r')
    pop_sample = defaultdict(list)
    for line in samplefile:
        sample, pop = line.strip().split()
        if pop in pops:
            pop_sample[pop].append(sample)
    log.debug(f"Filtered sample info to limit populations within model file. Populations: {pops}.")

    # check if pops without admixed is same as grabbed populations
    assert len(pops)-1 == len(list(set(pop_sample.keys())))

    # read VCF to get sample names and their order so we can easily call GTs via their index from their sample. 
    #      IE sample HG00097 is index 1 in a list of samples [HG00096 HG00097 HG00098]
    # create sample dictionary that holds sample name to the index in the vcf file for quick access 
    if variant_file.endswith(".pgen"):
        vcf = GenotypesPLINK(variant_file, chunk_size=chunk_size, log=log)
    else:
        vcf = GenotypesVCF(variant_file, log=log)
    
    if not region:
        vcf.read()
    else:
        vcf.read(region=f"{region['chr']}:{region['start']}-{region['end']}")
    vcf.check_missing()

    log.debug(f"Read in variants from {variant_file}")

    sample_dict = {}
    for ind, sample in enumerate(vcf.samples):
        sample_dict[sample] = ind

    # initialize hap_used array which should contain list of lists for each reference sample's genome and what segment has been used for each
    haps_used = [[] for samp in range(len(vcf.samples)*2)]

    log.debug(f"Created index array storing per sample which segment is being processed.")

    # Load populations and samples from sample info and breakpoints as 2 matrices with size variants x samples x 2
    # Use search_sorted() numpy function to find indices in populations/samples in order to determine population and sample info
    ref_vars = vcf.variants
    output_gts = np.empty((int(len(breakpoints)//2), len(vcf.variants), 2), dtype=vcf.data.dtype)
    if not out.endswith(".pgen"):
        if pop_field:
            output_pops = np.empty((int(len(breakpoints)//2), len(vcf.variants), 2), dtype=np.uint8)
        if sample_field:
            output_labels = np.empty((int(len(breakpoints)//2), len(vcf.variants), 2), dtype=object)

    # cover "chr" prefix cases
    if ref_vars["chrom"][0].startswith("chr"):
        cur_chrom = "chr"
    else:
        cur_chrom = ""

    # create samples x variants x 2 matrix of populations
    for hap_ind, haplotype in enumerate(breakpoints):
        cur_var = 0
        ref_gts = np.empty((len(ref_vars), ), dtype=np.uint8)
        if pop_field:
            ref_pops = np.empty((len(ref_vars), ), dtype=np.uint8)
        if sample_field:
            ref_labels = np.empty((len(ref_vars), ), dtype=object)
        # convert vcf variant chroms to ints 
        for chrom in chroms:
            # limit reference vcf variants to current chrom
            ref_vars_chrom = ref_vars["pos"][ref_vars["chrom"] == f"{cur_chrom}{chrom}"]

            # Convert haplotype to numpy array of segment position, populations, samples, and sample haplotypes
            hap_positions, hap_pops, hap_samples_name, hap_samples_ind, hap_sample_haps = \
                    _convert_haplotype(haplotype, chrom, pop_dict, pop_sample, sample_dict, haps_used, no_replacement)

            # if the variant position is = breakpoint end then we consider it part of that bkp
            bkp_pos = np.searchsorted(ref_vars_chrom, hap_positions, side='right')

            # Interval lengths to repeat each population and sample
            sub_array = np.insert(bkp_pos, 0, 0)
            inter_len = np.diff(sub_array)

            # create ref_gts from mapping hap_samples to their respective genotypes 
            ref_sample_inds_chrom = np.repeat(hap_samples_ind, inter_len)

            # select which alleles for each segment
            if no_replacement:
                ref_sample_haps_chrom = np.repeat(hap_sample_haps, inter_len)
            else:
                hap_sample_haps = np.random.randint(2, size=len(hap_samples_ind))
                ref_sample_haps_chrom = np.repeat(hap_sample_haps, inter_len)

            # grab random haplotype from samples and use as gts for our simulated samples
            end_var = cur_var+ref_sample_inds_chrom.shape[0]
            gt_vars = np.arange(cur_var, end_var, 1)

            ref_gts_chrom = vcf.data[ref_sample_inds_chrom,
                                     gt_vars,
                                     ref_sample_haps_chrom]
            ref_gts[cur_var:end_var] = ref_gts_chrom

            if pop_field:
                ref_pops_chrom = np.repeat(hap_pops, inter_len)
                ref_pops[cur_var:end_var] = ref_pops_chrom
            if sample_field:
                ref_labels_chrom = np.repeat(hap_samples_name, inter_len)
                ref_labels[cur_var:end_var] = ref_labels_chrom
            cur_var = end_var

        output_gts[int(hap_ind//2), : , hap_ind % 2] = ref_gts
        if not out.endswith(".pgen"):
            if pop_field:
                output_pops[int(hap_ind//2), : , hap_ind % 2] = ref_pops
            if sample_field:
                output_labels[int(hap_ind//2), : , hap_ind % 2] = ref_labels

    # output vcf header to new vcf file we create
    output_samples = [f"Sample_{hap+1}" for hap in range(int(len(breakpoints)/2))]

    # If PGEN use genotypesPLINK class otherwise use GenotypesAncestry to hold genotypes 
    gts = None
    if not out.endswith(".pgen"):
        # Setup Genotypes class to hold our genotype data
        if pop_field:
            # Initialize Ancestry pops matrix 
            gts = GenotypesAncestry(out, log=log)
            gts.popnum_ancestry = pop_dict
            gts.ancestry = output_pops

        # Setup Genotypes class to hold our genotype data
        if sample_field:
            if gts is None:
                gts = GenotypesAncestry(out, log=log)
                gts.valid_labels = output_labels

        if not pop_field and not sample_field:
            gts = GenotypesVCF(out, log=log)

    else:
        gts = GenotypesPLINK(out, chunk_size=chunk_size, log=log)
    
    gts.samples = output_samples
    gts.variants = vcf.variants
    gts.data = output_gts
    gts.write()
    log.debug("Writing Complete!")
    
    return

def _convert_haplotype(haplotype, chrom, pop_dict, pop_sample, sample_dict, haps_used, no_replacement):
    if chrom == 'X': chrom = 23
    # Do binary search to find beginning of chr in haplotype segments
    hap_start_ind = start_segment(0, int(chrom), haplotype)
    hap_subset = haplotype[hap_start_ind:]
    hap_pos = []
    hap_pops = []
    hap_inds = []
    hap_samples = []
    hap_samples_ind = []

    # collect all segments within chromosome
    for segment in hap_subset:
        if not segment.get_chrom() == int(chrom):
            break

        # Grab reference sample to take variants for current segment
        population = pop_dict[segment.get_pop()]

        # Sample without replacement by keeping track of all segments used for each sample
        if no_replacement:
            # Shuffle samples for each pop if no_replacement to have random selection of samples
            np.random.shuffle(pop_sample[population])
                
            # No segments have been collected yet so start at position 0
            if not hap_pos:
                sample_name, hap_ind = _find_random_sample(
                                                pop_sample[population],
                                                sample_dict,
                                                haps_used,
                                                chrom,
                                                0,
                                                segment.get_end_coord(),
                                                )
            else:
                sample_name, hap_ind = _find_random_sample(
                                                pop_sample[population],
                                                sample_dict,
                                                haps_used,
                                                chrom,
                                                hap_pos[-1]+1,
                                                segment.get_end_coord(),
                                                )
            hap_inds.append(hap_ind)
        else:
            sample_name = np.random.choice(pop_sample[population])
        
        hap_pos.append(segment.get_end_coord())
        hap_pops.append(segment.get_pop())
        hap_samples.append(sample_name)
        hap_samples_ind.append(sample_dict[sample_name])

    return np.asarray(hap_pos, dtype=np.int64), \
           np.asarray(hap_pops, dtype=np.uint8), \
           np.asarray(hap_samples, dtype=object), \
           np.asarray(hap_samples_ind, dtype=np.int32), \
           np.asarray(hap_inds, dtype=np.uint8)

def _find_random_sample(samples, sample_dict, haps_used, chrom, start_coord, end_coord):
    for sample in samples:
        for haplotype in range(2):
            # Determine whether the coord we want to sample is being used
            if _find_coord(haps_used[sample_dict[sample]*2+haplotype], chrom, start_coord, end_coord):
                continue
            else:
                return sample, haplotype

    raise Exception(f"No available sample for the current coords {start_coord}-{end_coord}.")

def _find_coord(cur_hap, chrom, start_coord, end_coord):
    if cur_hap:
        for coords in cur_hap:
            # check for whether the segment has already been used or partially been used
            if chrom == coords[0]:
                if start_coord <= coords[1] < end_coord or start_coord < coords[2] <= end_coord: 
                    return True
    # segment isn't found so add segment to current sample's haplotype
    cur_hap.append((chrom, start_coord, end_coord))
    return False

def _prepare_coords(coords_dir, chroms, region):
    # coord file structure chr variant cMcoord bpcoord
    # NOTE coord files in directory should have chr{1-22, X} in the name
    def numeric_alpha(x):
        chrom = re.search(r'(?<=chr)(X|\d+)', x).group()
        if chrom == 'X':
            return 23
        else:
            return int(chrom)

    # sort coordinate files to ensure coords read are in sorted order
    # remove all chr files not found in chroms list
    all_coord_files = glob.glob(f'{coords_dir}/*.map')
    all_coord_files = [coord_file for coord_file in all_coord_files \
                if re.search(r'(?<=chr)(X|\d+)', coord_file) and \
                   re.search(r'(?<=chr)(X|\d+)', coord_file).group() in chroms]
    all_coord_files.sort(key=numeric_alpha)

    if len(all_coord_files) != len(chroms):
        raise Exception(f"Unable to find all chromosomes {chroms} in map file directory.")
    
    # coords list has form chroms x coords
    coords = []
    for coords_file in all_coord_files:
        file_coords = []
        with open(coords_file, 'r') as cfile:
            prev_coord = None
            for line in cfile:
                # create marker from each line and append to coords
                data = line.strip().split()
                if len(data) != 4:
                    raise Exception(f"Map file contains an incorrect amount of fields {len(data)}. It should contain 4.")
                
                if data[0] == 'X':
                    chrom = 23
                else:
                    chrom = int(data[0])
                gen_mark = GeneticMarker(chrom, float(data[2]), 
                                         int(data[3]), prev_coord)
                prev_coord = gen_mark
                file_coords.append(gen_mark)
        coords.append(file_coords)

    # subset the coordinates to be within the region
    if region:
        start_ind = -1
        for ind, marker in enumerate(coords[0]):
            if marker.get_bp_pos() >= region['start'] and start_ind < 0:
                start_ind = ind

            if marker.get_bp_pos() >= region['end']:
                end_ind = ind+1
                break
            else:
                end_ind = len(coords[0])
        coords = [coords[0][start_ind:end_ind]]

    # Update end coords of each chromosome to have max int as the bp coordinate to
    #    prevent issues with variants in VCF file beyond the specified coordinate
    for chrom_coord in coords:
        chrom_coord[-1].bp_map_pos = np.iinfo(np.int32).max

    # store end coords 
    end_coords = [chrom_coord[-1] for chrom_coord in coords]

    # convert coords to numpy array for easy masking
    max_coords = max([len(chrom_coord) for chrom_coord in coords])
    np_coords = np.zeros((len(coords), max_coords)).astype(object)
    return coords, np_coords, max_coords, end_coords

def simulate_gt(model_file, coords_dir, chroms, region, popsize, log, seed=None):
    """
    Simulate admixed genotypes based on the parameters of model_file.

    Parameters
    ----------
    model_file: str
        File with the following structure. (Must be tab delimited)

        * Header: number of samples, Admixed, {all pop labels}
        * Below: generation number, frac, frac, frac

        For example,

        .. code-block::

                40    Admixed    CEU   YRI
                1       0        0.05  0.95
                2       0.20     0.05  0.75

    coords_dir: str
        Directory containing files ending in .map with genetic map coords in cM used
        for recombination points
    chroms: list[str]
        List of chromosomes to simulate admixture for.
    region: dict()
        Dictionary with the keys "chr", "start", and "end" holding chromosome,
        start position adn end position allowing the simulation process to only 
        within that region.
    popsize: int
        Size of population created for each generation. 
    log: log object
        Outputs messages to the appropriate channel.
    seed: int
        Seed used for randomization.

    Returns
    -------
    num_samples: int
        Total number of samples to output 
    pop_dict: dict(int->str)
        Dictionary that maps populations from their encoded version as integers
        to their population name as a string. ex: {1:CEU, 2:YRI}
    next_gen_samples: list[list[HaplotypeSegment]]
        Each list is a person containing a variable number of Haplotype Segments
        based on how many recombination events occurred throughout the generations
        of ancestors for this person.
    """
    # initialize seed used for breakpoints
    if seed:
        np.random.seed(seed)
        log.info(f"Using seed {seed}")

    # load population samples and labels to be simulated 
    mfile = open(model_file, 'r')
    num_samples, *pops = mfile.readline().strip().split()
    num_samples = int(num_samples)
    pop_dict = {}
    for pop_ind, pop in enumerate(pops):
        pop_dict[pop_ind] = pop

    # Load coordinates to use for simulating
    coords, np_coords, max_coords, end_coords = _prepare_coords(coords_dir, chroms, region)
    
    # precalculate recombination probabilities (given map pos in cM and we want M)
    #     shape: len(chroms) x max number of coords (max so not uneven)
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

    # starting generation is 0
    prev_gen = 0
    next_gen_samples = []

    # iterate over generations in model file
    for gen in mfile:
        # setup population proportions and generations to simulate
        cur_gen, *pop_fracs = gen.strip().split()
        cur_gen = int(cur_gen)
        pop_fracs = np.array(pop_fracs).astype(np.float32) 
        sim_gens = cur_gen - prev_gen
        
        # sim generation
        log.info(f"Simulating generation {prev_gen+1}")
        next_gen_samples = _simulate(popsize, pops, pop_fracs, prev_gen, chroms,
                                     coords, end_coords, recomb_probs, next_gen_samples)

        # simulate remaining generations
        for i in range(1, sim_gens):
            log.info(f"Simulating generation {prev_gen+i+1}")

            # update pop_fracs to have 100% admixture since this generation has not been specified in model file
            pop_fracs = [0]*len(pops)
            pop_fracs[0] = 1

            # simulate next generations using previous generations to sample from for admixture
            next_gen_samples = _simulate(popsize, pops, pop_fracs, prev_gen+i, chroms,
                                         coords, end_coords, recomb_probs, next_gen_samples)

        prev_gen = cur_gen 

    mfile.close()
    return num_samples, pop_dict, next_gen_samples

def write_breakpoints(samples, pop_dict, breakpoints, out, log):
    """
    Write out a subsample of breakpoints to out determined by samples.

    Parameters
    ----------
    samples: int
        Number of samples to output
    pop_dict: dict(int->str)
        Maps population codes in integers to their names. ex: {1:CEU, 2:YRI}
    breakpoints: list[list[HaplotypeSegment]]
        Each list is a person containing a variable number of Haplotype Segments
        based on how many recombination events occurred throughout the generations
        of ancestors for this person.
    out: str
        output prefix used to output the breakpoint file
    log: log object
        Outputs messages to the appropriate channel.

    Returns
    -------
    breakpoints: list[list[HaplotypeSegment]]
        subsampled breakpoints only containing number of samples
    """
    breakpt_file = out + '.bp'
    log.info(f"Outputting breakpoint file {breakpt_file}")

    # randomly sample breakpoints to get the correct amount of samples to output
    breakpoints = np.array(breakpoints, dtype=object)
    breakpoints_ind = np.random.choice(range(breakpoints.shape[0]), size=2*samples, replace=False)
    breakpoints = breakpoints[breakpoints_ind]

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
                pop = pop_dict[segment.get_pop()]
                chrom = segment.get_chrom()
                end_coord = segment.get_end_coord()
                end_pos = segment.get_end_pos()
                output.write(f"{pop}\t{chrom}\t{end_coord}\t{end_pos}\n")
    return breakpoints

def _simulate(samples, pops, pop_fracs, pop_gen, chroms, coords, end_coords, recomb_probs, prev_gen_samples=None):
    """
    Simulate a single generation of creating a population.

    Parameters
    ----------
    samples: int
        Number of samples to output
    pops: list[str]
        List of populations to be simulated
    pop_fracs: list[float]
        Fraction of the populations that contribute 
        ex: 40    Admixed    CEU   YRI
            1       0        0.05  0.95 -> pop_fracs=[0, 0.05, 0.95]
            2       0.20     0.05  0.75 -> pop_fracs=[0.2, 0.05, 0.75]
        Generation 2 has 20% admixed from the prior generation, 5% pure CEU,
        and 75% pure YRI that contribute to creation of generation 2
    pop_gen: int
        Current generation being simulated
    chroms: list[str]
        sorted list of chromosomes used to generate samples. 
    coords: list[list[GeneticMarker]]
        Each list of markers (cM) corresponds to a chromosome in chroms.
        ie if our list of chroms is [3,4,6] then the first list is to chrom 3, second to 4, etc.
    end_coords: list[GeneticMarker]
        List of the last genetic markers (cM) for each chromosome specified in chroms.
        The indices of the list correspond to those in chroms. 
    recomb_probs: Numpy 2d array
        rows = chroms
        cols = number of genetic markers-1
        Holds probabilities for each marker for whether a recombination event will occur.
        prob of event = 1-np.exp(-dist/100) where dist is in cM and is calculated via the
            current and prior genetic markers
    prev_gen_samples: list[list[HaplotypeSegment]], optional
        Prior generation of samples used to choose parents and swap markers when recombination
        events occur. Each list is a person's haplotype of segments having a distinct population label.

    Returns
    -------
    hap_samples: list[list[HaplotypeSegment]]
        Current generation of samples of the same format to prev_gen_samples. 
    """
    # convert chroms to integer and change X to 23
    chroms = [int(chrom) if chrom != 'X' else 23 for chrom in chroms]

    # generate all samples
    hap_samples = []
    
    # pre compute haplotypes and parent population 
    # if there is no previous generation randomly choose population based on frac
    parent_pop = np.random.choice(np.arange(len(pops)), size=samples, p=pop_fracs)

    # If the individual is admixed find parent chromosomes
    haplotypes = np.random.randint(samples, size=2*samples)
    for i, pop in enumerate(parent_pop):
        if not pop:
            # ensure parent haplotypes are not the same
            while haplotypes[2*i] == haplotypes[2*i+1]:
                haplotypes[2*i+1] = np.random.randint(samples)

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

        # generate haplotype blocks over all chromosomes in chroms
        prev_chrom = chroms[0]
        prev_ind = 0
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
                if segments and segments[-1].get_end_coord() == end_coords[prev_ind].get_bp_pos():
                    prev_ind += 1
                    prev_chrom = chroms[prev_ind]
                    start_bp = 0

                # for each chromosome in between recombination events
                # want every chromosome between cur chrom and prev chrom in chroms list
                # find index of cur_chrom since prev_chrom is chrom_ind
                cur_ind = chroms.index(cur_chrom)
                for i in range(cur_ind - prev_ind):
                    # end_bp = end of chromosome
                    end_bp = end_coords[prev_ind+i].get_bp_pos()
                    prev_map_pos = end_coords[prev_ind+i].get_map_pos()

                    # output segments of prev_chrom+i
                    segments.extend(get_segment(p_pop, haps[homolog], 
                                                chroms[prev_ind+i], start_bp, 
                                                end_bp, prev_map_pos,
                                                prev_gen_samples))
                    
                    # change homolog
                    homolog = np.random.randint(2)
                    start_bp = 0

                prev_ind = cur_ind
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
        if not segments:
            start_bp = 0
        elif segments[-1].get_end_coord() == end_coords[prev_ind].get_bp_pos():
            prev_ind += 1
        else:
            start_bp = segments[-1].get_end_coord()+1

        # output remaining chromosomes
        for i in range(len(chroms)-(prev_ind)):
            # end_bp = end of chromosome
            end_bp = end_coords[prev_ind+i].get_bp_pos()
            prev_map_pos = end_coords[prev_ind+i].get_map_pos()

            # output segments of prev_chrom+i
            segments.extend(get_segment(p_pop, haps[homolog], 
                                        chroms[prev_ind+i], start_bp, 
                                        end_bp, prev_map_pos, prev_gen_samples))
            
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

    Parameters
    ----------
    pop: int
        index of population. Can recover population name from pop_dict
    haplotype: int
        index of range [0, len(prev_gen_samples)] to identify the parent haplotype
        to copy segments from
    chrom: int
        chromosome the haplotype segment lies on
    start_coord: int
        starting coordinate from where to begin taking segments from previous
        generation samples
    end_coord: int
        ending coordinate of haplotype segment
    end_pos: float
        ending coordinate in centimorgans
    prev_gen_samples: list[list[HaplotypeSegment]]
        the previous generation simulated used as the parents for the current
        generation

    Returns
    -------
    segments: list[HaplotypeSegment]
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

    Parameters
    ----------
    start: int
        Coordinate in bp for the start of the segment to output
    chrom: int
        Chromosome that the segments lie on. 
    segments: list[HaplotypeSegment]
        List of the hapltoype segments to search from for a starting point.

    Returns
    -------
    mid: int
        Index of the first genetic segment to collect for output. 
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

def validate_params(model, mapdir, chroms, popsize, invcf, sample_info, no_replacement, region=None, only_bp=False):
    # validate model file
    mfile = open(model, 'r')
    num_samples, *pops = mfile.readline().strip().split()

    try:
        num_samples = int(num_samples)
    except:
        raise Exception("Can't convert samples number to an integer.")

    num_pops = len(pops)

    if num_pops < 3:
        raise Exception(f"Invalid number of populations given: {num_pops}. We require at least 2.")

    if num_samples < 1:
        raise Exception("Number of samples is less than 1.")
    
    # ensure the number of pops = number in pop_fracs
    prev_gen = 0
    for gen in mfile:
        cur_gen, *pop_fracs = gen.strip().split()
        try:
            cur_gen = int(cur_gen)
        except:
            raise Exception("Can't convert generation to integer.")

        try:
            pop_fracs = np.array(pop_fracs).astype(np.float32)
        except:
            raise Exception("Can't convert population fractions to type float.")

        sim_gens = cur_gen - prev_gen
        
        if len(pop_fracs) != num_pops:
            raise Exception("Total fractions given to populations do not match number of populations in the header.")
        if sim_gens < 1:
            raise Exception(f"Current generation {cur_gen} - previous generation {prev_gen} = {sim_gens} is less than 1. "
                            "Please ensure the generations given in the first column are correct.")
        if np.absolute(np.sum(pop_fracs)-1) > 1e-6:
            raise Exception(f"Population fractions for generation {cur_gen} do not sum to 1.")

        prev_gen = cur_gen 

    # Check if mapdir is a valid path
    if not os.path.isdir(mapdir):
        raise Exception("Map directory given is not a valid path.")
    
    # validate chroms given are correctly named
    valid_chroms = [str(x) for x in range(1,23)] + ['X']
    for chrom in chroms:
        if chrom not in valid_chroms:
            raise Exception(f"Chromosome {chrom} in the list given is not valid.")

    # Validate mapdir ensuring it contains proper files.
    try:
        all_coord_files = glob.glob(f'{mapdir}/*.map')
        all_coord_files = [coord_file for coord_file in all_coord_files \
            if re.search(r'(?<=chr)(X|\d+)', coord_file) and \
               re.search(r'(?<=chr)(X|\d+)', coord_file).group() in chroms]
    except:
        raise Exception("No valid coordinate files found. Must contain chr{1-22,X} in the file name"
                        " and end in .map")
    
    if not all_coord_files:
        raise Exception("No valid coordinate files found. Must contain chr{1-22,X} in the file name"
                        " and end in .map")
    
    # validate popsize
    if not isinstance(popsize, int):
        raise Exception("Popsize is not an Integer.")
    if popsize <= 0:
        raise Exception("Popsize must be greater than 0.")

    # Ensure popsize is proper size to sample from
    popsize = max(popsize, 10*num_samples)

    # Complete check if we're only outputting a breakpoint
    if only_bp:
        return popsize

    # Collect samples from vcf
    try:
        vcf_samples = None
        if invcf.endswith(".pgen"):
            vcf = GenotypesPLINK(str(invcf))
            vcf.read_samples()
            vcf_samples = vcf.samples
        else:
            vcf = VCF(str(invcf), lazy=True)
            vcf_samples = tuple(vcf.samples)

    except:
        raise Exception("Unable to collect vcf samples.")

    # validate sample_info file (ensure pops given are in model file and samples in vcf file)
    # ensure sample_info col 2 in pops
    total_per_pop = {}
    sample_pops = set()
    for line in open(sample_info, 'r'):
        sample = line.split()[0]
        info_pop = line.split()[1]
        sample_pops.add(info_pop)
        if not total_per_pop.get(info_pop, False):
            total_per_pop[info_pop] = 1
        else:
            total_per_pop[info_pop] += 1

        if sample not in vcf_samples and info_pop in pops:
            raise Exception(f"Sample {sample} from population {info_pop} in sampleinfo file is not present in the vcf file.")
    
    # Ensure that all populations from the model file are listed in the sample info file 
    #    and that there is sufficient enough samples when no_replacement is specified
    for model_pop in pops[1:]:
        if model_pop not in list(sample_pops):
            raise Exception(f"Population {model_pop} in model file is not present in the sample info file.")
        if no_replacement and total_per_pop[model_pop] < num_samples:
            raise Exception(f"Population {model_pop} does not have enough samples to sample without replacement. Please ensure that each population specified has >= total number of samples output: {num_samples}")

    # Ensure that the region parameter can be properly interpreted
    if region:
        if region['start'] > region['end']:
            raise Exception(f"End coordinates in region {region['end']} are less than the starting coordinates {region['start']}.")

    return popsize
