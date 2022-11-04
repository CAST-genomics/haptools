from __future__ import annotations

import os
import re
import sys
import glob
import time
import numpy as np
from cyvcf2 import VCF
from pysam import VariantFile
from collections import defaultdict
from .admix_storage import GeneticMarker, HaplotypeSegment
from .data import GenotypesRefAlt, GenotypesPLINK


def output_vcf(breakpoints, chroms, model_file, vcf_file, sampleinfo_file, region, out):
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
                
    vcf_file: str
        file path that contains samples and respective variants
    sampleinfo_file: str
        file path that contains mapping from sample name in vcf to population
    region: dict()
        Dictionary with the keys "chr", "start", and "end" holding chromosome,
        start position adn end position allowing the simulation process to only 
        within that region.
    out: str
        output prefix
    """

    print(f"Outputting VCF file {out}.vcf")

    # details to know
    # vcf file: how to handle samples and which sample is which haplotype block randomly choose out of current population types
    # need to go line by line of the vcf when creating the new vcf

    # read on populations from model file
    mfile = open(model_file, 'r')
    num_samples, admix, *pops = mfile.readline().strip().split()

    # filter sampleinfo so only populations from model file are there
    samplefile = open(sampleinfo_file, 'r')
    pop_sample = defaultdict(list)
    for line in samplefile:
        sample, pop = line.strip().split()
        if pop in pops:
            pop_sample[pop].append(sample)

    # check if pops are there and if they aren't throw an error
    assert sorted(pops) == sorted(list(set(pop_sample.keys())))

    # read VCF to get sample names and their order so we can easily call GTs via their index from their sample. 
    #      IE sample HG00097 is index 1 in a list of samples [HG00096 HG00097 HG00098]
    # create sample dictionary that holds sample name to the index in the vcf file for quick access 
    vcf = VCF(vcf_file)
    sample_dict = {}
    for ind, sample in enumerate(vcf.samples):
        sample_dict[sample] = ind

    # create index array to store for every sample which haplotype 
    # block we are currently processing and choose what samples 
    # will be used for each haplotype block
    current_bkps = np.zeros(len(breakpoints), dtype=np.int64)
    hapblock_samples = []
    for haplotype in breakpoints:
        hap_samples = []
        for block in haplotype:
            sample_name = np.random.choice(pop_sample[block.get_pop()])
            sample_ind = sample_dict[sample_name]
            hap_samples.append(sample_ind)
        hapblock_samples.append(hap_samples)

    # output vcf header to new vcf file we create
    output_samples = [f"Sample_{hap+1}" for hap in range(int(len(hapblock_samples)/2))]

    # Process
    # Choose starting samples (random choice) in VCF from respective population for each haplotype segment
    #     Also precalculate (random choice) all samples that will be switched too once the current local ancestry block ends.
    # Iterate over VCF and output variants to file(in the beginning write out the header as well) until end of haplotype block for a sample (have to iterate over all samples each time to check)
    # Once VCF is complete we've output everything we wanted
    # VCF output have a FORMAT field where under format is GT:POP and our sample output is GT:POP ie 1|1:YRI|CEU
    # Note: comment out the code below to enable (very experimental!) PGEN support
    # curr_bkps = current_bkps.copy()
    # _write_pgen(breakpoints, chroms, region, hapblock_samples, curr_bkps, output_samples, vcf_file, out+".pgen")
    _write_vcf(breakpoints, chroms, region, hapblock_samples, vcf.samples, current_bkps, output_samples, vcf, out+".vcf")
    return

def _write_vcf(breakpoints, chroms, region, hapblock_samples, vcf_samples, current_bkps, out_samples, in_vcf, out_vcf):
    """
    in_vcf = cyvcf2 variants we are reading in
    out_vcf = output vcf file we output too
    """
    # output vcf file
    write_vcf = VariantFile(out_vcf, mode="w")

    # make sure the header is properly structured with contig names from ref VCF
    for contig in in_vcf.seqnames:
        # remove chr in front of seqname if present and compare
        if contig.startswith('chr'):
            if contig[3:] in chroms:
                write_vcf.header.contigs.add(contig)
        if contig in chroms:
            write_vcf.header.contigs.add(contig)

    try:
        write_vcf.header.add_samples(out_samples)
    except AttributeError:
        print(
            "Upgrade to pysam >=0.19.1 to reduce the time required to create "
            "VCFs. See https://github.com/pysam-developers/pysam/issues/1104",
            file = sys.stderr,
        )
        for sample in out_samples:
            write_vcf.header.add_sample(sample)

    write_vcf.header.add_meta(
        "FORMAT",
        items=[
            ("ID", "GT"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Genotype"),
        ],
    )
    write_vcf.header.add_meta(
        "FORMAT",
        items=[
            ("ID", "POP"),
            ("Number", 2),
            ("Type", "String"),
            ("Description", "Origin Population of each respective allele in GT"),
        ],
    )
    write_vcf.header.add_meta(
        "FORMAT",
        items=[
            ("ID", "SAMPLE"),
            ("Number", 2),
            ("Type", "String"),
            ("Description", "Origin sample and haplotype of each respective allele in GT"),
        ],
    )
    for var in in_vcf:
        # parse chromosome
        chrom = re.search(r'X|\d+', var.CHROM).group()
        if chrom not in chroms: continue

        # limit output to region
        if region:
            if var.start < region["start"]: continue
            if var.end > region["end"]: break
        
        rec = {
            "contig": var.CHROM,
            "start": var.start,
            "stop": var.start + len(var.REF),
            "qual": None,
            "alleles": (var.REF, *var.ALT),
            "id": var.ID,
            "filter": None,
        }
        
        # parse the record into a pysam.VariantRecord
        record = write_vcf.new_record(**rec)

        for hap in range(len(hapblock_samples)):
            sample_num = hap // 2

            # If breakpoint end coord is < current variant update breakpoint
            bkp = breakpoints[hap][current_bkps[hap]]
            if chrom == 'X': 
                chrom = 23
            while bkp.get_chrom() < int(chrom) or (bkp.get_chrom() == int(chrom) and bkp.get_end_coord() < int(var.start)):
                current_bkps[hap] += 1
                bkp = breakpoints[hap][current_bkps[hap]]
            
            var_sample = hapblock_samples[hap][current_bkps[hap]]
            if hap % 2 == 0:
                # store variant
                if hap > 0:
                    record.samples[f"Sample_{sample_num}"]["GT"] = tuple(gt)
                    record.samples[f"Sample_{sample_num}"]["POP"] = tuple(pops)
                    record.samples[f"Sample_{sample_num}"]["SAMPLE"] = tuple(samples)
                    record.samples[f"Sample_{sample_num}"].phased = True
                gt = []
                pops = []
                samples = []
                hap_var = var.genotypes[var_sample][hap % 2]
                gt.append(hap_var)
                pops.append(bkp.get_pop())
                samples.append(vcf_samples[var_sample] + f"-{hap_var}")
            else:
                hap_var = var.genotypes[var_sample][hap % 2]
                gt.append(hap_var)
                pops.append(bkp.get_pop())
                samples.append(vcf_samples[var_sample] + f"-{hap_var}")

        sample_num = hap // 2
        record.samples[f"Sample_{sample_num+1}"]["GT"] = tuple(gt)
        record.samples[f"Sample_{sample_num+1}"]["POP"] = tuple(pops)
        record.samples[f"Sample_{sample_num+1}"]["SAMPLE"] = tuple(samples)
        record.samples[f"Sample_{sample_num+1}"].phased = True

        # write the record to a file
        write_vcf.write(record)
    write_vcf.close()
    return

def _write_pgen(breakpoints, chroms, region, hapblock_samples, current_bkps, out_samples, in_vcf, out):
    """
    in_vcf = GenotypesRefAlt object we are reading in
    out = pgen file we output to
    """
    # initialize input reader
    if in_vcf.endswith(".pgen"):
        in_vcf = GenotypesPLINK(in_vcf)
    else:
        in_vcf = GenotypesRefAlt(in_vcf)
    in_vcf.read(region=f"{region['chr']}:{region['start']}-{region['end']}")
    # TODO: check with someone, do we need to do this QC?
    in_vcf.check_missing(discard_also=True)
    in_vcf.check_biallelic(discard_also=True)
    in_vcf.check_phase()

    # initialize output writer
    gts = GenotypesPLINK(out)
    gts.samples = out_samples
    gts.variants = in_vcf.variants
    gts.data = np.empty((len(out_samples), len(gts.variants), 2), dtype=in_vcf.data.dtype)

    # Now we just fill out gts.data
    # TODO: figure out if there's a way to optimize the following lines of code so that
    # this for loop is performed in numpy rather than in python? This is especially
    # relevant for situations in which the input is PGEN b/c it's probably just as fast
    # to stream the VCF like this otherwise
    for var_idx, var in enumerate(gts.variants):
        # parse chromosome
        chrom = re.search(r'X|\d+', var["chrom"]).group()
        if chrom not in chroms: continue
        if chrom == 'X':
            chrom = 23
        for hap in range(len(hapblock_samples)):
            sample_num = hap // 2
            # If breakpoint end coord is < current variant update breakpoint
            bkp = breakpoints[hap][current_bkps[hap]]
            # Note: for some reason, var.start in _write_vcf() is always equal to var["pos"]-1 in _write_pgen()
            # We should probably investigate at some point
            while bkp.get_chrom() < int(chrom) or (bkp.get_chrom() == int(chrom) and bkp.get_end_coord() < int(var["pos"])-1):
                current_bkps[hap] += 1
                bkp = breakpoints[hap][current_bkps[hap]]
            var_sample = hapblock_samples[hap][current_bkps[hap]]
            if hap % 2 == 0:
                # store variant
                if hap > 0:
                    gts.data[sample_num-1, var_idx] = tuple(gt)
                gt = []
            gt.append(int(in_vcf.data[var_sample, var_idx, hap % 2]))
        sample_num = hap // 2
        gts.data[sample_num, var_idx] = tuple(gt)

    gts.write()

def simulate_gt(model_file, coords_dir, chroms, region, popsize, seed=None):
    """
    Simulate admixed genotypes based on the parameters of model_file.

    Parameters
    ----------
    model: str
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
    seed: int
        Seed used for randomization.

    Returns
    -------
    num_samples: int
        Total number of samples to output 
    next_gen_samples: list[list[HaplotypeSegment]]
        Each list is a person containing a variable number of Haplotype Segments
        based on how many recombination events occurred throughout the generations
        of ancestors for this person.
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
    def numeric_alpha(x):
        chrom = re.search(r'(?<=chr)(X|\d+)', x).group()
        if chrom == 'X':
            return 23
        else:
            return int(chrom)

    # sort coordinate files to ensure coords read are in sorted order
    # remove all chr files not found in chroms list
    all_coord_files = glob.glob(f'{coords_dir}/*.map')
    if region:
        try:
            all_coord_files = [coord_file for coord_file in all_coord_files \
                        if f"chr{region['chr']}" in coord_file]
        except:
            raise Exception(f"Unable to find region chromosome {region['chr']} in map file directory.")
    else:
        all_coord_files = [coord_file for coord_file in all_coord_files \
                       if re.search(r'(?<=chr)(X|\d+)', coord_file).group() in chroms]
        all_coord_files.sort(key=numeric_alpha)
    
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

            if marker.get_bp_pos() >= region['end'] and coords[0][ind-1].get_bp_pos() < region['end']:
                end_ind = ind+1
                break
        coords = [coords[0][start_ind:end_ind]]

    # store end coords 
    end_coords = [chrom_coord[-1] for chrom_coord in coords]

    # convert coords to numpy array for easy masking
    max_coords = max([len(chrom_coord) for chrom_coord in coords])
    np_coords = np.zeros((len(coords), max_coords)).astype(object)
    
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
        print(f"Simulating generation {prev_gen+1}")
        next_gen_samples = _simulate(popsize, pops, pop_fracs, prev_gen, chroms,
                                     coords, end_coords, recomb_probs, next_gen_samples)

        # simulate remaining generations
        for i in range(1, sim_gens):
            print(f"Simulating generation {prev_gen+i+1}")

            # update pop_fracs to have 100% admixture since this generation has not been specified in model file
            pop_fracs = [0]*len(pops)
            pop_fracs[0] = 1

            # simulate next generations using previous generations to sample from for admixture
            next_gen_samples = _simulate(popsize, pops, pop_fracs, prev_gen+i, chroms,
                                         coords, end_coords, recomb_probs, next_gen_samples)

        prev_gen = cur_gen 

    mfile.close()
    return num_samples, next_gen_samples

def write_breakpoints(samples, breakpoints, out):
    """
    Write out a subsample of breakpoints to out determined by samples.

    Parameters
    ----------
    samples: int
        Number of samples to output
    breakpoints: list[list[HaplotypeSegment]]
        Each list is a person containing a variable number of Haplotype Segments
        based on how many recombination events occurred throughout the generations
        of ancestors for this person.
    out: str
        output prefix used to output the breakpoint file

    Returns
    -------
    breakpoints: list[list[HaplotypeSegment]]
        subsampled breakpoints only containing number of samples
    """
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
                    segments.extend(get_segment(p_pop, pops, haps[homolog], 
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
            segments.extend(get_segment(p_pop, pops, haps[homolog], cur_chrom,
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
            segments.extend(get_segment(p_pop, pops, haps[homolog], 
                                        chroms[prev_ind+i], start_bp, 
                                        end_bp, prev_map_pos, prev_gen_samples))
            
            # change homolog
            homolog = np.random.randint(2)
            start_bp = 0

        hap_samples.append(segments)
    return hap_samples

def get_segment(pop, str_pops, haplotype, chrom, start_coord, end_coord, end_pos, prev_gen_samples):
    """
    Create a segment or segments for an individual of the current generation
    using either a population label (>0) or the previous generation's samples if 
    the admix pop type (0) is used. 

    Parameters
    ----------
    pop
        index of population corresponding to the population in str_pops
    str_pops
        array of population names
    haplotype
        index of range [0, len(prev_gen_samples)] to identify the parent haplotype
        to copy segments from
    chrom
        chromosome the haplotype segment lies on
    start_coord
        starting coordinate from where to begin taking segments from previous
        generation samples
    end_coord
        ending coordinate of haplotype segment
    end_pos
        ending coordinate in centimorgans
    prev_gen_samples
        the previous generation simulated used as the parents for the current
        generation

    Returns
    -------
    list
        A list of HaplotypeSegments storing the population type and end coordinate
    """
    # Take from population data not admixed data
    if pop:
        return [HaplotypeSegment(str_pops[pop], chrom, end_coord, end_pos)]
    
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
    segments: list[HaplotypeSegments]
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

def validate_params(model, mapdir, chroms, popsize, invcf, sample_info, region=None, only_bp=False):
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
            if re.search(r'(?<=chr)(X|\d+)', coord_file).group() in chroms]
    except:
        raise Exception("Could not parse map directory files.")
    
    if not all_coord_files:
        raise Exception("No valid coordinate files found. Must contain chr{1-22,X} in the file name.")
    
    # validate popsize
    if not isinstance(popsize, int):
        raise Exception("Popsize is not an Integer.")
    if popsize <= 0:
        raise Exception("Popsize must be greater than 0.")

    # Ensure popsize is proper size to sample from
    popsize = max(popsize, 10*num_samples)

    # Complete check if we're only outputting a breakpoint
    if only_bp:
        return

    # Collect samples from vcf
    try:
        vcf = VCF(invcf)
        vcf_samples = vcf.samples
    except:
        raise Exception("Unable to collect vcf samples.")

    # validate sample_info file (ensure pops given are in model file and samples in vcf file)
    # ensure sample_info col 2 in pops
    sample_pops = set()
    for line in open(sample_info, 'r'):
        sample = line.split()[0]
        info_pop = line.split()[1]
        sample_pops.add(info_pop)

        if sample not in vcf_samples:
            raise Exception(f"Sample {sample} in sampleinfo file is not present in the vcf file.")
    
    # Ensure that all populations from the model file are listed in the sample info file
    for model_pop in pops[1:]:
        if model_pop not in list(sample_pops):
            raise Exception(f"Population {model_pop} in model file is not present in the sample info file.")

    # Ensure that the region parameter can be properly interpreted
    if region:
        if region['start'] > region['end']:
            raise Exception(f"End coordinates in region {region['end']} are less than the starting coordinates {region['start']}.")

    return popsize
