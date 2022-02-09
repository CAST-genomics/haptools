import vcf
import numpy as np
from .admix_storage import GeneticMarker, HaplotypeSegment

# TODO at a certain point we are going to need to ensure populations in model file are also in invcf files
#      This is only required for outputting the haplotypes in a vcf 

def simulate_gt(model_file, coords_file, seed=None):
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
        coords_file
            File containing genetic map coords in cM used for recombination points
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
    num_samples, not_used, *pops = mfile.readline().strip().split()

    # coord file structure chr variant cMcoord bpcoord
    coords = []
    with open(coords_file, 'r') as cfile:
        for line in cfile:
            # create marker from each line and append to coords
            data = line.strip().split()
            gen_mark = GeneticMarker(int(data[0]), float(data[2]), int(data[3])) 
            coords.append(gen_mark)

    # number of haplotypes simulated per generation
    haps_per_gen = max(10000, 20 * num_samples)

    # starting generation is 0
    prev_gen = 0
    next_gen_samples = []

    # iterate over generations in model file
    for gen in mfile:
        # setup population proportions and generations to simulate
        cur_gen, *pop_fracs = gen.strip.split()
        cur_gen = int(cur_gen)
        pop_fracs = np.array(pop_fracs).astype(np.float) 
        sim_gens = cur_gen - prev_gen
        
        assert sim_gens > 0
        assert np.sum(pop_fracs) == 1

        # sim generation
        next_gen_samples = _simulate(num_samples, pop_fracs, pre_gen, coords, next_gen_samples)

        # simulate remaining generations
        for generations in range(1, sim_gens):
            # simulate next generations using previous generations to sample from for admixture
            next_gen_samples = _simulate(num_samples, pop_fracs, pre_gen+i, coords, next_gen_samples)

        prev_gen = cur_gen 

    mfile.close()
    return next_gen_samples


def _simulate(samples, pop_fracs, pop_gen, coords, prev_gen_samples=None):
    # generate all samples
    hap_samples = []
    for sample in range(samples):
        segments = []
        prev_chrom = -1
        prev_map_pos = 0.0
        start_bp = 0
        end_bp = 0  

        # if there is no previous generation randomly choose population based on frac
        parent_pop = 0
        if not prev_gen_samples:
            parent_pop = np.random.choice(np.arange(len(pop_fracs)), p=pop_fracs)

        # choose a haplotype to copy
        # if previous generation of samples randomize two haplotypes
        # otherwise only choose a single haplotype
        haplotypes = -1*np.ones(2)
        if not parent_pop:
            # choose two of previous samples as haplotypes to choose from
            haplotypes = np.random.choice(np.arange(samples), size=2, replace=False)
        
        # haplotypes only matters when pop == 0 otherwise we are sampling from a population which has a label 
        # choose our starting chromosome
        homolog = np.random.randint(2)

        # iterate over all coords in the coords file

        for coord in coords:
            # get current chromosome and position
            cur_chrom = coord.get_chrom()
            cur_map_pos = coord.get_map_pos()
            
            # start of new segment used for storing previous generation samples  
            if segments:
                start_bp = segments[-1].get_end_coord() + 1

            if prev_chrom < 0:
                prev_chrom = cur_chrom

            # Ended chromosome write segment
            if cur_chrom != prev_chrom:
                # Store haplotype
                segments.extend(get_segment(parent_pop, haplotypes[homolog], 
                                            start_bp, end_bp, prev_gen_samples))
                
                # update homolog, previous chrom, and previous position
                prev_chrom = cur_chrom
                homolog = np.random.randint(2)
                prev_map_pos = cur_map_pos
                continue
        
            # check if a recombination event occurs and if so write segment and swap homolog
            dist = cur_map_pos - prev_map_pos
            recomb_prob = 1-np.exp(-dist)
            if np.random.choice([0, 1], 1, p=[1-recomb_prob, recomb_prob]):
                # Store haplotype segments switching homologs
                segments.extend(get_segment(parent_pop, haplotypes[homolog], 
                                            start_bp, end_bp, prev_gen_samples))
                homolog = 1-homolog

            # update positions and chromosome
            prev_map_pos = cur_map_pos
            prev_chrom = cur_chrom
            end_bp = coord.get_bp_pos()

        # Record remaining segments
        segments.extend(get_segment(parent_pop, haplotypes[homolog], start_bp,
                                    end_bp, prev_gen_samples))

        # append segments of 1 sample
        hap_samples.append(segments)
    return hap_samples

def get_segment(pop, haplotype, start_coord, end_coord, prev_gen_samples):
    """
    Create a segment or segments for an individual of the current generation
    using either a population label (>0) or the previous generation's samples if 
    the admix pop type (0) is used. 
    Arguments
        pop - population to generate haplotype segments
        haplotype - index of range [0, len(prev_gen_samples)] to identify
                    the parent haplotype to copy segments from
        start_coord - starting coordinate from where to begin taking segments
                      from previous generation samples
        end_coord - ending coordinate of haplotype segment
        prev_gen_samples - the previous generation simulated used as the parents
                           for the current generation
    Returns
        A list of HaplotypeSegments storing the population type and end coordinate
    """
    # Take from population data not admixed data
    if pop:
        return [HaplotypeSegment(pop, end_coord)]
    
    # Take from previous admixed data
    else:
        segments = []
        prev_gen_segments = prev_gen_samples[haplotype]

        # iterate over haplotype segments to collect relevant ones
        for prev_segment in prev_gen_segments:
            # check if current segment is within boundaries
            if prev_segment.get_end_coord() < start_coord:
                continue
            if prev_segment.get_end_coord() > end_coord:
                break

            segments.append(prev_segment)
            
        # Append last segment using previous segments population
        segments.append(HaplotypeSegment(segments[-1].get_pop(), end_coord))
        return segments

# TODO remove
def test():
    model_file = '/storage/mlamkin/data/simwas/AA.dat' 
    coords_file = '/storage/mlamkin/data/simwas/plink.chr10.GRCh38.map'
    simulate_gt(model_file, coords_file)

if __name__ == '__main__':
    test()
