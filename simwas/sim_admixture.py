import vcf
import numpy as np
from classes import GeneticMarker, HaplotypeSegment

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

    # iterate over generations in model file
    for gen in mfile:
        # setup population proportions and generations to simulate
        cur_gen, *pop_fracs = gen.strip.split()
        cur_gen = int(cur_gen)
        pop_fracs = np.array(pop_fracs).astype(np.float) 
        sim_gens = cur_gen - pre_gen
        
        assert sim_gens > 0
        assert np.sum(pop_fracs) == 1

        # sim first generation
        prev_gen_samples = _simulate(num_samples, pop_fracs, pre_gen, coords)

        # simulate remaining generations
        for generations in range(1, sim_gens):
            # simulate next generations using previous generations to sample from for admixture
            next_gen_samples = _simulate(num_samples, pop_fracs, pre_gen+i, coords, prev_gen_samples)
            prev_gen_samples = next_gen_samples

        prev_gen = cur_gen 

    mfile.close()
    return


def _simulate(samples, pop_fracs, pop_gen, coords, prev_gen_samples=None):
    # generate all samples
    for sample in range(samples):
        prev_chrom = -1
        prev_map_pos = 0.0

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
            cur_chrom = coord.get_chrom()
            cur_map_pos = coord.get_map_pos()
            if prev_chrom < 0:
                prev_chrom = cur_chrom
            if cur_chrom != prev_chrom:
                # TODO end of chromosome save haplotype segments
                
                prev_chrom = cur_chrom
                homolog = np.random.randint(2)
                prev_map_pos = cur_map_pos
        
            # check if a recombination event occurs
            # TODO function or logic for recombination event
            dist = cur_map_pos - prev_map_pos
            if recomb:
                # TODO save haplotype segments since were moving to a new homolog after this breakpoint
                homolog = 1-homolog
            prev_map_pos = cur_map_pos
            prev_chrom = cur_chrom

        # TODO record remaining segments

        # TODO output debug message for individual simulated

# TODO remove
def test():
    model_file = './AA.dat' 
    coords_file = '/storage/mlamkin/data/simwas/plink.chr10.GRCh38.map'
    simulate_gt(model_file, coords_file)

if __name__ == '__main__':
    test()
