import sys
import numpy as np
from pathlib import Path
from haptools.sim_genotype import output_vcf
from haptools.admix_storage import HaplotypeSegment 
DATADIR = Path(__file__).parent.joinpath("data")

def main():
    # load breakpoints file and convert to to breakpoints
    bkp_file = "/storage/mlamkin/data/haptools/5gen.bp"
    model_file = "/storage/mlamkin/data/haptools/models/5gen.dat"
    vcf_file = "/storage/mlamkin/workspace/haptools/tests/data/simple.vcf"
    sampleinfo_file = "/storage/mlamkin/workspace/haptools/tests/data/simpleinfo.tab"
    out = "/storage/mlamkin/data/haptools/test_out"
    breakpoints = []

    with open(bkp_file, 'r') as bkp:
        sample = []
        # create breakpoints list 
        for line in bkp:
            info = line.strip().split()
            if len(info) > 1:
                # gather segments
                segment = HaplotypeSegment(info[0], int(info[1]), 
                                           int(info[2]), float(info[3]))
                sample.append(segment)
            else:
                # write out sample
                if sample:
                    breakpoints.append(sample)
                sample = []

        # end of file write out last sample
        breakpoints.append(sample)

    breakpoints = np.array(breakpoints, dtype=object)

    # evaluate output_vcf
    # TODO create test sampleinfo file and use vcf file from test sets
    output_vcf(breakpoints, model_file, vcf_file, sampleinfo_file, out)
    return


if __name__ == '__main__':
    main()

