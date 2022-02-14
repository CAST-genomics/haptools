[![Snakemake](https://img.shields.io/badge/snakemake-≥6.7.0-brightgreen.svg?style=flat-square)](https://snakemake.bitbucket.io)

### Note: This repository is still under-construction!
Please wait until we have published our first tagged release before using our code.

# haptools

Haptools is a collection of tools for simulating and analyzing genotypes and phenotypes while taking into account haplotype information. It is particularly designed for analysis of individuals with admixed ancestries, although the tools can also be used for non-admixed individuals.

Homepage: https://haptools.readthedocs.io/

## Installation

UNDER CONSTRUCTION

## Haptools utilities

Haptools consists of multiple utilities listed below. Click on a utility to see more detailed usage information.

* [`haptools simgenome`](haptools/simgenotype/README.md): Simulate genotypes for admixed individuals under user-specified demographic histories. 

`haptools simgenome` takes as input a reference set of ancestry-labeled haplotypes and a demographic model and outputs a VCF file with local ancestry information annotated for each variant. It also outputs a list of local ancestry breakpoints which can be visualized using `haptools karyogram`. The output VCF file can be used as input to downstream tools such as `haptools simphenotype` to simulate phenotype information.

* [`haptools simphenotype`](haptools/simphenotype/README.md): Simulate a complex trait, taking into account local ancestry- or haplotype- specific effects.

'haptools simphenotype' takes as input a VCF file and outputs simulated phenotypes for each sample.

* [`haptools karyogram`](haptools/karyogram/README.md): Visualize a "chromosome painting" of local ancestry labels based on breakpoints output by `haptools simgenome`.


## Contributing

If you are interested in contributing to `haptools`, please get in touch by submitting a Github issue or contacting us at mlamkin@ucsd.edu.



