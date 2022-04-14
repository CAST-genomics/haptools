### Note: This repository is still under-construction!
Please wait until we have published our first tagged release before using our code.

# haptools

Haptools is a collection of tools for simulating and analyzing genotypes and phenotypes while taking into account haplotype information. It is particularly designed for analysis of individuals with admixed ancestries, although the tools can also be used for non-admixed individuals.

Homepage: https://haptools.readthedocs.io/

## Installation

UNDER CONSTRUCTION

## Haptools utilities

Haptools consists of multiple utilities listed below. Click on a utility to see more detailed usage information.

* [`haptools simgenome`](docs/commands/simgenotype.md): Simulate genotypes for admixed individuals under user-specified demographic histories. 

* [`haptools simphenotype`](docs/commands/simphenotype.md): Simulate a complex trait, taking into account local ancestry- or haplotype- specific effects. `haptools simphenotype` takes as input a VCF file and outputs simulated phenotypes for each sample.

* [`haptools karyogram`](docs/commands/karyogram.md): Visualize a "chromosome painting" of local ancestry labels based on breakpoints output by `haptools simgenome`.

Outputs produced by these utilities are compatible with each other. For example
`haptools simgenome` outputs a VCF file with local ancestry information annotated for each variant. The output VCF file can be used as input to `haptools simphenotype` to simulate phenotype information. `haptools simgenome` also outputs a list of local ancestry breakpoints which can be visualized using `haptools karyogram`. 

## Contributing

If you are interested in contributing to `haptools`, please get in touch by submitting a Github issue or contacting us at mlamkin@ucsd.edu.



