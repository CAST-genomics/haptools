### Note: This repository is still under-construction!
Please wait until we have published our first tagged release before using our code.

# haptools

Haptools is a collection of tools for simulating and analyzing genotypes and phenotypes while taking into account haplotype information. It is particularly designed for analysis of individuals with admixed ancestries, although the tools can also be used for non-admixed individuals.

Homepage: [https://haptools.readthedocs.io/](https://haptools.readthedocs.io/)

## Installation

We have not officially published `haptools` yet, but in the meantime, you can install it directly from our Github repository.
```bash
pip install git+https://github.com/gymrek-lab/haptools.git
```
Installing `haptools` with the "files" extra requirements enables automatic support for a variety of additional file formats, like PLINK2 PGEN files.
```bash
pip install git+https://github.com/gymrek-lab/haptools.git#egg=haptools[files]
````

## Haptools utilities

Haptools consists of multiple utilities listed below. Click on a utility to see more detailed usage information.

* [`haptools simgenotype`](docs/commands/simgenotype.md): Simulate genotypes for admixed individuals under user-specified demographic histories. 

* [`haptools simphenotype`](https://haptools.readthedocs.io/en/latest/commands/simphenotype.html): Simulate a complex trait, taking into account local ancestry- or haplotype- specific effects. `haptools simphenotype` takes as input a VCF file and outputs simulated phenotypes for each sample.

* [`haptools karyogram`](docs/commands/karyogram.md): Visualize a "chromosome painting" of local ancestry labels based on breakpoints output by `haptools simgenotype`.

* [`haptools transform`](https://haptools.readthedocs.io/en/latest/commands/transform.html): Transform a set of genotypes via a list of haplotypes. Create a new VCF containing haplotypes instead of variants.

* [`haptools ld`](https://haptools.readthedocs.io/en/latest/commands/ld.html): Compute Pearson's correlation coefficient between a target haplotype and a set of haplotypes.

Outputs produced by these utilities are compatible with each other. For example
`haptools simgenome` outputs a VCF file with local ancestry information annotated for each variant. The output VCF file can be used as input to `haptools simphenotype` to simulate phenotype information. `haptools simgenome` also outputs a list of local ancestry breakpoints which can be visualized using `haptools karyogram`. 

## Contributing

We gladly welcome any contributions to `haptools`!

Please read [our contribution guidelines](https://haptools.readthedocs.io/en/latest/project_info/contributing.html) and then submit a [Github issue](https://github.com/gymrek-lab/haptools/issues).
