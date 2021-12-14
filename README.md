[![Snakemake](https://img.shields.io/badge/snakemake-≥6.7.0-brightgreen.svg?style=flat-square)](https://snakemake.bitbucket.io)

### Note: This repository is still under-construction!
Please wait until we have published our first tagged release before using our code.

# admixtools
Simulate phenotypes for fine-mapping. Use real variants to simulate real, biological LD patterns.
The pipeline also uses the results of the simulation to test several fine-mapping methods, including FINEMAP and SuSiE.

# download
Execute the following command.
```
git clone https://github.com/aryarm/simulate_gwas
```
You can also download example data for the pipeline. See [the config file](config/config.yml) for links and instructions.

# setup
The pipeline is written as a Snakefile which can be executed via [Snakemake](https://snakemake.readthedocs.io). For reproduciblity, we recommend installing the version that we used (6.7.0):
```
conda create -n snakemake -c conda-forge --no-channel-priority 'bioconda::snakemake==6.7.0'
```
`snakemake` will [automatically install all dependencies](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management) of the pipeline upon its first execution using `conda`.

# execution
1. Activate snakemake via `conda`:
    ```
    conda activate snakemake
    ```
2. Execute the pipeline on the example data

    Locally:
    ```
    ./run.bash &
    ```
    __or__ on a TORQUE cluster:
    ```
    qsub run.bash
    ```
### Output
All output of the pipeline will be placed in a new directory (`out/`, by default).
Log files describing the output of the pipeline will also be created there. The `log` file contains a basic description of the progress of each rule, while the `qlog` file is more detailed and will contain any errors or warnings.

### Executing the pipeline on your own data
You must modify [the config.yaml file](config/config.yml) to specify paths to your data before you perform step 2 above. Currently, the pipeline is configured to run on our example data.

## If this is your first time using Snakemake
Snakemake is a tool for creating and running pipelines. The contents of the pipeline are specified as a set of steps (called _rules_) in a special file, called a Snakefile. Each step consists of inputs, outputs, a shell command that uses the inputs to create the outputs, and a conda environment with software dependencies for the shell command. When you execute Snakemake, it will assemble the rules into a pipeline according to each rule's inputs and outputs, along with the outputs that you request from Snakemake.

We recommend that you run `snakemake --help` to learn about Snakemake's options. For example, to check that the pipeline will be executed correctly before you run it, you can call Snakemake with the `-n -p -r` flags. This is also a good way to familiarize yourself with the steps of the pipeline and their inputs and outputs (the latter of which are inputs to the first rule in each workflow -- ie the `all` rule).

Note that Snakemake will not recreate output that it has already generated, unless you request it. If a job fails or is interrupted, subsequent executions of Snakemake will just pick up where it left off. This can also apply to files that *you* create and provide in place of the files it would have generated.

By default, the pipeline will automatically delete some files it deems unnecessary (ex: unsorted copies of a VCF). You can opt to keep these files instead by providing the `--notemp` flag to Snakemake when executing the pipeline.

