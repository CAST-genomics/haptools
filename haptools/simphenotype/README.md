# Haptools simphenotype

Haptools simphenotype simulates a complex trait, taking into account haplotype- or local-ancestry- specific effects as well as traditional variant-level effects. It takes causal effects and genotypes as input and outputs simulated phenotypes.

Usage is modeled based on the [GCTA GWAS Simulation](https://yanglab.westlake.edu.cn/software/gcta/#GWASSimulation) utility.

## Usage

Below is a basic `haptools simphenotype` command:

```
haptools simphenotype \
  --vcf <gzipped vcf file> \
  --hap <gzipped hap file> \
  --out <outprefix> \
  < --simu-qt | --simu-cc > \
  [simulation options]
```

Required parameters:

* `--vcf <string>`: A bgzipped, tabix-indexed, phased VCF file. If you are simulating local-ancestry effects, the VCF file must contain the `FORMAT/LA` tag included in output of `haptools simgenotype`. See [haptools file formats](../../docs/project_info/haptools_file_formats.rst) for more details.

* `--hap <string>`: A bgzipped, tabix-indexed HAP file, which specifies causal effects. This is a custom format described in more detail in [haptools file formats](../../docs/project_info/haptools_file_formats.rst). The HAP format enables flexible specification of a range of effect types including traditional variant-level effects, haplotype-level effects, associations with repeat lengths at short tandem repeats, and interaction of these effects with local ancestry labels. See [Examples](#examples) below for detailed examples of how to specify effects.

* `--out <string>`: Prefix to name output files.

* `--simu-qt` or `simu-cc` indicate whether to simulate a quantitative or case control trait. One of these two options must be specified.

Additional parameters:

* `--simu-rep <int>`: Number of phenotypes to simulate. Default: 1.

* `--simu-hsq <float>`: Trait heritability. Default: 0.1.

* `--simu-k <float>`: Disease prevalence (for case-control). Default: 0.1

## Output files

`haptools simphenotypes` outputs the following files:

* `<outprefix>.phen`: Based on the phenotype files accepted by [plink](https://www.cog-genomics.org/plink/1.9/input#pheno). Tab-delimited file with one row per sample. The first and second columns give the sample ID. The third column gives the simulated phenotype. If `--simu-rep` was set to greater than 1, additional columns are output for each simulated trait. Example file:

```
HG00096	HG00096	0.058008375906919506
HG00097	HG00097	0.03472768471423458
HG00099	HG00099	-0.20850127859705808
HG00100	HG00100	-0.21206803352471154
HG00101	HG00101	0.3157913763428323
```

* `<outprefix>.par`: summarizes the frequencies and effects of simulated haplotypes. The columns are: haplotype ID (from the HAP file), haplotype frequency, and effect. Example file:

```
Haplotype	Frequency	Effect
H-001	0.6	-0.2
```

<a name="examples"></a>
## Examples

1. Simulate a single haplotype-effect based on a 2 SNP haplotype:

```
haptools simphenotype \
  --vcf tests/data/simple.vcf.gz \
  --hap tests/data/simple.hap.gz \
  --out test \
  --simu-qt --simu-hsq 0.8 --simu-rep 1
```

based on this HAP file (available in `tests/data`)

```
H-001	1	10114	10116	1:10114:T:C,1:10116:A:G	T,G	*	-0.2
```