# Haptools karyogram

`haptools karyogram` takes as input a breakpoints file (e.g. as output by `haptools simgenotype`) and a sample name, and plots a karyogram depicting local ancestry tracks.


## Basic usage

```
haptools karyogram \
  --bp tests/data/5gen.bp \
  --sample Sample_1 \
  --out karyogram.png
```

See details of the breakpoints file [here](simgenotype.md). If you specify `--sample $SAMPLE`, the breakpoints file must have breakpoints for `$SAMPLE_1` and `$SAMPLE_2` (the two haplotypes of `$SAMPLE`).

## Additional options

You may also specify the following options:

* `--centromeres <FILE>`: path to a file describing the locations of chromosome ends and centromeres. An example file is given here: `tests/data/centromeres_hg19.txt`. The columns are: chromosome, chrom_start, centromere, chrom_end. For acrocentric chromosomes, the centromere field is ommitted. This file format was taken from [here](https://github.com/armartin/ancestry_pipeline).
* `--colors "pop1:color1,pop2:color2..."`: You can optionally specify which colors should be used for each population. If colors are not given, the script chooses reasonable defaults.

## Example command

The following example can be run with files in this repository:

```
haptools karyogram --bp tests/data/5gen.bp --sample Sample_1 \
       --out test_karyogram.png --centromeres tests/data/centromeres_hg19.txt \
       --colors 'CEU:blue,YRI:red'
```

This will output a file `test_karyogram.png`. The example is shown below.

![Example karyogram](../images/test_karyogram.png)