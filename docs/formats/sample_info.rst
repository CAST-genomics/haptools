.. _formats-sample_info:


Sample Info
===========

1000 Genomes sample_info file format
------------------------------------
Within the subcommand ``haptools simgenotype`` we use a file to map samples in the 
reference to their population listed in the model file. This code is expected to work
out of the box with 1000 genomes data and we have pre-computed this mapping file as 
well as given the command to how to compute it if you desire another as well.

``cut -f 1,4 igsr-1000\ genomes\ on\ grch38.tsv | sed '1d' | sed -e 's/ /\t/g' > 1000genomes_sampleinfo.tsv``

See ``example-files/1000genomes_sampleinfo.tsv`` for an example of the 1000genomes 
GRCh38 samples mapped to their subpopulations.
