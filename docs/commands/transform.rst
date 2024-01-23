.. _commands-transform:


transform
=========

Transform a set of genotypes via a list of haplotypes. Create a new VCF containing haplotypes instead of variants.

The ``transform`` command takes as input a set of **phased** genotypes and a list of haplotypes and outputs a set of haplotype *pseudo-genotypes*, where each haplotype is encoded as a bi-allelic variant record in the output. In other words, each sample will have a genotype of ``0|0``, ``1|0``, ``0|1``, or ``1|1`` indicating whether each of their two chromosome copies contains the alleles of a haplotype.

.. figure:: https://github.com/CAST-genomics/haptools/assets/23412689/fb3accd9-4b15-4ba7-a09c-022b405aa26f
  :figwidth: 600
  :align: center
  :alt: Transforming genotypes via haplotypes

Users may also specify an ancestral population label for each haplotype. See the :ref:`ancestry section <commands-transform-input-ancestry>` for more details.

Usage
~~~~~
.. code-block:: bash

	haptools transform \
	--region TEXT \
	--sample SAMPLE --sample SAMPLE \
	--samples-file FILENAME \
	--id ID --id ID \
	--ids-file FILENAME \
	--chunk-size INT \
	--discard-missing \
	--ancestry \
	--output PATH \
	--verbosity [CRITICAL|ERROR|WARNING|INFO|DEBUG|NOTSET] \
	GENOTYPES HAPLOTYPES

Input
~~~~~
Genotypes must be specified in VCF and haplotypes must be specified in the :doc:`.hap file format </formats/haplotypes>`.

Alternatively, you may specify genotypes in PLINK2 PGEN format. Just use the appropriate ".pgen" file extension in the input. See the documentation for genotypes in :ref:`the format docs <formats-genotypesplink>` for more information.

.. _commands-transform-input-ancestry:

Ancestry
--------
If your ``.hap`` file contains :ref:`an "ancestry" extra field <formats-haplotypes-extrafields-transform>` and your VCF contains a "POP" format field (as output by :doc:`simgenotype </commands/simgenotype>`), you should specify the ``--ancestry`` flag.
This will enable us to match the population labels of each haplotype against those in the genotypes output by :doc:`simgenotype </commands/simgenotype>`.
In other words, a sample is said to contain a haplotype only if all of the alleles of the haplotype are labeled with the haplotype's ancestry.

.. figure:: https://github.com/CAST-genomics/haptools/assets/23412689/f00553c9-8a82-4b9e-9929-042da6d95f02
  :figwidth: 600
  :align: center
  :alt: Transforming via ancestry labels

Alternatively, you may specify a :doc:`breakpoints file </formats/breakpoints>` accompanying the genotypes file. It must have the same name as the genotypes file but with a ``.bp`` file ending. If such a file exists, ``transform`` will ignore any "POP" format fields in the genotypes file and instead obtain the ancestry labels from the breakpoints file. This is primarily a speed enhancement, since it's faster to load ancestral labels from the breakpoints file.

Output
~~~~~~
Transform outputs *psuedo-genotypes* in VCF, but you may request genotypes in PLINK2 PGEN format, instead. Just use the appropriate ".pgen" file extension in the output path. See the documentation for genotypes in :ref:`the format docs <formats-genotypesplink>` for more information.

Examples
~~~~~~~~
.. code-block:: bash

	haptools transform tests/data/simple.vcf.gz tests/data/simple.hap

Let's try transforming just two samples and let's output to PGEN format:

.. code-block:: bash

	haptools transform -o output.pgen -s HG00097 -s NA12878 tests/data/apoe.vcf.gz tests/data/apoe4.hap

To get progress information, increase the verbosity to "INFO":

.. code-block:: bash

	haptools transform --verbosity INFO -o output.vcf.gz tests/data/example.vcf.gz tests/data/basic.hap.gz

To match haplotypes as well as their ancestral population labels, use the ``--ancestry`` flag:

.. code-block:: bash

	haptools transform --ancestry tests/data/simple-ancestry.vcf tests/data/simple.hap

All files used in these examples are described :doc:`here </project_info/example_files>`.


Detailed Usage
~~~~~~~~~~~~~~

.. click:: haptools.__main__:main
   :prog: haptools
   :nested: full
   :commands: transform
