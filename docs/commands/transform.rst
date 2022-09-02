.. _commands-transform:


transform
=========

Transform a set of genotypes via a list of haplotypes. Create a new VCF containing haplotypes instead of variants.

The ``transform`` command takes as input a set of genotypes in VCF and a list of haplotypes (specified as a :doc:`.hap file </formats/haplotypes>` without any extra fields) and outputs a set of haplotype "genotypes" in VCF.

You may also specify genotypes in PLINK2 PGEN format. Just use the appropriate ".pgen" file extension in the input and/or output. See the documentation for genotypes in :ref:`the format docs <formats-genotypesplink>` for more information.

If your ``.hap`` file contains an "ancestry" extra field and your VCF contains a "POP" format field (as output by ``simgenotype``), you should specify the ``--ancestry`` flag. This will enable us to match the population labels of each haplotype against those in the genotypes output by ``simgenotype``. See :ref:`this section <formats-haplotypes-extrafields-simphenotype>` of the ``.hap`` format spec for more details.

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

Examples
~~~~~~~~
.. code-block:: bash

	haptools transform tests/data/simple.vcf.gz tests/data/simple.hap

Let's try transforming just two samples and let's output to PGEN format:

.. code-block:: bash

	haptools transform -o output.pgen -s NA12878 -s HG01503 tests/data/apoe.vcf.gz tests/data/apoe4.hap

To get progress information, increase the verbosity to "INFO":

.. code-block:: bash

	haptools transform --verbosity INFO -o output.vcf.gz tests/data/example.vcf.gz tests/data/basic.hap.gz

To match haplotypes as well as their ancestral population labels, use the ``--ancestry`` flag:

.. code-block:: bash

	haptools transform --ancestry tests/data/simple-ancestry.vcf tests/data/simple.hap

If your VCF has multi-allelic variants, they must be split into bi-allelic records before you can use ``transform``. After splitting, you should rename the IDs in your file to ensure they remain unique:

.. code-block:: bash

	bcftools norm -m- -Ou input.vcf.gz | \
	bcftools annotate -Ov --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' | \
	haptools transform -o output.vcf.gz /dev/stdin file.hap

..
	To include ancestral population labels in the transformation, use the ``--ancestry`` flag:

	.. code-block:: bash

		haptools transform --ancestry tests/data/example.vcf.gz tests/data/simphenotype.hap

Detailed Usage
~~~~~~~~~~~~~~

.. click:: haptools.__main__:main
   :prog: haptools
   :nested: full
   :commands: transform
