.. _commands-transform:


transform
=========

Transform a set of genotypes via a list of haplotypes. Create a new VCF containing haplotypes instead of variants.

The ``transform`` command takes as input a set of genotypes in VCF and a list of haplotypes (specified as a :doc:`.hap file </formats/haplotypes>` without any extra fields) and outputs a set of haplotype "genotypes" in VCF.

You may also specify genotypes in PLINK2 PGEN format. Just use the appropriate ".pgen" file extension in the input and/or output. See the documentation for genotypes in :ref:`the format docs <formats-genotypesplink>` for more information.

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
	--output PATH \
	--verbosity [CRITICAL|ERROR|WARNING|INFO|DEBUG|NOTSET] \
	GENOTYPES HAPLOTYPES

Examples
~~~~~~~~
.. code-block:: bash

	haptools transform tests/data/example.vcf.gz tests/data/basic.hap.gz | less

Let's try transforming just two samples and let's output to PGEN format:

.. code-block:: bash

	haptools transform -o output.pgen -s NA12878 -s HG01503 tests/data/apoe.vcf.gz tests/data/apoe4.hap

To get progress information, increase the verbosity to "INFO":

.. code-block:: bash

	haptools transform --verbosity INFO -o output.vcf.gz tests/data/apoe.vcf.gz tests/data/apoe4.hap

Detailed Usage
~~~~~~~~~~~~~~

.. click:: haptools.__main__:main
   :prog: haptools
   :nested: full
   :commands: transform
