.. _commands-transform:


transform
=========

Transform a set of genotypes via a list of haplotypes. Create a new VCF containing haplotypes instead of variants.

The ``transform`` command takes as input a set of genotypes in VCF and a list of haplotypes (specified as a :doc:`.hap file </formats/haplotypes>`) and outputs a set of haplotype "genotypes" in VCF.

You may also specify genotypes in PLINK2 PGEN format. See the documentation for genotypes in :ref:`the format docs <formats-genotypesplink>` for more information.

Usage
~~~~~
.. code-block:: bash

	haptools transform \
	--region TEXT \
	--sample SAMPLE \
	--samples-file FILENAME \
	--hap-id ID \
	--chunk-size INT \
	--discard-missing \
	--output PATH \
	--verbosity [CRITICAL|ERROR|WARNING|INFO|DEBUG|NOTSET] \
	GENOTYPES HAPLOTYPES

Examples
~~~~~~~~
.. code-block:: bash

	haptools transform tests/data/example.vcf.gz tests/data/basic.hap.gz | less

.. code-block:: bash

	haptools transform -o output.vcf.gz -s NA12878 tests/data/apoe.vcf.gz tests/data/apoe4.hap

Detailed Usage
~~~~~~~~~~~~~~

.. click:: haptools.__main__:main
   :prog: haptools
   :show-nested:
   :commands: transform
