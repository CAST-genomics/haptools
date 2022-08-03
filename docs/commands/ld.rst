.. _commands-ld:


ld
=========

Compute the pair-wise LD (Pearson's correlation) between haplotypes and a *TARGET* variant or haplotype.

The ``ld`` command takes as input a set of genotypes in VCF and a list of haplotypes (specified as a :doc:`.hap file </formats/haplotypes>`) and outputs a new :doc:`.hap file </formats/haplotypes>` with the computed LD values in an extra field.

You may also specify genotypes in PLINK2 PGEN format. See the documentation for genotypes in :ref:`the format docs <formats-genotypesplink>` for more information.

Usage
~~~~~
.. code-block:: bash

	haptools ld \
	--region TEXT \
	--sample SAMPLE \
	--samples-file FILENAME \
	--hap-id ID \
	--chunk-size INT \
	--discard-missing \
	--output PATH \
	--verbosity [CRITICAL|ERROR|WARNING|INFO|DEBUG|NOTSET] \
	TARGET GENOTYPES HAPLOTYPES

Examples
~~~~~~~~
*TARGET* can either be a haplotype or a variant.

For example, let's compute LD with the haplotype 'chr21.q.3365*1'.

.. code-block:: bash

	haptools ld 'chr21.q.3365*1' tests/data/example.vcf.gz tests/data/basic.hap.gz | less

Or, let's compute LD with the variant 'rs429358'.

.. code-block:: bash

	haptools ld -o apoe4_ld.hap rs429358 tests/data/apoe.vcf.gz tests/data/apoe4.hap

Detailed Usage
~~~~~~~~~~~~~~

.. click:: haptools.__main__:main
   :prog: haptools
   :nested: full
   :commands: ld
