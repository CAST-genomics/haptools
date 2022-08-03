.. _commands-ld:


ld
=========

Compute the pair-wise LD (`Pearson's correlation coefficient <https://numpy.org/doc/stable/reference/generated/numpy.corrcoef.html>`_) between haplotypes (or genotypes) and a single *TARGET* haplotype (or variant).

The ``ld`` command takes as input a set of genotypes in VCF and a list of haplotypes (specified as a :doc:`.hap file </formats/haplotypes>`) and outputs a new :doc:`.hap file </formats/haplotypes>` with the computed LD values in an extra field.

By default, LD is computed with each haplotype in the **.hap** file. To compute LD with the variants in the genotypes file instead, you should use the `--from-gts <#cmdoption-haptools-ld-from-gts>`_ switch. When this mode is enabled, the **.hap** output will be replaced by a tab-delimited text file similar to `PLINK 1.9's .ld file format <https://www.cog-genomics.org/plink/1.9/formats#ld>`_. It will have a header denoting the following columns:

1. CHR - Chromosome code for the variant
2. BP	- Base-pair coordinate of the variant
3. SNP - ID of the variant
4. R - Pearson's correlation coefficient between the variant and the *TARGET*

You may also specify genotypes in PLINK2 PGEN format instead of VCF format. See the documentation for genotypes in :ref:`the format docs <formats-genotypesplink>` for more information.

Usage
~~~~~
.. code-block:: bash

	haptools ld \
	--region TEXT \
	--sample SAMPLE \
	--samples-file FILENAME \
	--id ID \
	--ids-file FILENAME \
	--chunk-size INT \
	--discard-missing \
	--from-gts \
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

Alternatively, we can compute LD between the APOe4 haplotype and all genotypes in the VCF by using the ``--from-gts`` switch.

.. code-block:: bash

	haptools ld --from-gts -o apoe4.ld APOe4 tests/data/apoe.vcf.gz tests/data/apoe4.hap

Detailed Usage
~~~~~~~~~~~~~~

.. click:: haptools.__main__:main
   :prog: haptools
   :nested: full
   :commands: ld
