.. _commands-ld:


ld
===

Compute the pair-wise LD (`Pearson's correlation coefficient <https://numpy.org/doc/stable/reference/generated/numpy.corrcoef.html>`_) between haplotypes (or genotypes) and a single *TARGET* haplotype (or variant).

The ``ld`` command takes as input a set of genotypes in VCF and a list of haplotypes (specified as a :doc:`.hap file </formats/haplotypes>`) and outputs a new :doc:`.hap file </formats/haplotypes>` with the computed LD values in an extra field.

By default, LD is computed with each haplotype in the ``.hap`` file. To compute LD with the variants in the genotypes file instead, you should use the `--from-gts <#cmdoption-haptools-ld-from-gts>`_ switch. When this mode is enabled, the ``.hap`` output will be replaced by an :doc:`.ld file </formats/ld>`.

.. note::
	Repeats are not currently supported by the ``ld`` command. Any repeats in your ``.hap`` file will be ignored.

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

Alternatively, we can compute LD between the APOe4 haplotype and all genotypes in the VCF by using the ``--from-gts`` switch. Note that we should use a different extension for the output file now.

.. code-block:: bash

	haptools ld --from-gts -o apoe4.ld APOe4 tests/data/apoe.vcf.gz tests/data/apoe4.hap

You can select a subset of variants (or haplotypes) using the ``--id`` parameter multiple times (or the ``--ids-file`` parameter).

.. code-block:: bash

	haptools ld --from-gts -i rs543363163 -i rs7412 APOe4 tests/data/apoe.vcf.gz tests/data/apoe4.hap

All files used in these examples are described :doc:`here </project_info/example_files>`.

Detailed Usage
~~~~~~~~~~~~~~

.. click:: haptools.__main__:main
   :prog: haptools
   :nested: full
   :commands: ld
