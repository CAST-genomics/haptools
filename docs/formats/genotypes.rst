.. _formats-genotypes:


Genotypes
=========

.. figure:: https://drive.google.com/uc?export=view&id=1_JARKJQ0LX-DzL0XsHW1aiQgLCOJ1ZvC

	The time required to load various genotype file formats.

VCF/BCF
~~~~~~~

Genotype files must be specified as VCF or BCF files. They can be bgzip-compressed.

.. _formats-genotypesplink:

PLINK2 PGEN
~~~~~~~~~~~

There is also experimental support for `PLINK2 PGEN <https://github.com/chrchang/plink-ng/blob/master/pgen_spec/pgen_spec.pdf>`_ files (accomponied by PVAR and PSAM files) in some commands. These files can be loaded and created much more quickly than VCFs, so we highly recommend using them if you're working with large datasets. See the documentation for the :class:`GenotypesPLINK` class in :ref:`the API docs <api-data-genotypesplink>` for more information.

If you run out memory when using PGEN files, consider reading/writing variants from the file in chunks via the ``--chunk-size`` parameter.

Converting from VCF to PGEN
---------------------------
To convert a VCF containing only biallelic SNPs to PGEN, use the following command.

.. code-block:: bash

	plink2 --snps-only 'just-acgt' --max-alleles 2 --vcf input.vcf --make-pgen --out output

To convert a VCF containing tandem repeats to PGEN, use this command, instead.

.. code-block:: bash

	plink2 --vcf-half-call m --make-pgen 'pvar-cols=vcfheader,qual,filter,info' --vcf input.vcf --make-pgen --out output
