.. _formats-genotypes:


Genotypes
=========

.. figure:: https://drive.google.com/uc?export=view&id=1_JARKJQ0LX-DzL0XsHW1aiQgLCOJ1ZvC

	The time required to load various genotype file formats.

VCF/BCF
-------

Genotype files must be specified as VCF or BCF files. They can be bgzip-compressed.

.. _formats-genotypesplink:

PLINK2 PGEN
-----------

There is also experimental support for `PLINK2 PGEN <https://github.com/chrchang/plink-ng/blob/master/pgen_spec/pgen_spec.pdf>`_ files in some commands. These files can be loaded and created much more quickly than VCFs, so we highly recommend using them if you're working with large datasets. See the documentation for the :class:`GenotypesPLINK` class in :ref:`the API docs <api-data-genotypesplink>` for more information.

If you run out memory when using PGEN files, consider reading/writing variants from the file in chunks via the ``--chunk-size`` parameter.

.. note::
	PLINK2 support depends on the ``Pgenlib`` python library. This can be installed automatically with ``haptools`` if you specify the "files" extra requirements during installation.

	.. code-block:: bash

		pip install haptools[files]

.. warning::
	At the moment, only biallelic SNPs can be encoded in PGEN files because of limitations in the ``Pgenlib`` python library. It doesn't properly support multiallelic variants yet (`source <https://github.com/chrchang/plink-ng/blob/c4b8d4361de74c58f0cc11361062eca4f34210d3/2.0/Python/python_api.txt#L88-L89>`_). To ensure your PGEN files only contain SNPs, we recommend use the following command to convert from VCF to PGEN.

	.. code-block:: bash

		plink2 --snps-only 'just-acgt' --max-alleles 2 --vcf input.vcf --make-pgen --out output
