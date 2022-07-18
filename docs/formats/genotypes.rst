.. _formats-genotypes:


Genotypes
=========

Genotype files must be specified as VCF or BCF files.

.. _formats-genotypesplink:

There is also experimental support for `PLINK2 PGEN <https://github.com/chrchang/plink-ng/blob/master/pgen_spec/pgen_spec.pdf>`_ files in some commands. These files can be loaded much more quickly than VCFs, so we highly recommend using them if you're working with large datasets. See the documentation for the :class:`GenotypesPLINK` class in :ref:`the API docs <api-data-genotypesplink>` for more information.

.. figure:: https://drive.google.com/uc?export=view&id=1_JARKJQ0LX-DzL0XsHW1aiQgLCOJ1ZvC

	The time required to load various genotype file formats.

.. warning::
	PGEN files are not officially supported yet because our code relies upon the as-yet-unpublished ``pgenlib`` python library. See `issue #16 <https://github.com/gymrek-lab/haptools/pull/16>`_ for current progress on this challenge. In the meantime, you must install the library manually from Github via ``pip``.

	.. code-block:: bash

		pip install git+https://github.com/chrchang/plink-ng.git#subdirectory=2.0/Python
