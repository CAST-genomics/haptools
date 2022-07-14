.. _formats-genotypes:


Genotypes
=========

Genotype files must be specified as VCF or BCF files.

There is also experimental support for PLINK2 PGEN files in some commands. See the documentation for the :class:`GenotypesPLINK` class in :ref:`the API docs <api-data-genotypesplink>` for more information.

.. warning::
	PGEN files are not officially supported yet because our code relies upon the as-yet-unpublished ``pgenlib`` python library. See `issue #16 <https://github.com/gymrek-lab/haptools/pull/16>`_ for current progress on this challenge. In the meantime, you must install the library manually from Github via ``pip``.

	.. code-block:: bash

		pip install git+https://github.com/chrchang/plink-ng.git#subdirectory=2.0/Python
