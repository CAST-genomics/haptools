.. _commands-transform:


transform
=========

Transform a set of genotypes via a list of haplotypes. Create a new VCF containing haplotypes instead of variants.

The ``transform`` command takes as input a set of genotypes in VCF and a list of haplotypes (specified as a :doc:`.hap file </formats/haplotypes>`) and outputs a set of haplotype "genotypes" in VCF.

You may also specify genotypes in PLIN2 PGEN format. See the documentation for the :class:`GenotypesPLINK` class in :ref:`the API docs <api-data-genotypesplink>` for more information.

.. warning::
	PGEN files are not officially supported yet because our code relies upon the as-yet-unpublished ``pgenlib`` python library. See `issue #16 <https://github.com/gymrek-lab/haptools/pull/16>`_ for current progress on this challenge. In the meantime, you must install the library manually from Github via ``pip``.

	.. code-block:: bash

		pip install git+https://github.com/chrchang/plink-ng.git#subdirectory=2.0/Python

Usage
~~~~~
.. code-block:: bash

	haptools transform \
	--region TEXT \
	--sample SAMPLE \
	--samples-file FILENAME \
	--output PATH \
	--verbosity [CRITICAL|ERROR|WARNING|INFO|DEBUG|NOTSET] \
	GENOTYPES HAPLOTYPES

Examples
~~~~~~~~
.. code-block:: bash

	haptools transform tests/data/example.vcf.gz tests/data/basic.hap.gz | less

Detailed Usage
~~~~~~~~~~~~~~

.. click:: haptools.__main__:main
   :prog: haptools
   :show-nested:
   :commands: transform
