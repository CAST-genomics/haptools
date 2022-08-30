.. _commands-index:


index
=========

Index a list of haplotypes.

The ``index`` command takes as input a list of haplotypes (specified as a :doc:`.hap file </formats/haplotypes>`) and outputs the list as a .gz and a .tbi file.

Usage
~~~~~
.. code-block:: bash

	haptools index \
	--output PATH \
	--verbosity [CRITICAL|ERROR|WARNING|INFO|DEBUG|NOTSET] \
	HAPLOTYPES

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
   :commands: index
