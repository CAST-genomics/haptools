.. _commands-index:


index
=========

Index a list of haplotypes.

The ``index`` command takes as input a list of haplotypes (specified as a :doc:`.hap file </formats/haplotypes>`) and outputs the list as a .gz and a .tbi file. It also sorts the file in the process.

Usage
~~~~~
.. code-block:: bash

	haptools index \
	--output PATH \
	--verbosity [CRITICAL|ERROR|WARNING|INFO|DEBUG|NOTSET] \
	HAPLOTYPES

Example
~~~~~~~~
.. code-block:: bash

	haptools index tests/data/basic.hap

You may also specify a custom output path for the compressed file to be written to.

.. code-block:: bash

	haptools index --output tests/data/sorted.basic.hap.gz tests/data/basic.hap

Use the ``--no-sort`` flag to skip the sorting step if your file is already sorted.

.. code-block:: bash

	haptools index --no-sort --output tests/data/basic.hap.gz tests/data/basic.hap


Detailed Usage
~~~~~~~~~~~~~~

.. click:: haptools.__main__:main
   :prog: haptools
   :show-nested:
   :commands: index
