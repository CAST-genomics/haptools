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

Example
~~~~~~~~
.. code-block:: bash

	haptools index tests/data/basic.hap


Detailed Usage
~~~~~~~~~~~~~~

.. click:: haptools.__main__:main
   :prog: haptools
   :show-nested:
   :commands: index
