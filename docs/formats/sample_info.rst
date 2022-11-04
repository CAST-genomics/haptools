.. _formats-sample_info:


Sample Info
===========

A *samples info* file maps samples in a reference to their population listed in a :doc:`model file </formats/models>`. This file is used by the :doc:`simgenotype </commands/simgenotype>` command.

1000 Genomes sample_info file format
------------------------------------
You can download a *samples info* file compatible with the 1000G reference by executing the following.

.. code-block:: bash

	wget https://raw.githubusercontent.com/CAST-genomics/haptools/main/example-files/1000genomes_sampleinfo.tsv

If you'd like to compute this mapping file yourself, execute the following:

.. code-block:: bash

	cut -f 1,4 "igsr-1000 genomes on grch38.tsv" | \
	sed '1d' | \
	sed -e 's/ /\t/g' > 1000genomes_sampleinfo.tsv

Examples
--------------------------

.. code-block::

	HG00372	FIN
	HG00132	GBR
	HG00237	GBR
	HG00404	CHS

See `example-files/1000genomes_sampleinfo.tsv <https://github.com/cast-genomics/haptools/blob/main/example-files/1000genomes_sampleinfo.tsv>`_ for an example of the 1000genomes 
GRCh38 samples mapped to their subpopulations.

.. include:: ../../example-files/1000genomes_sampleinfo.tsv
    :literal:
    :start-line: 13
    :end-line: 25
