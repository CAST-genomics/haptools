.. _formats-models:


Models
======

Structure of a model.dat file
-----------------------------

* ``num_samples`` - Total number of samples to be output by the simulator (``num_samples*2`` haplotypes)  
* ``num_generations`` - Number of generations to simulate admixture, must be > 0  
* ``*_freq`` - Frequency of populations to be present in the simulated samples

.. code-block::

  {num_samples} Admixed Pop_label1 Pop_label2 ... Pop_labeln
  {num_generations} {admixed_freq} {pop_label1_freq} {pop_label2_freq} ... {pop_labeln_freq}

Example model.dat file
----------------------

.. code-block::

  40   Admixed   CEU    YRI
  6    0         0.2    0.8

Simulating 6 generations in this case implies the first generation has population freqs ``Admixed=0, CEU=0.2, YRI=0.8`` and the remaining 2-6 generations have population frequency ``Admixed=1, CEU=0, YRI=0``
