.. _formats-models:


Models
======

This model file format is based on `admix simu's <https://github.com/williamslab/admix-simu>`_. 

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

Simulating 40 samples for 6 generations in this case implies the first generation has population freqs ``Admixed=0, CEU=0.2, YRI=0.8`` and the remaining 2-6 generations have population frequency ``Admixed=1, CEU=0, YRI=0``

Example pulse event model.dat file
----------------------------------

.. code-block::

  40   Admixed   CEU    YRI
  1    0         0.2    0.8
  2    1         0      0
  3    0.5       0.5    0
  4    1         0      0

Simulating 40 samples for 4 generations in this case implies the first generation has population freqs ``Admixed=0, CEU=0.2, YRI=0.8`` the second generation is purely admixed, the third has an event where a pure CEU population is introduced again at freqs ``Admixed=0.5, CEU=0.5, YRI=0`` and finally we end with pure admixture. 

More Example files
------------------
We have generated example model files that simulate current population structures within different populations in America as well as the Caribbean that can be found in the haptools repository here: `haptools/example-files/models/ <https://github.com/CAST-genomics/haptools/tree/main/example-files/models>`_
There are additional example model files that can be found in the haptools repository under `haptools/tests/data/ <https://github.com/CAST-genomics/haptools/tree/main/tests/data>`_