.. _commands-simgenotype:


simgenotype
===========

Takes as input a reference set of haplotypes in VCF format and a user-specified admixture model.

Outputs a VCF file with simulated genotype information for admixed genotypes, as well as a breakpoints file that can be used for visualization.

Basic Usage
~~~~~~~~~~~
.. code-block:: bash

  haptools simgenotype \
  --model MODELFILE \
  --mapdir GENETICMAPDIR \
  --chroms LIST,OF,CHROMS \
  --invcf REFVCF \
  --sample_info SAMPLEINFOFILE \
  --out OUTPREFIX

Detailed information about each option, and example commands using publicly available files, are shown below.

Parameter Descriptions
~~~~~~~~~~~~~~~~~~~~~~
* ``--invcf`` - Input VCF file used to simulate specifiic haplotypes for resulting samples
* ``--sample_info`` - File used to map samples in ``REFVCF`` to populations found in ``MODELFILE``
* ``--model`` - Parameters for simulating admixture across generations
* ``--map`` - .map file used to determine recombination events during the simulation
* ``--out`` - Output prefix of the structure ``/path/to/output`` which results in the vcf file ``output.vcf.gz`` and breakpoints file ``output.bp``

File Formats
~~~~~~~~~~~~
Model Format
------------

**Structure of model.dat file**

* ``num_samples`` - Total number of samples to be output by the simulator (``num_samples*2`` haplotypes)  
* ``num_generations`` - Number of generations to simulate admixture, must be > 0  
* ``*_freq`` - Frequency of populations to be present in the simulated samples

.. code-block::

  {num_samples} Admixed Pop_label1 Pop_label2 ... Pop_labeln
  {num_generations} {admixed_freq} {pop_label1_freq} {pop_label2_freq} ... {pop_labeln_freq}

**Example model.dat file**

.. code-block::

  40   Admixed   CEU    YRI
  6    0         0.2    0.8

Simulating 6 generations in this case implies the first generation has population freqs ``Admixed=0, CEU=0.2, YRI=0.8`` and the remaining 2-6 generations have population frequency ``Admixed=1, CEU=0, YRI=0``

Map Format
----------

* ``chr`` - chromosome of coordinate (1-22, X)  
* ``var`` - variant identifier   
* ``pos cM`` - Position in centimorgans   
* ``pos bp`` - Base pair coordinate  

.. code-block::python

  {chr}\t{var}\t{pos cM}\t{pos bp}

Beagle Genetic Maps used in simulation (GRCh38): http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/


Outfile Format
--------------

* ``Sample Header`` - Name of sample following the structure ``Sample_{number}_{hap}`` eg. ``Sample_10_1`` for sample number 10 haplotype 1  
* ``pop`` - Population label corresponding to the index of the population in the dat file so in the example above CEU = 1, YRI = 2  
* ``chr`` - chromosome (1-22, X)  

.. code-block::

  Sample Header
  {pop}\t{chr}\t{pos bp}
  ...
  Sample Header 2
  ...

Examples
~~~~~~~~

.. code-block:: bash

  haptools simgenotype \
  --model tests/data/outvcf_gen.dat \
  --mapdir tests/data/map/ \
  --chroms 1,2 \
  --invcf tests/data/outvcf_test.vcf \
  --sample_info tests/data/outvcf_info.tab \
  --out tests/data/example_simgenotype


Detailed Usage
~~~~~~~~~~~~~~

.. click:: haptools.__main__:main
   :prog: haptools
   :nested: full
   :commands: simgenotype
