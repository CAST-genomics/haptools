.. _commands-simgenotype:


simgenotype
===========

Takes as input a reference set of haplotypes in VCF format and a user-specified admixture model.

Outputs a VCF file with simulated genotype information for admixed genotypes, as well as a breakpoints file that can be used for visualization. For example, you could simulate a 50/50 mixture of CEU and YRI for 10 generations. Other more complex models such as involving pulse events of new populations can also be simulated.

Basic Usage
~~~~~~~~~~~
.. code-block:: bash

  haptools simgenotype \
  --model MODELFILE \
  --mapdir GENETICMAPDIR \
  --chroms LIST,OF,CHROMS \
  --region CHR:START-END \
  --ref_vcf REFVCF \
  --sample_info SAMPLEINFOFILE \
  --pop_field \
  --out /PATH/TO/OUTPUT.VCF.GZ
  
Detailed information about each option, and example commands using publicly available files, are shown below.

Parameter Descriptions
~~~~~~~~~~~~~~~~~~~~~~
* ``--model`` - Parameters for simulating admixture across generations including sample size, population fractions, and number of generations.
* ``--mapdir`` - Directory containing all .map files with this `structure <https://www.cog-genomics.org/plink/1.9/formats#map>`_ where the third position is in centiMorgans
* ``--out`` - Full output path to file of the structure ``/path/to/output.(vcf|bcf|vcf.gz|pgen)`` which if ``vcf.gz`` is chosen outputs ``/path/to/output.vcf.gz`` and breakpoints file ``/path/to/output.bp``
* ``--chroms`` - List of chromosomes to be simulated. The map file directory must contain the "chr<CHR>" where <CHR> is the chromosome identifier eg. 1,2,...,X
* ``--seed`` - Seed for randomized calculations during simulation of breakpoints. [Optional]
* ``--popsize`` - Population size for each generaetion that is sampled from to create our simulated samples. Default = max(10000, 10*samples) [Optional]
* ``--ref_vcf`` - Input VCF or PGEN file used to simulate specifiic haplotypes for resulting samples
* ``--sample_info`` - File used to map samples in ``REFVCF`` to populations found in ``MODELFILE``
* ``--region`` - Limit the simulation to a region within a single chromosome. Overwrites chroms with the chrom listed in this region. eg 1:1-10000 [Optional]
* ``--pop_field`` - Flag for ouputting population field in VCF output. Note this flag does not work when your output is in PGEN format. [Optional]
* ``--sample_field`` - Flag for ouputting sample field in VCF output. Note this flag does not work when your output is in PGEN format. Should only be used for debugging. [Optional]
* ``--no_replacement`` - Flag for deteremining during the VCF generation process whether we grab samples' haplotypes with or without replacement from the reference VCF file. Default = False (With replacement) [Optional]
* ``--verbosity`` - What level of output the logger should print to stdout. Please see `logging levels <https://docs.python.org/3/library/logging.html>`_ for output levels. Default = INFO [Optional]
* ``--only_breakpoint`` - Flag which when provided only outputs the breakpoint file. Note you will not need to provide a ``--ref_vcf`` or ``--sample_info`` file and can instead put NA. eg.  ``--ref_vcf NA`` and ``--sample_info NA`` [Optional]

File Formats
~~~~~~~~~~~~
* :doc:`sampleinfo format </formats/sample_info>`
* :doc:`model format </formats/models>`
* :doc:`map format </formats/maps>`
* :doc:`breakpoint format </formats/breakpoints>`

Examples
~~~~~~~~

.. code-block:: bash

  haptools simgenotype \
  --model tests/data/outvcf_gen.dat \
  --mapdir tests/data/map/ \
  --region 1:1-83000 \
  --ref_vcf tests/data/outvcf_test.vcf.gz \
  --sample_info tests/data/outvcf_info.tab \
  --pop_field \
  --out tests/data/example_simgenotype.vcf

If speed is important, it's generally faster to use PGEN files than VCFs.

.. code-block:: bash

  haptools simgenotype \
  --model tests/data/outvcf_gen.dat \
  --mapdir tests/data/map/ \
  --region 1:1-83000 \
  --ref_vcf tests/data/outvcf_test.pgen \
  --sample_info tests/data/outvcf_info.tab \
  --pop_field \
  --out tests/data/example_simgenotype.pgen

.. warning::
  Writing PGEN files will require more memory than writing VCFs. The memory will depend on the number of simulated samples and variants.
  You can reduce the memory required for this step by writing the variants in chunks. Just specify a ``--chunk-size`` value.

All files used in these examples are described :doc:`here </project_info/example_files>`.


Detailed Usage
~~~~~~~~~~~~~~

.. click:: haptools.__main__:main
   :prog: haptools
   :nested: full
   :commands: simgenotype
