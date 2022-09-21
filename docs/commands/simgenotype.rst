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
  --region CHR:START-END \
  --invcf REFVCF \
  --sample_info SAMPLEINFOFILE \
  --out OUTPREFIX
  
Detailed information about each option, and example commands using publicly available files, are shown below.

Parameter Descriptions
~~~~~~~~~~~~~~~~~~~~~~
`--from-gts <#cmdoption-haptools-ld-from-gts>`_
* ``--model`` - Parameters for simulating admixture across generations including sample size, population fractions, and number of generations.
* ``--mapdir`` - Directory containing all .map files with the `structure <https://www.cog-genomics.org/plink/1.9/formats#map>`_ where the third position is in centiMorgans
* ``--chroms`` - List of chromosomes to be simulated. The map file directory must contain the "chr<CHR>" where <CHR> is the chromosome identifier eg. 1,2,...,X
* ``--region`` - Limit the simulation to a region within a single chromosome. Overwrites chroms with the chrom listed in this region. eg 1:1-10000 [Optional]
* ``--invcf`` - Input VCF file used to simulate specifiic haplotypes for resulting samples
* ``--sample_info`` - File used to map samples in ``REFVCF`` to populations found in ``MODELFILE``
* ``--out`` - Output prefix of the structure ``/path/to/output`` which results in the vcf file ``output.vcf.gz`` and breakpoints file ``output.bp``

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
