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
