.. _formats-ld:


Linkage disequilibrium
======================

Linkage disequilibrium for SNPs can be specified in a tab-delimited format similar to `PLINK 1.9's .ld file format <https://www.cog-genomics.org/plink/1.9/formats#ld>`_. Our format contains only a subset of the full set of columns from PLINK's:

.. list-table::
   :widths: 15 15 25
   :header-rows: 1

   * - Name
     - Type
     - Description
   * - CHR
     - string
     - The name of the chromosome to which this SNP belongs (ex: 1)
   * - BP
     - integer
     - The position of this SNP on the chromosome (ex: 10114)
   * - SNP
     - string
     - A unique identifier for this SNP in the file (ex: 'rs1234')
   * - R
     - float
     - Pearson's correlation coefficient between the variant and some *TARGET* (ex: 0.42)


Examples
--------

.. code-block:: text

  CHR	BP	SNP	R
  19	45411941	rs429358	0.999
  19	45411947	rs11542041	0.027
  19	45411962	rs573658040	-0.012
  19	45411965	rs543363163	-0.012
  19	45412006	rs563140413	-0.012
  19	45412007	rs531939919	-0.012
  19	45412040	rs769455	0.006
  19	45412079	rs7412	-0.098
