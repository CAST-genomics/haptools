.. _formats-linear:


Summary Statistics
==================

Linear files (``.linear`` files) store summary statistics from a linear model. We follow the `PLINK2 .glm.linear file format <https://www.cog-genomics.org/plink/2.0/formats#glm_linear>`_ and use the following columns.

.. list-table::
   :widths: 15 15 25
   :header-rows: 1

   * - Name
     - Type
     - Description
   * - CHROM
     - string
     - The name of the chromosome to which this SNP belongs (ex: 1)
   * - POS
     - integer
     - The position of this SNP on the chromosome (ex: 10114)
   * - ID
     - string
     - A unique identifier for this SNP in the file (ex: 'rs1234')
   * - P
     - float
     - The p-value assigned to this SNP via association testing (ex: 43.078)

Examples
--------

See `tests/data/test_snpstats.linear <https://github.com/cast-genomics/haptools/blob/main/tests/data/test_snpstats.linear>`_ for an example of a short ``.linear`` file:

.. include:: ../../tests/data/test_snpstats.linear
  :literal:
