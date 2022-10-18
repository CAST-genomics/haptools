.. _formats-breakpoints:


Breakpoints
===========

Breakpoints files (``.bp`` files) can be used to label a chromosomal strand of each individual by a set of ancestral populations (ex: YRI or CEU).

Each chromosomal strand of each sample has a set of haplotype blocks in the file delimited by a sample header of the form ``{sample}_1`` (for the first chromosomal strand) or ``{sample}_2`` (for the second chromosomal strand).

Each set of haplotype blocks follow a tab-delimited format with the following fields.

.. list-table::
   :widths: 15 15 25
   :header-rows: 1

   * - Name
     - Type
     - Description
   * - pop
     - string
     - The population label of this haplotype block (ex: CEU or YRI)
   * - chrom
     - string
     - The name of the chromosome to which this haplotype block belongs (ex: chr1)
   * - bp
     - integer
     - The base-pair position of the end of the haplotype block (ex: 1001038)
   * - cm
     - float
     - The centimorgan position of the end of the haplotype block (ex: 43.078)

Examples
--------

See `tests/data/outvcf_test.bp <https://github.com/cast-genomics/haptools/blob/main/tests/data/outvcf_test.bp>`_ for an example of a short breakpoint file:

.. include:: ../../tests/data/outvcf_test.bp
   :literal:

See `tests/data/simple.bp <https://github.com/cast-genomics/haptools/blob/main/tests/data/simple.bp>`_ for an example of a longer example:

.. include:: ../../tests/data/simple.bp
   :literal:
