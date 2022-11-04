.. _formats-breakpoints:


Breakpoints
===========

Breakpoints files (``.bp`` files) store your samples' local ancestry labels. Each line in the file denotes the ancestral population (ex: YRI or CEU) of a portion of a chromosomal strand (or *haplotype block*) of an individual.

The set of haplotype blocks for an individual are delimited by a sample header of the form ``{sample}_1`` (for the first chromosomal strand) or ``{sample}_2`` (for the second chromosomal strand). Blocks from ``{sample}_1`` must be directly followed by blocks from ``{sample}_2``.

Each set of haplotype blocks follows a tab-delimited format with the following fields. Lines within a sample's set of blocks must be sorted according to ``chrom``, ``bp``, and ``cm`` - in that order.

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

See `tests/data/simple.bp <https://github.com/cast-genomics/haptools/blob/main/tests/data/simple.bp>`_ for a longer example:

.. include:: ../../tests/data/simple.bp
  :literal:
