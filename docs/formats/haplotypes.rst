.. _formats-haplotypes:


haplotypes
==========

This document describes our custom file format specification for haplotypes: the ``.hap`` file.

This is a tab-separated file composed of different types of lines. The first field of each line is a single, uppercase character denoting the type of line.

.. list-table::
   :widths: 25 25
   :header-rows: 1

   * - Type
     - Description
   * - #
     - Comment
   * - H
     - Haplotype
   * - V
     - Variant

``#`` Comment line
~~~~~~~~~~~~~~~~~~
Comment lines begin with ``#`` and are ignored. Extra fields are also declared here. The following extra fields should be declared for your ``.hap`` file to be compatible with ``haptools``:

- TODO: add extra fields for LA and beta

``H`` Haplotype
~~~~~~~~~~~~~~~
Haplotypes contain the following attributes:

.. list-table::
   :widths: 25 25 25 50
   :header-rows: 1

   * - Column
     - Field
     - Type
     - Description
   * - 1
     - Chromosome
     - string
     - The contig that this haplotype belongs on
   * - 2
     - Start Position
     - int
     - The start position of this haplotype on this contig
   * - 3
     - End Position
     - int
     - The end position of this haplotype on this contig
   * - 4
     - Haplotype ID
     - string
     - Uniquely identifies a haplotype
   * - 5
     - Local Ancestry
     - string
     - A population code denoting this haplotype's ancestral origins
   * - 6
     - Effect Size
     - float
     - The effect size of this haplotype; for use in ``simphenotype``

``V`` Variant
~~~~~~~~~~~~~
Each variant line belongs to a particular haplotype. These lines contain the following attributes:

.. list-table::
   :widths: 25 25 25 50
   :header-rows: 1

   * - Column
     - Field
     - Type
     - Description
   * - 1
     - Haplotype ID
     - string
     - Identifies the haplotype to which this variant belongs
   * - 2
     - Start Position
     - int
     - The start position of this variant on its contig
   * - 3
     - End Position
     - int
     - The end position of this variant on its contig

       Usually the same as the Start Position
   * - 4
     - Variant ID
     - string
     - The unique ID for this variant, as defined in the genotypes file
   * - 5
     - Allele
     - int
     - The allele of this variant within the haplotype
