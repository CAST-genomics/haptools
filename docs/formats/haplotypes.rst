.. _formats-haplotypes:


Haplotypes
==========

This document describes our custom file format specification for haplotypes: the ``.hap`` file.

This is a tab-separated file composed of different types of lines. The first field of each line is a single, uppercase character denoting the type of line. The following line types are supported.

.. list-table::
   :widths: 25 25
   :header-rows: 1

   * - Type
     - Description
   * - #
     - Comment/Header
   * - H
     - Haplotype
   * - V
     - Variant

Each line type (besides ``#``) has a set of mandatory fields described below. Additional "extra" fields can be appended to these to customize the file.

``#`` Comment line
~~~~~~~~~~~~~~~~~~
Comment lines begin with ``#`` and are ignored. They can appear anywhere in a ``.hap`` file.

It is best practice to immediately follow all comment lines with a space. Otherwise, the line may be at risk of being interpreted as part of the header, especially if the file is sorted.

Header line
~~~~~~~~~~~
Header lines begin with ``#`` and must precede all other line types, except for comment lines. Comment lines can be interleaved with the header lines in any order.

There are two types of header lines: those with file metadata and those that declare extra fields.

Header lines themselves can appear in any order. We recommend putting all of the metadata lines before the extra field declarations, as a best practice.

Metadata lines in the header
----------------------------
Metadata lines have the following, tab-separated fields:

1. A header symbol ``#``
2. A unique metadata name
3. Value(s)

It is best practice to include a metadata line declaring the version of the haplotype format that your file uses. Otherwise, your file will be assumed to use the latest version of the specification.

.. code-block::

  #	version	0.1.0

If you are declaring any extra fields (see the next section), then you should include a metadata line that declares the order of the extra fields. For example, if the ``H`` haplotype lines in your file have two extra fields, "ancestry" and "beta", that appear in that order, you would write

.. code-block::

  #	orderH	ancestry	beta

Declaring extra fields in the header
------------------------------------
Any extra fields in the file must be declared in the header. To declare an extra field, create a tab-separated line containing the following fields:

1. A header symbol followed by a line type symbol (ex: ``#H`` or ``#V``)
2. Field name
3. Python format string (ex: 'd' for int, 's' for string, or '.3f' for a float with 3 decimals)
4. Description

Note that the first field must follow the ``#`` symbol immediately (ex: ``#H`` or ``#V``).

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

.. note::
   It is not currently possible to encode haplotypes that span more than one contig.

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
     - string
     - The allele of this variant within the haplotype

Examples
~~~~~~~~
You can find an example of a ``.hap`` file without any extra fields in `tests/data/basic.hap <https://github.com/cast-genomics/haptools/blob/main/tests/data/basic.hap>`_:

.. include:: ../../tests/data/basic.hap
   :literal:

You can find an example with extra fields added within `tests/data/simphenotype.hap <https://github.com/cast-genomics/haptools/blob/main/tests/data/simphenotype.hap>`_:

.. include:: ../../tests/data/simphenotype.hap
   :literal:


Compressing and indexing
~~~~~~~~~~~~~~~~~~~~~~~~
We encourage you to sort, bgzip compress, and index your ``.hap`` file whenever possible. This will reduce both disk usage and the time required to parse the file, but it is entirely optional.

.. code-block:: bash

  LC_ALL=C sort -k1,4 -o file.hap file.hap
  bgzip file.hap
  tabix -s 2 -b 3 -e 4 file.hap.gz

In order to properly index the file, the set of IDs in the haplotype lines must be distinct from the set of chromosome names. This is a best practice in unindexed ``.hap`` files but a requirement for indexed ones. In addition, you must sort on the first field (ie the line type symbol) in addition to the latter three.

Querying an indexed file
~~~~~~~~~~~~~~~~~~~~~~~~
You can query an indexed ``.hap`` file on both the haplotype and variant levels with the following syntax.

.. code-block:: bash

  tabix file.hap.gz REGION

For example, to extract all haplotypes between positions 100 and 200 on chromosome ``chr19``:

.. code-block:: bash

  tabix file.hap.gz chr19:100-200

Or to get all alleles between positions 100 and 200 on the haplotype with ID ``hap1``:

.. code-block:: bash

  tabix file.hap.gz hap1:100-200

Extra fields
~~~~~~~~~~~~
Additional fields can be appended to the ends of the haplotype and variant lines as long as they are declared in the header.

.. _formats-haplotypes-extrafields-transform:

transform
---------
If you would like to simulate an ancestry-based effect, you should run ``transform`` with an *ancestry* extra field declared in your ``.hap`` file.

You can download an example header with an *ancestry* extra field from `tests/data/simphenotype.hap <https://github.com/cast-genomics/haptools/blob/main/tests/data/simphenotype.hap>`_

.. code-block:: bash

  curl https://raw.githubusercontent.com/cast-genomics/haptools/main/tests/data/simphenotype.hap 2>/dev/null | head -n4

``H`` Haplotype
+++++++++++++++

.. list-table::
   :widths: 25 25 25 50
   :header-rows: 1

   * - Column
     - Field
     - Type
     - Description
   * - 5
     - Local Ancestry
     - string
     - A population code denoting this haplotype's ancestral origins

``V`` Variant
+++++++++++++
No extra fields are required here.

.. _formats-haplotypes-extrafields-simphenotype:

simphenotype
------------
The *beta* extra field should be declared for your ``.hap`` file to be compatible with the ``simphenotype`` subcommand.

You can download an example header with a *beta* extra field from `tests/data/simphenotype.hap <https://github.com/cast-genomics/haptools/blob/main/tests/data/simphenotype.hap>`_

.. code-block:: bash

  curl https://raw.githubusercontent.com/cast-genomics/haptools/main/tests/data/simphenotype.hap 2>/dev/null | head -n4


..
  _TODO: figure out how to tab this code block so that the tabs get copied when someone copies from it


``H`` Haplotype
+++++++++++++++

.. list-table::
   :widths: 25 25 25 50
   :header-rows: 1

   * - Column
     - Field
     - Type
     - Description
   * - 5
     - Effect Size
     - float
     - The effect size of this haplotype; for use in ``simphenotype``

``V`` Variant
+++++++++++++
No extra fields are required here.

Changelog
~~~~~~~~~
v0.1.0
------
Updates to the header lines in the specification. See `PR #80 <https://github.com/cast-genomics/haptools/pull/80>`_.

We've created a new type of metadata line for specifying the "order" of the extra fields in each line.
In the absence of this metadata line, the extra fields will be assumed to appear in the order of the extra-field declarations in the header. Unfortunately, sorting can change that. By specifying the order of the extra fields up-front, you can ensure that the file will be parsed the same regardless of whether it is sorted.

In addition, we now allow you to have additional extra fields besides the ones that are used by the specific tool you are using. For example, the ``transform`` subcommand used to complain if it found any extra fields in your ``.hap`` file. But now, it will gracefully ignore those extra fields and load only the fields that it might need.

If your ``.hap`` file does not have any extra fields, you can safely bump the version number without changing the rest of your file.

v0.0.1
------
Initialized the spec! See `PR #43 <https://github.com/cast-genomics/haptools/pull/43>`_.
