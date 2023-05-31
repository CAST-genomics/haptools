.. _formats-snplist:


Causal SNPs
===========

You can specify causal SNPs in a tab-delimited ``.snplist`` file. We follow `GCTA's .snplist format <https://yanglab.westlake.edu.cn/software/gcta/#GWASSimulation>`_ for this type of file. It has just two columns:

.. list-table::
   :widths: 15 15 25
   :header-rows: 1

   * - Name
     - Type
     - Description
   * - ID
     - string
     - A unique identifier for this variant in the file (ex: 'rs1234')
   * - BETA
     - float
     - The effect size assigned to this variant (ex: 0.08)

.. note::
  You should not include a header in this file. The file format does not have one.

Examples
--------

Refer to `tests/data/apoe.snplist <https://github.com/cast-genomics/haptools/blob/main/tests/data/apoe.snplist>`_ for an example containing just two SNPs.

.. include:: ../../tests/data/apoe.snplist
  :literal:

Converting to a ``.hap`` file
-----------------------------
The capabilities of the ``.snplist`` format are limited. For example, it does not allow users to specify a causal allele (REF vs ALT) for each SNP. You can use :ref:`the haptools API to upgrade to a .hap file <api-examples-snps2hap>` if needed.
