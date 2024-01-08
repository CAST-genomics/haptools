.. _commands-clump:


clump
=====

Clump a set of variants specified in a :doc:`.linear file </formats/linear>`.

The ``clump`` command creates a clump file joining SNPs or STRs in LD with one another.

Usage
~~~~~
.. code-block:: bash

  haptools clump \
  --verbosity [CRITICAL|ERROR|WARNING|INFO|DEBUG|NOTSET] \
  --summstats-snps PATH \
  --gts-snps PATH \
  --summstats-strs PATH \
  --gts-strs PATH \
  --clump-field TEXT \
  --clump-id-field TEXT \
  --clump-chrom-field TEXT \
  --clump-pos-field TEXT \
  --clump-p1 FLOAT \
  --clump-p2 FLOAT \
  --clump-r2 FLOAT \
  --clump-kb FLOAT \
  --ld [Exact|Pearson] \
  --out PATH

Examples
~~~~~~~~
.. code-block:: bash

  haptools clump \
    --summstats-snps tests/data/test_snpstats.linear \
    --gts-snps tests/data/simple.vcf \
    --clump-id-field ID \
    --clump-chrom-field CHROM \
    --clump-pos-field POS \
    --out test_snps.clump

You can use ``--ld [Exact|Pearson]`` to indicate which type of LD calculation you'd like to perform. ``Exact`` utilizes an exact cubic solution adopted from `CubeX <https://github.com/t0mrg/cubex>`_ whereas ``Pearson`` utilizes a Pearson R calculation. Note ``Exact`` only works on SNPs and not any other variant type eg. STRs. The default value is ``Pearson``.

.. code-block:: bash

  haptools clump \
    --summstats-snps tests/data/test_snpstats.linear \
    --gts-snps tests/data/simple.vcf \
    --clump-id-field ID \
    --clump-chrom-field CHROM \
    --clump-pos-field POS \
    --ld Exact \
    --out test_snps.clump

You can modify thresholds and values used in the clumping process. ``--clump-p1`` is the largest value of a p-value to consider being an index variant for a clump. ``--clump-p2`` dictates the maximum p-value any variant can have to be considered when clumping. ``--clump-r2`` is the R squared threshold where being greater than this value implies the candidate variant is in LD with the index variant. ``--clump-kb`` is the maximum distance upstream or downstream from the index variant to collect candidate variants for LD comparison. For example, ``--clump-kb 100`` implies all variants 100 Kb upstream and 100 Kb downstream from the variant will be considered.

.. code-block:: bash

  haptools clump \
    --summstats-snps tests/data/test_snpstats.linear \
    --gts-snps tests/data/simple.vcf \
    --clump-id-field ID \
    --clump-chrom-field CHROM \
    --clump-pos-field POS \
    --clump-p1 0.001 \
    --clump-p2 0.05 \
    --clump-r2 0.7 \
    --clump-kb 200.5 \
    --out test_snps.clump

You can also input STRs when calculating clumps. They can be used together with SNPs or alone.

.. code-block:: bash

  haptools clump \
    --summstats-strs tests/data/test_strstats.linear \
    --gts-strs tests/data/simple_tr.vcf \
    --summstats-snps tests/data/test_snpstats.linear \
    --gts-snps tests/data/simple.vcf \
    --clump-id-field ID \
    --clump-chrom-field CHROM \
    --clump-pos-field POS \
    --ld Exact \
    --out test_snps.clump

.. code-block:: bash

  haptools clump \
    --summstats-strs tests/data/test_strstats.linear \
    --gts-strs tests/data/simple_tr.vcf \
    --clump-id-field ID \
    --clump-chrom-field CHROM \
    --clump-pos-field POS \
    --ld Exact \
    --out test_snps.clump

All files used in these examples are described :doc:`here </project_info/example_files>`.


Detailed Usage
~~~~~~~~~~~~~~

.. click:: haptools.__main__:main
   :prog: haptools
   :show-nested:
   :commands: clump