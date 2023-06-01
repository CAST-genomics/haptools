.. _commands-simphenotype:


simphenotype
============

Simulates a complex trait, taking into account haplotype- or local-ancestry- specific effects as well as traditional variant-level effects. The user denotes causal variants or haplotypes by specifying them in a :doc:`.snplist file </formats/snplist>` or :doc:`.hap file </formats/haplotypes>`. Phenotypes are simulated from genotypes output by the :doc:`transform command </commands/transform>`.

The implementation is based on the `GCTA GWAS Simulation <https://yanglab.westlake.edu.cn/software/gcta/#GWASSimulation>`_ utility.

Usage
~~~~~
.. code-block:: bash

   haptools simphenotype \
   --replications INT \
   --environment FLOAT \
   --heritability FLOAT \
   --prevalence FLOAT \
   --normalize \
   --region TEXT \
   --sample SAMPLE --sample SAMPLE \
   --samples-file FILENAME \
   --id ID --id ID \
   --ids-file FILENAME \
   --chunk-size INT \
   --repeats PATH \
   --seed INT \
   --output PATH \
   --verbosity [CRITICAL|ERROR|WARNING|INFO|DEBUG|NOTSET] \
   GENOTYPES HAPLOTYPES

Model
~~~~~
Each normalized haplotype :math:`\vec{Z_j}` is encoded as an independent causal variable in a linear model:

.. math::

   \vec{y} = \sum_j \beta_j \vec{Z_j} + \vec \epsilon

where

.. math::

   \epsilon_i \sim N(0, \sigma^2)

.. math::

   \sigma^2 = v (\frac 1 {h^2} - 1)

The variable :math:`v` can be specified via the ``--environment`` parameter. When not provided, :math:`v` is inferred from the variance of the genotypes:

.. math::

   v = Var[\sum_j \beta_j \vec{Z_j}]

The heritability :math:`h^2` can be specified via the ``--heritability`` parameter and defaults to 0.5 when not provided.

When both :math:`v` and :math:`h^2` aren't provided, :math:`\sigma^2` is computed purely from the effect sizes, instead:

.. math::

   \sigma^2 = \Biggl \lbrace {1 - \sum \beta_j^2 \quad \quad {\sum \beta_j^2 \le 1} \atop 0 \quad \quad \quad \quad \quad \text{ otherwise }}

If a prevalence for the disease is specified via the ``--prevalence`` parameter, the final :math:`\vec{y}` is thresholded to produce a binary case/control trait with the desired fraction of diseased individuals.

Input
~~~~~
Genotypes must be specified in VCF and haplotypes must be specified in the :doc:`.snplist </formats/snplist>` or :doc:`.hap file format </formats/haplotypes>`.

.. note::
   Your ``.hap`` files must contain a "beta" extra field. See :ref:`this section <formats-haplotypes-extrafields-simphenotype>` of the ``.hap`` format spec for more details.

Alternatively, you may also specify genotypes in PLINK2 PGEN format. Just use the appropriate ".pgen" file extension in the input. See the documentation for genotypes in :ref:`the format docs <formats-genotypesplink>` for more information.

Output
~~~~~~
Phenotypes are output in the PLINK2-style ``.pheno`` file format. If ``--replications`` was set to greater than 1, additional columns are output for each simulated trait.

.. note::
   Case/control phenotypes are encoded as 0 (control) + 1 (case) **not** 1 (control) + 2 (case). The latter is assumed by PLINK2 unless the ``--1`` flag is used (see `the PLINK2 docs <https://www.cog-genomics.org/plink/2.0/input#input_missing_phenotype>`_). Therefore, you must use ``--1`` when providing our ``.pheno`` files to PLINK.

Examples
~~~~~~~~
In its simplest usage, ``simphenotype`` can be used to simulate traits arising from SNPs in a :doc:`.snplist file </formats/snplist>`.

.. code-block:: bash

   haptools simphenotype tests/data/apoe.vcf.gz tests/data/apoe.snplist

However, if you want to simulate haplotype-based effects, you will need to ``transform`` your SNPs into haplotypes first. You can pass the same ``.hap`` file to both commands.

.. code-block:: bash

   haptools transform tests/data/simple.vcf tests/data/simple.hap | \
   haptools simphenotype -o simulated.pheno /dev/stdin tests/data/simple.hap

By default, all of the effects in the ``.hap`` file will be encoded as causal variables. Alternatively, you can select the causal variables manually via the ``--id`` or ``--ids-file`` parameters.

.. code-block:: bash

   haptools transform tests/data/simple.vcf tests/data/simple.hap | \
   haptools simphenotype --id 'H1' /dev/stdin tests/data/simple.hap

To simulate ancestry-specific effects from a genotypes file with population labels, use the ``--ancestry`` switch when running ``transform``.

.. code-block:: bash

   haptools transform --ancestry tests/data/simple-ancestry.vcf tests/data/simple.hap | \
   haptools simphenotype --id 'H1' /dev/stdin tests/data/simple.hap

If speed is important, it's generally faster to use PGEN files than VCFs.

.. code-block:: bash

   haptools transform -o simple-haps.pgen tests/data/simple.pgen tests/data/simple.hap
   haptools simphenotype --id 'H1' simple-haps.pgen tests/data/simple.hap

To simulate causal tandem repeats we require an 'R' line in the **.hap** file and a genotypes file with repeats instead of haplotypes.

.. code-block:: bash

   haptools simphenotype --id 1:10114:GTT tests/data/simple_tr.vcf tests/data/simple_tr.hap

.. note::
   If you would like to simulate from a mix of both haplotypes and repeats, you should specify your repeats in a separate file via the ``--repeats`` argument.

Let's simulate two replicates of a case/control trait that occurs in 60% of samples with a heritability of 0.8. We'll encode only two of the haplotypes in ``tests/data/simphenotype.hap`` as independent causal variables.

.. code-block:: bash

   haptools transform tests/data/example.vcf.gz tests/data/simphenotype.hap | \
   haptools simphenotype \
   --replications 2 \
   --heritability 0.8 \
   --prevalence 0.6 \
   --id 'chr21.q.3365*10' \
   --id 'chr21.q.3365*11' \
   --output simulated.pheno \
   /dev/stdin tests/data/simphenotype.hap

All files used in these examples are described :doc:`here </project_info/example_files>`.

Detailed Usage
~~~~~~~~~~~~~~

.. click:: haptools.__main__:main
   :prog: haptools
   :nested: full
   :commands: simphenotype
