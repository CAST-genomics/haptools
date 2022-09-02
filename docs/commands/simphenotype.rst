.. _commands-simphenotype:


simphenotype
============

Simulates a complex trait, taking into account haplotype- or local-ancestry- specific effects as well as traditional variant-level effects. The user denotes causal variables to use within the simulation by specifying them in a :doc:`.hap file </formats/haplotypes>`. Phenotypes are simulated from genotypes output by the :doc:`transform command </commands/transform>`.

To encode simple SNPs as causal variants within a ``.hap`` file, use the haptools API like in :ref:`this example <api-examples-snps2hap>`.

The implementation is based on the `GCTA GWAS Simulation <https://yanglab.westlake.edu.cn/software/gcta/#GWASSimulation>`_ utility.

.. note::
   Your ``.hap`` files must contain a "beta" extra field. See :ref:`this section <formats-haplotypes-extrafields-simphenotype>` of the ``.hap`` format spec for more details.

Usage
~~~~~
.. code-block:: bash

   haptools simphenotype \
   --replications INT \
   --heritability FLOAT \
   --prevalence FLOAT \
   --region TEXT \
   --sample SAMPLE --sample SAMPLE \
   --samples-file FILENAME \
   --id ID --id ID \
   --ids-file FILENAME \
   --chunk-size INT \
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

   \sigma^2 = Var[\sum_j \beta_j \vec{Z_j}] * (\frac 1 {h^2} - 1)

The heritability :math:`h^2` is user-specified, but if it is not provided, then :math:`\sigma^2` will be computed purely from the effect sizes, instead:

.. math::

   \sigma^2 = \Biggl \lbrace {1 - \sum \beta_j^2 \quad \quad {\sum \beta_j^2 \le 1} \atop 0 \quad \quad \quad \quad \quad \text{ otherwise }}

If a prevalence for the disease is specified, the final :math:`\vec{y}` value will be thresholded to produce a binary case/control trait with the desired fraction of diseased individuals.

Output
~~~~~~
Phenotypes are output in the PLINK2-style ``.pheno`` file format. If ``--replications`` was set to greater than 1, additional columns are output for each simulated trait.

Note that case/control phenotypes are encoded as 0 (control) + 1 (case) **not** 1 (control) + 2 (case). The latter is used by PLINK2 unless the ``--1`` flag is used (see `the PLINK2 docs <https://www.cog-genomics.org/plink/2.0/input#input_missing_phenotype>`_). Therefore, you must use ``--1`` when providing our ``.pheno`` files to PLINK.

Examples
~~~~~~~~
.. code-block:: bash

   haptools transform tests/data/example.vcf.gz tests/data/simphenotype.hap | \
   haptools simphenotype -o simulated.pheno /dev/stdin tests/data/simphenotype.hap

By default, all of the haplotypes in the ``.hap`` file will be encoded as causal variables. Alternatively, you can select the causal variables manually via the ``--id`` or ``--ids-file`` parameters.

.. code-block:: bash

   haptools transform tests/data/example.vcf.gz tests/data/simphenotype.hap | \
   haptools simphenotype --id 'chr21.q.3365*1' /dev/stdin tests/data/simphenotype.hap

Simulate two replicates of a case/control trait that occurs in 60% of your samples with a heritability of 0.8. Encode all of the haplotypes in ``tests/data/example.hap.gz`` as independent causal variables.

.. code-block:: bash

   haptools transform tests/data/example.vcf.gz tests/data/simphenotype.hap | \
   haptools simphenotype \
   --replications 2 \
   --heritability 0.8 \
   --prevalence 0.6 \
   --output simulated.pheno \
   /dev/stdin tests/data/example.hap.gz

Detailed Usage
~~~~~~~~~~~~~~

.. click:: haptools.__main__:main
   :prog: haptools
   :nested: full
   :commands: simphenotype
