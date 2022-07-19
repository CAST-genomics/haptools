.. _commands-simphenotype:


simphenotype
============

Simulates a complex trait, taking into account haplotype- or local-ancestry- specific effects as well as traditional variant-level effects. The user denotes causal variables to use within the simulation by specifying them in a ``.hap`` file.

The implementation is based on the `GCTA GWAS Simulation <https://yanglab.westlake.edu.cn/software/gcta/#GWASSimulation>`_ utility.

Usage
~~~~~
.. code-block:: bash

   haptools simphenotype \
   --replications INT \
   --heritability FLOAT \
   --prevalence FLOAT \
   --region TEXT \
   --sample SAMPLE \
   --samples-file FILENAME \
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

Examples
~~~~~~~~
.. code-block:: bash

   haptools simphenotype -o simulated.pheno tests/data/example.vcf.gz tests/data/example.hap.gz

Simulate two replicates of a case/control trait that occurs in 60% of your samples with a heritability of 0.8. Encode all of the haplotypes in ``tests/data/example.hap.gz`` as independent causal variables.

.. code-block:: bash

   haptools simphenotype \
   --replications 2 \
   --heritability 0.8 \
   --prevalence 0.6 \
   --output simulated.pheno \
   tests/data/example.vcf.gz tests/data/example.hap.gz

Detailed Usage
~~~~~~~~~~~~~~

.. click:: haptools.__main__:main
   :prog: haptools
   :show-nested:
   :commands: simphenotype
