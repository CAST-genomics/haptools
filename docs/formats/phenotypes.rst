.. _formats-phenotypes:


Phenotypes and Covariates
=========================

Phenotype file format
---------------------
Phenotypes are expected to follow the PLINK2 ``.pheno`` file format. This is a
tab-separated format where the first column corresponds to the sample ID, and
subsequent columns contain each of your phenotypes.

The first line of the file corresponds with the header and must begin with ``#IID``.
The names of each of your phenotypes belong in the subbsequent columns of the header.

See `tests/data/simple.pheno <https://github.com/gymrek-lab/haptools/blob/main/tests/data/simple.pheno>`_ for an example of a phenotype file:

.. include:: ../../tests/data/simple.pheno
   :literal:

Covariate file format
---------------------
Covariates follow the same format as phenotypes.

See `tests/data/simple.covar <https://github.com/gymrek-lab/haptools/blob/main/tests/data/simple.covar>`_ for an example of a covariate file:

.. include:: ../../tests/data/simple.covar
   :literal:
