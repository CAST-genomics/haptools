.. _manual-main:

haptools
========

Haptools is a collection of tools for simulating and analyzing genotypes and phenotypes while taking into account haplotype information. It is particularly designed for analysis of individuals with admixed ancestries, although the tools can also be used for non-admixed individuals.

Installation
~~~~~~~~~~~~
.. note::
   To reduce the likelihood of errors, we recommend installing ``haptools`` within a new conda environment using a recent version of pip:

   .. code-block:: bash

      conda create -y -n haptools -c conda-forge 'pip>=22.2.2'
      conda activate haptools

We have not officially published ``haptools`` yet, but in the meantime, you can install it directly from our Github repository.

.. code-block:: bash

   pip install git+https://github.com/cast-genomics/haptools.git

Installing ``haptools`` with the "files" extra requirements enables automatic support for a variety of additional file formats, like PLINK2 PGEN files.

.. code-block:: bash

   pip install git+https://github.com/cast-genomics/haptools.git#egg=haptools[files]

Summary of Commands
~~~~~~~~~~~~~~~~~~~

``haptools`` consists of multiple utilities listed below. Click on a utility to see more detailed usage information.

* `haptools simgenotype </commands/simgenotype>`_: Simulate genotypes for admixed individuals under user-specified demographic histories.

* `haptools simphenotype </commands/simphenotype>`_: Simulate a complex trait, taking into account local ancestry- or haplotype- specific effects. ``haptools simphenotype`` takes as input a VCF file (usually from ``haptools transform``) and outputs simulated phenotypes for each sample.

* `haptools karyogram </commands/karyogram>`_: Visualize a "chromosome painting" of local ancestry labels based on breakpoints output by ``haptools simgenotype``.

* `haptools transform </commands/transform>`_: Transform a set of genotypes via a list of haplotypes. Create a new VCF containing haplotypes instead of variants.

* `haptools ld </commands/ld>`_: Compute Pearson's correlation coefficient between a target haplotype and a set of haplotypes.

Outputs produced by these utilities are compatible with each other.
For example ``haptools simgenotype`` outputs a VCF file with local ancestry information annotated for each variant. The output VCF file can be used as input to ``haptools transform`` and ``haptools simphenotype`` to simulate phenotype information. ``haptools simgenotype`` also outputs a list of local ancestry breakpoints which can be visualized using ``haptools karyogram``.

Contributing
~~~~~~~~~~~~

We gladly welcome any contributions to ``haptools``!

Please read our :doc:`contribution guidelines </project_info/contributing>` and then submit a `Github issue <https://github.com/cast-genomics/haptools/issues>`_.


.. toctree::
   :caption: File Formats
   :name: formats
   :hidden:
   :maxdepth: 1

   formats/genotypes.rst
   formats/haplotypes.rst
   formats/phenotypes.rst
   formats/breakpoints.rst
   formats/sample_info.rst
   formats/models.rst
   formats/maps.rst

.. toctree::
   :caption: Commands
   :name: commands
   :hidden:
   :maxdepth: 1

   commands/simgenotype.rst
   commands/simphenotype.rst
   commands/karyogram.rst
   commands/transform.rst
   commands/ld.rst

.. toctree::
   :caption: API
   :name: api
   :hidden:
   :maxdepth: 1

   api/data
   api/modules
   api/examples

.. toctree::
   :caption: Project Info
   :name: project-info
   :hidden:
   :maxdepth: 1

   project_info/contributing
