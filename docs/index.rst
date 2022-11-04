.. _manual-main:

haptools
========

Haptools is a collection of tools for simulating and analyzing genotypes and phenotypes while taking into account haplotype and ancestry information.

We support fast simulation of admixed genomes, visualization of admixture tracks, simulating haplotype- and local ancestry-specific phenotype effects, and computing a variety of common file operations and statistics in a haplotype-aware manner.

At the core of haptools lies the :doc:`.hap file </formats/haplotypes>`, our new file format for haplotypes designed for speed, extensibility, and ease-of-use.

Commands
~~~~~~~~

* :doc:`haptools simgenotype </commands/simgenotype>`: Simulate genotypes for admixed individuals under user-specified demographic histories.

* :doc:`haptools simphenotype </commands/simphenotype>`: Simulate a complex trait, taking into account local ancestry- or haplotype- specific effects. ``haptools simphenotype`` takes as input a VCF file (usually from ``haptools transform``) and outputs simulated phenotypes for each sample.

* :doc:`haptools karyogram </commands/karyogram>`: Visualize a "chromosome painting" of local ancestry labels based on breakpoints output by ``haptools simgenotype``.

* :doc:`haptools transform </commands/transform>`: Transform a set of genotypes via a list of haplotypes. Create a new VCF containing haplotypes instead of variants.

* :doc:`haptools index </commands/index>`: Sort, compress, and index our custom file format for haplotypes.

* :doc:`haptools ld </commands/ld>`: Compute Pearson's correlation coefficient between a target haplotype and a set of haplotypes.

.. figure:: https://drive.google.com/uc?id=1c0i_Hjms7579s24zRsKp5yMs7BxNHed_
  :figwidth: 600
  :align: center
  :alt: Overview of haptools commands

Outputs produced by these utilities are compatible with each other.
For example ``haptools simgenotype`` outputs a VCF file with local ancestry information annotated for each variant.
The VCF and breakpoints file output by ``haptools simgenotype`` can be used as input to ``haptools transform``, which is then used by ``haptools simphenotype`` to simulate phenotypes for a list of haplotypes.
The local ancestry breakpoints from ``haptools simgenotype`` can also be visualized using ``haptools karyogram``.

Detailed information about each command can be found in `the commands section </commands>`_. Examples there utilize files from the `tests directory <https://github.com/CAST-genomics/haptools/tree/main/tests/data>`_.

Contributing
~~~~~~~~~~~~

We gladly welcome any contributions to ``haptools``!

Please read our :doc:`contribution guidelines </project_info/contributing>` and then submit a `Github issue <https://github.com/cast-genomics/haptools/issues>`_.


.. toctree::
   :caption: Overview
   :name: overview
   :hidden:
   :maxdepth: 1

   project_info/installation
   project_info/contributing

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
   commands/index.rst
   commands/ld.rst

.. toctree::
   :caption: API
   :name: api
   :hidden:
   :maxdepth: 1

   api/data
   api/modules
   api/examples
