.. _project_info-installation:

============
Installation
============

Using pip
---------

You can install ``haptools`` from PyPI using ``pip``.

.. code-block:: bash

   pip install haptools

.. warning::
   We recommend using ``pip >= 22.2.2`` because of `an issue in pysam <https://github.com/pysam-developers/pysam/issues/1132>`_.

Installing ``haptools`` with the "files" extra requirements enables automatic support for a variety of additional file formats, like PLINK2 PGEN files.

.. code-block:: bash

   pip install 'haptools[files]'

.. note::
   Installing ``haptools`` with the "files" extra requirement requires ``gcc`` and a few other compiler tools. Please make sure that they are installed first. To install with conda, for example, please execute the following:

   .. code-block:: bash

      conda install -c conda-forge gxx_linux-64

   Alternatively, you can use the following on Ubuntu:

   .. code-block:: bash

      sudo apt install build-essential

   See `issue 217 <https://github.com/chrchang/plink-ng/issues/217>`_ for current progress on this problem.

Using conda
-----------

We also support installing ``haptools`` from bioconda using ``conda``.

.. code-block:: bash

   conda install -c conda-forge -c bioconda haptools

.. note::
   Installing ``haptools`` from bioconda with PGEN support is not yet possible. See `issue 228 <https://github.com/chrchang/plink-ng/issues/228>`_ for current progress on this challenge.

Installing the latest, unreleased version
-----------------------------------------
Can't wait for us to tag and release our most recent updates? You can install ``haptools`` directly from the ``main`` branch of our Github repository using ``pip``.

.. code-block:: bash

   pip install --upgrade --force-reinstall git+https://github.com/cast-genomics/haptools.git
