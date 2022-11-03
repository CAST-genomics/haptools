.. _project_info-installation:

============
Installation
============

You can install ``haptools`` from PyPI using ``pip``. We recommend using ``pip >= 22.2.2`` because of `an issue in pysam <https://github.com/pysam-developers/pysam/issues/1132>`_.

.. code-block:: bash

   pip install haptools

Installing ``haptools`` with the "files" extra requirements enables automatic support for a variety of additional file formats, like PLINK2 PGEN files.

.. code-block:: bash

   pip install haptools[files]

.. note::
   Installing ``haptools`` with the "files" extra requirement requires ``gcc`` and a few other compiler tools. Please make sure that they are installed first. See `issue 217 <https://github.com/chrchang/plink-ng/issues/217>`_ for current progress on this problem. To install with conda, for example, please execute the following:

   .. code-block:: bash

      conda install -c conda-forge gxx_linux-64

   Alternatively, you can use the following on Ubuntu:

   .. code-block:: bash

      sudo apt install build-essential

We also support installing ``haptools`` from Bioconda.

.. code-block:: bash

   conda install -c conda-forge -c bioconda haptools

.. note::
   Installing ``haptools`` from bioconda with PGEN support is not yet possible. See `issue 228 <https://github.com/chrchang/plink-ng/issues/228>`_ for current progress on this challenge.
