.. _project_info-installation:

============
Installation
============

You can install ``haptools`` from PyPI.

.. code-block:: bash

   pip install haptools

Installing ``haptools`` with the "files" extra requirements enables automatic support for a variety of additional file formats, like PLINK2 PGEN files.

.. code-block:: bash

   pip install haptools[files]

We also support installing ``haptools`` from Bioconda.

.. code-block:: bash

   conda install -c conda-forge -c bioconda haptools

.. note::
   Installing ``haptools`` from bioconda with PGEN support is not yet possible. See `this thread <https://github.com/chrchang/plink-ng/issues/228>`_ for current progress on this issue.
