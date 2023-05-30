.. _commands-index:


index
=====

Index a set of haplotypes specified as a :doc:`.hap file </formats/haplotypes>`.

The ``index`` command creates a sorted ``.hap.gz`` and a ``.hap.gz.tbi`` index file from a ``.hap`` (or ``.hap.gz``) file.

By default, the ``index`` command will also sort your haplotypes, since this is a prerequisite for indexing.

Usage
~~~~~
.. code-block:: bash

  haptools index \
  --sort \
  --output PATH \
  --verbosity [CRITICAL|ERROR|WARNING|INFO|DEBUG|NOTSET] \
  HAPLOTYPES

Examples
~~~~~~~~
.. code-block:: bash

  haptools index tests/data/basic.hap

You may also specify a custom output path for the compressed file to be written to.

.. code-block:: bash

  haptools index --output tests/data/sorted.basic.hap.gz tests/data/basic.hap

You can use the ``--no-sort`` flag to skip the sorting step if your file is already sorted.

.. code-block:: bash

  haptools index --no-sort --output tests/data/basic.hap.gz tests/data/basic.hap.gz

.. warning::
  If the ``--no-sort`` flag *isn't* used, the ``index`` command will ignore all extra fields when processing your ``.hap`` file. To retain them, just sort the file manually first.

  .. code-block:: bash

    awk '$0 ~ /^#/ {print; next} {print | "sort -k2,4"}' tests/data/simphenotype.hap | \
    haptools index --no-sort --output tests/data/simphenotype.hap.gz /dev/stdin

All files used in these examples are described :doc:`here </project_info/example_files>`.


Detailed Usage
~~~~~~~~~~~~~~

.. click:: haptools.__main__:main
   :prog: haptools
   :show-nested:
   :commands: index
