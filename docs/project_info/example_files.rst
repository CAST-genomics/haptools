.. _project_info-example_files:

=============
Example files
=============

Locating the files
------------------

The examples throughout our documentation make use of two sets of files: test files and example files.

Test files can be found in the `tests/data/ directory <https://github.com/CAST-genomics/haptools/tree/main/tests/data>`_ of our Github repository. These are short, simplified files used exclusively by our automated test suite.

Example files can be found in the `example-files/ directory <https://github.com/CAST-genomics/haptools/tree/main/example-files>`_ of our Github repository. Unlike test files, we expect example files to be useful in your own commands. For example, if you use simgenotype with the 1000 Genomes dataset, you can use our `1000G sample_info file  <https://github.com/cast-genomics/haptools/blob/main/example-files/1000genomes_sampleinfo.tsv>`_. We have also included a set of `model files <https://github.com/cast-genomics/haptools/blob/main/example-files/models>`_ that you can use to create pre-configured admixed populations.

.. _running-an-example-command:

Running an example command
--------------------------
To run any of the example code or commands in our documentation, follow these steps.

1. :doc:`Install haptools </project_info/installation>`
2. Clone our Github repository

    .. code-block:: bash

    	git clone -b v$(haptools --version) https://github.com/CAST-genomics/haptools.git

3. Change to the cloned directory

    .. code-block:: bash

    	cd haptools

4. Execute the example command

Running all examples
--------------------
All of our examples are included within our test suite, which is executed regularly by our continuous integration system. To check that all of the examples work on your system, you can just have ``pytest`` automatically run all of our tests.

1. Follow the :ref:`first three steps above <running-an-example-command>`
2. Install ``pytest``

    .. code-block:: bash

    	pip install 'pytest>=6.2.5'

3. Run our tests

    .. code-block:: bash

    	pytest tests/
