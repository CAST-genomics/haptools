.. _project_info-contributing:

============
Contributing
============

Contributions are welcome and greatly appreciated!


----------------------
Types of Contributions
----------------------
~~~~~~~~~~~~
Report a bug
~~~~~~~~~~~~
If you have found a bug, please report it on `our issues page <https://github.com/aryarm/haptools/issues>`_ rather than emailing us directly. Others may have the same issue and this helps us get that information to them.

Before you submit a bug, please search through our issues to ensure it hasn't already been reported. If you encounter an issue that has already been reported, please upvote it by reacting with a thumbs-up emoji. This helps us prioritize the issue.

The most helpful Github issues include
    - the version of haptools you are using, although it's best to use the latest version
    - detailed steps to help us reproduce your error, ideally with the example datasets in the :code:`tests/data` directory

~~~~~~~~~
Fix a bug
~~~~~~~~~
Look through our issues page for bugs. We especially need help with bugs labeled "help wanted". If you want to start working on a bug then please write a message within the thread for that issue on our issues page, so that no one is duplicating work.

~~~~~~~~~~~~~~~~~~~~~~~
Implement a new feature
~~~~~~~~~~~~~~~~~~~~~~~
Our issues page will almost always have features on our wishlist. Once again, if you want to start working on a feature then please write a message within the thread for that feature on our issues page, so that no one is duplicating work.

Have an idea for a new feature that isn't on our wishlist? We'd love to hear about it! Before getting to work, please create a Github issue for it, so that you can make sure we're in agreement about what it should do.

-------------------------------------------
How to fix a bug or implement a new feature
-------------------------------------------
Please create a pull request! A PR is a collection of changes that you have made to the code that we can review and potentially integrate into haptools.

To create a pull request you need to do these steps:
    1. Create a Github account
    2. `Fork the repository <https://docs.github.com/en/get-started/quickstart/fork-a-repo#forking-a-repository>`_
        - Click the "Fork" button in the top, right corner
        - Or, if you had already forked the repository a while ago, `sync your fork <https://docs.github.com/en/github/collaborating-with-pull-requests/working-with-forks/syncing-a-fork>`_ to make sure you're working with the latest version of haptools
    3. `Clone your fork locally <https://docs.github.com/en/get-started/quickstart/fork-a-repo#cloning-your-forked-repository>`_
    4. :code:`cd haptools` into the new directory
    5. Create a new branch off of the :code:`main` branch with :code:`git checkout -b <descriptive_branch_name>`. Please follow `best practices <https://www.conventionalcommits.org/>`_ when naming your branch
    6. Setup our development environment by following the instructions in :ref:`dev-setup-instructions` below
    7. Make your changes to the code
    8. Add additional tests to the :code:`tests/` directory and add comments to the documentation to explain how to use your new code. We use pytest for testing and sphinx/numpydoc for documentation. If you add example code or an example command to the documentation, you should make sure to create an automated test that executes it, as well.
    9. Run the automated code-checking steps detailed in :ref:`code-check-instructions` below
    10. Commit your changes. Please use informative commit messages and do your best to ensure the commit history is clean and easy to interpret
    11. Now you can push your changes to your Github copy of haptools by running :code:`git push origin <descriptive_branch_name>`
    12. Go to your Github copy of haptools in your browser and create a pull request titled according to the `conventional commits spec <https://www.conventionalcommits.org/>`_. Be sure to change the pull request target branch to :code:`main` on this original repository.
    13. Please write an informative pull request detailing the changes you have made and why you made them. Tag any related issues by referring to them by a hashtag followed by their ID


.. _dev-setup-instructions:

------------
Dev Setup
------------

Follow these steps to set up a development environment.

1. Create a conda environment with ``poetry``

    .. code-block:: bash

        conda create -n haptools-dev -c conda-forge 'poetry>=1.1.14' 'python=3.7' 'gxx_linux-64'
2. Activate the environment

    .. code-block:: bash

        conda activate haptools-dev
3. Install our dependencies

    .. code-block:: bash

        poetry install -E docs -E tests -E files

Now, try importing ``haptools`` or running it on the command line.

.. note::
   If you run into an attribute error ``module 'distutils' has no attribute 'util'`` see `this workaround <https://github.com/python-poetry/poetry/issues/3336#issuecomment-831789763>`_.


---------------------
Managing Dependencies
---------------------
Run ``poetry help`` to read about the suite of commands it offers for managing dependencies.

For example, to add a pypi dependency to our list and install it, just run

    .. code-block:: bash

        poetry add <dependency>

You should specify a `version constraint <https://python-poetry.org/docs/master/dependency-specification>`_ when adding a dependency. Use the oldest version compatible with your code. Don't worry if you're not sure at first -- you can (and should!) always update it later. For example, to specify a version of ``click`` >= 8.0.4:

    .. code-block:: bash

        poetry add 'click>=8.0.4'


------------------------------------------
Modifying our command line interface (CLI)
------------------------------------------
We use the `click library <https://click.palletsprojects.com/>`_ to define ``haptools``'s command line interface as `nested commands <https://click.palletsprojects.com/quickstart/#nesting-commands>`_. All of the CLI logic is defined in `__main__.py <https://github.com/CAST-genomics/haptools/blob/main/haptools/__main__.py>`_.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Add or modify a command-line option
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
First, locate the definition of the command in `__main__.py <https://github.com/CAST-genomics/haptools/blob/main/haptools/__main__.py>`_

You can add a ``@click.option`` or ``@click.argument`` line if you want to add a new option or argument. Please follow `click's convention <https://click.palletsprojects.com/parameters/#parameters>`_ and only use ``@click.argument`` for required arguments and ``@click.option`` for optional ones. See `the click documentation <https://click.palletsprojects.com/#documentation>`_ for directions on modifying or adding parameters like options/arguments.

~~~~~~~~~~~~~~~~~
Add a new command
~~~~~~~~~~~~~~~~~
To add a new command, you only have to define a new function in `__main__.py <https://github.com/CAST-genomics/haptools/blob/main/haptools/__main__.py>`_. Within that function, you can import and call the rest of your code. For example, to add a command called ``mycommand`` which takes a single required file called ``arg1``, you might do the following.

.. code-block:: python

    @main.command(short_help="A short description of my command")
    @click.argument("arg1", type=click.Path(exists=True, path_type=Path))
    @click.option(
        "-o",
        "--output",
        type=click.Path(path_type=Path),
        default=Path("/dev/stdout"),
        show_default="stdout",
        help="The output of my command",
    )
    @click.option(
        "-v",
        "--verbosity",
        type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]),
        default="INFO",
        show_default=True,
        help="The level of verbosity desired",
    )
    def mycommand(
        arg1: Path,
        output: Path = None,
        verbosity: str = "INFO",
    ):
        """
        A longer description of mycommand
        """

        from .mycommand import run_things
        from .logging import getLogger

        log = getLogger(name="mycommand", level=verbosity)

        run_things(arg1, output, log)

Notice that we usually define a logging object here to use throughout our code. For more information about logging, see the :ref:`section about it below <contributing-style-errors>`. All ``haptools`` commands should use a default verbosity of ``INFO``.

.. _code-check-instructions:

-----------
Code Checks
-----------
Before creating your pull request, please run each of our code checks.

1. Format the code correctly

    .. code-block:: bash

        black .

2. If you made changes to the docs, check that they appear correctly.

    .. code-block:: bash

        sphinx-build docs docs/_build
        open docs/_build/html/index.html

3. Run all of the tests

    .. code-block:: bash

        pytest tests/

    You can also build the package and run the tests from the built version using ``nox``. This will fully simulate installing the package from PyPI.

    .. code-block:: bash

        nox --session=tests

---------------------
Publish a new version
---------------------
To publish a new version of haptools:

1. First, locate `the most recent haptools PR <https://github.com/CAST-genomics/haptools/pulls>`_ prefixed "chore(main)" created by our Github actions bot
2. Comment "/check" on the PR to force the required checks to run and then merge the PR
3. The bot will automatically create a new version on PyPI and tag a release on Github
4. A few hours later, bioconda will automatically detect the new release on PyPI and create a PR in `their repository <https://github.com/bioconda/bioconda-recipes/pulls>`_
5. Check that all of the dependencies in the recipe have been updated properly. If they are, you should comment on the bioconda PR with "@BiocondaBot please add label"
6. After 1-2 days, someone from the bioconda team will merge our PR and the version will get updated on bioconda. Otherwise, ping them a reminder on `Gitter <https://gitter.im/bioconda/Lobby>`_

-----
Style
-----
~~~~
Code
~~~~

    1. Please type-hint all function parameters
    2. Please adhere to PEP8 whenever possible. :code:`black` will help you with this.
    3. Please use relative imports whenever importing modules from the code base
    4. For readability, please separate imports into three paragraph blocks:
        i. from the python standard library
        ii. from external, third party packages
        iii. from our own internal code

.. _contributing-style-errors:

~~~~~~
Errors
~~~~~~
We use the `Python logging module <https://coralogix.com/blog/python-logging-best-practices-tips/>`_ for all messages, including warnings, debugging info, and otherwise. For example, all classes in the ``data`` module have a ``log`` property that stores a logger object. If you are creating a new command, you can use our custom logging module to retrieve a suitable object.

.. code-block:: python

    from .logging import getLogger

    # the level of verbosity desired by the user
    # can be: CRITICAL, ERROR, WARNING, INFO, DEBUG, or NOTSET
    verbosity = "DEBUG"

    # create a new logger object for the transform command
    log = getLogger(name="transform", level=verbosity)

    # log a warning message to the logger
    log.warning("This is a warning")

This way, the user can choose their level of verbosity among *CRITICAL*, *ERROR*, *WARNING*, *INFO*, *DEBUG*, and *NOTSET*. However, for critical errors (especially for those in the ``data`` module), our convention is to raise exceptions, usually with a custom ``ValueError``.

~~~~~~~~~~~~~~~~~~~
Git commit messages
~~~~~~~~~~~~~~~~~~~

    1. Use the present tense ("Add feature" not "Added feature")
    2. Use the imperative mood ("Move cursor to..." not "Moves cursor to...")
    3. Reference issues and pull requests liberally after the first line
