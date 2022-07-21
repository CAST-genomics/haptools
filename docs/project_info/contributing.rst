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
    5. Create a new branch off of master with :code:`git checkout -b <descriptive_branch_name>`. Please follow best practices when naming your branch
    6. Setup our development environment by following the instructions in :ref:`dev-setup-instructions` below
    7. Make your changes to the code
    8. Add additional tests to the :code:`tests/` directory and add comments to the documentation to explain how to use your new code. We use pytest for testing and sphinx/numpydoc for documentation
    9. Run the automated code-checking steps detailed in :ref:`code-check-instructions` below
    10. Commit your changes. Please use informative commit messages and do your best to ensure the commit history is clean and easy to interpret
    11. Now you can push your changes to your Github copy of haptools by running :code:`git push origin <descriptive_branch_name>`
    12. Go to your Github copy of haptools in your browser and create a pull request. Be sure to change the pull request target branch to :code:`main` on this original repository!
    13. Please write an informative pull request detailing the changes you have made and why you made them. Tag any related issues by referring to them by a hashtag followed by their ID


.. _dev-setup-instructions:

------------
Dev Setup
------------

Follow these steps to set up a development environment.

1. Create a conda environment with ``poetry``

    .. code-block:: bash

        conda create -n haptools-dev -c conda-forge 'poetry>=1.1.14' 'python=3.7'
2. Activate the environment

    .. code-block:: bash

        conda activate haptools-dev
3. Install our dependencies

    .. code-block:: bash

        poetry install -E docs -E test

---------------------
Managing Dependencies
---------------------
Run ``poetry help`` to read about the suite of commands it offers for managing dependencies.

For example, to add a pypi dependency to our list and install it, just run

    .. code-block:: bash

        poetry add <dependency>

.. _code-check-instructions:

-----------
Code Checks
-----------
Before creating your pull request, please do the following.

1. Format the code correctly

    .. code-block:: bash

        black .

2. If you made changes to the docs, check that they appear correctly.

    .. code-block:: bash

        ( cd docs && sphinx-build -M html . _build )
        open docs/_build/html/index.html

3. Run all of the tests

    .. code-block:: bash

        pytest tests/

-----
Style
-----
~~~~
Code
~~~~

    1. Please type-hint all function parameters
    2. Please adhere to PEP8 whenever possible. :code:`black` will help you with this.
    3. For readability, please separate imports into three paragraph blocks:
        i. from the python standard library
        ii. from external, third party packages
        iii. from our own internal code

~~~~~~~~~~~~~~~~~~~
Git commit messages
~~~~~~~~~~~~~~~~~~~

    1. Use the present tense ("Add feature" not "Added feature")
    2. Use the imperative mood ("Move cursor to..." not "Moves cursor to...")
    3. Reference issues and pull requests liberally after the first line
