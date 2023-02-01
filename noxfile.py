"""Nox sessions."""
import os
import sys
import shutil
from pathlib import Path

import nox
from nox_poetry import Session
from nox_poetry import session


package = "haptools"
python_versions = ["3.8", "3.9", "3.10", "3.11"]
nox.needs_version = ">= 2021.6.6"
nox.options.sessions = (
    "docs",
    "lint",
    "tests",
)


@session(python=python_versions[0])
def docs(session: Session) -> None:
    """Build the documentation."""
    args = session.posargs or ["docs", "docs/_build"]
    if not session.posargs and "FORCE_COLOR" in os.environ:
        args.insert(0, "--color")

    build_dir = Path("docs", "_build")
    if build_dir.exists():
        shutil.rmtree(build_dir)

    session.run("sphinx-build", *args)


@session(python=python_versions[0])
def lint(session: Session) -> None:
    """Lint our code."""
    session.install("black")
    session.run("black", "--check", ".")


# detect whether conda/mamba is installed
if os.getenv("CONDA_EXE"):
    conda_cmd = "conda"
    if (Path(os.getenv("CONDA_EXE")).parent / "mamba").exists():
        conda_cmd = "mamba"
    conda_args = ["-c", "conda-forge"]

    @session(venv_backend=conda_cmd, venv_params=conda_args, python=python_versions)
    def tests(session: Session) -> None:
        """Run the test suite."""
        session.conda_install(
            "coverage[toml]", "pytest", "numpy>=1.20.0", channel="conda-forge"
        )
        # TODO: change this to ".[files]" once plink-ng Alpha 3.8 is released
        # https://github.com/chrchang/plink-ng/releases
        session.install(".")

        try:
            session.run(
                "coverage", "run", "--parallel", "-m", "pytest", *session.posargs
            )
        finally:
            if session.interactive:
                session.notify("coverage", posargs=[])

else:

    @session(python=python_versions)
    def tests(session: Session) -> None:
        """Run the test suite."""
        session.install("coverage[toml]", "pytest")
        # TODO: change this to ".[files]" once plink-ng Alpha 3.8 is released
        # https://github.com/chrchang/plink-ng/releases
        session.install(".")

        try:
            session.run(
                "coverage", "run", "--parallel", "-m", "pytest", *session.posargs
            )
        finally:
            if session.interactive:
                session.notify("coverage", posargs=[])


@session(python=python_versions[0])
def coverage(session: Session) -> None:
    """Produce the coverage report."""
    args = session.posargs or ["report"]

    if not session.posargs and any(Path().glob(".coverage.*")):
        session.run("coverage", "combine")

    session.run("coverage", *args)
