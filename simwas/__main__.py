#!/usr/bin/env python

import click

from pathlib import Path
from typing import Union, Tuple

from . import data, tree


@click.group()
@click.version_option()
def main():
    """
    simwas: Simulate phenotypes for GWAS and subsequent fine-mapping

    Use real variants to simulate real, biological LD patterns.
    """
    pass


@main.command()
def simulate():
    """
    Use the tool to simulate genotypes
    """
    pass


if __name__ == "__main__":
    # run the CLI if someone tries 'python -m simwas' on the command line
    main(prog_name="simwas")
