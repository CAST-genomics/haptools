from __future__ import annotations
from csv import reader
from pathlib import Path
from collections import namedtuple
from logging import getLogger, Logger
from fileinput import hook_compressed

import numpy as np

from .data import Data
from .phenotypes import Phenotypes


class Covariates(Phenotypes):
    """
    A class for processing covariates from a file

    Attributes
    ----------
    data : np.array
        The covariates in an n (samples) x m (covariates) array
    fname : Path
        The path to the read-only file containing the data
    samples : tuple[str]
        The names of each of the n samples
    names : tuple[str]
        The names of the covariates
    log: Logger
        A logging instance for recording debug statements.

    Examples
    --------
    >>> covariates = Covariates.load('tests/data/covars.tsv')
    """

    def __init__(self, fname: Path, log: Logger = None):
        super(Phenotypes, self).__init__(fname, log)
        self._ext = 'covar'
