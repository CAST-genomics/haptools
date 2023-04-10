from __future__ import annotations
from pathlib import Path
from logging import getLogger, Logger

from .phenotypes import Phenotypes


class Covariates(Phenotypes):
    """
    A class for processing covariates from a file

    Attributes
    ----------
    data : np.array
        The covariates in an n (samples) x m (covariates) array
    fname : Path | str
        The path to the read-only file containing the data
    samples : tuple[str]
        The names of each of the n samples
    names : tuple[str]
        The names of the covariates
    log: Logger
        A logging instance for recording debug statements.

    Examples
    --------
    >>> covariates = Covariates.load('tests/data/simple.covar')
    """

    def __init__(self, fname: Path | str, log: Logger = None):
        super().__init__(fname, log)
        self._ext = "covar"
