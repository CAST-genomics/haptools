from __future__ import annotations
from csv import reader
from pathlib import Path
from fileinput import hook_compressed

import numpy as np

from .data import Data


class Covariates(Data):
    """
    A class for processing covariates from a file

    Attributes
    ----------
    data : np.array
        The covariates in an n (samples) x 1 (covariate value) array
    fname : Path
        The path to the read-only file containing the data
    samples : tuple[str]
        The names of each of the n samples
    log: Logger
        A logging instance for recording debug statements.

    Examples
    --------
    >>> covariates = Covariates.load('tests/data/covars.tsv')
    """

    def __init__(self, fname: Path, log: Logger = None):
        super().__init__(fname, log)
        self.samples = tuple()
        self.names = tuple()

    @classmethod
    def load(cls: Covariates, fname: Path, samples: list[str] = None) -> Covariates:
        """
        Load covariates from a TSV file

        Read the file contents and standardize the covariates

        Parameters
        ----------
        fname
            See documentation for :py:attr:`~.Data.fname`
        samples : list[str], optional
            See documentation for :py:meth:`~.Data.Covariates.read`

        Returns
        -------
        covariates
            A Covariates object with the data loaded into its properties
        """
        covariates = cls(fname)
        covariates.read(samples)
        return covariates

    def read(self, samples: list[str] = None):
        """
        Read covariates from a TSV file into a numpy matrix stored in :py:attr:`~.Covariates.data`

        Parameters
        ----------
        samples : list[str], optional
            A subset of the samples from which to extract covariates

            Defaults to loading covariates from all samples

        Raises
        ------
        AssertionError
            If the provided file doesn't follow the expected format
        """
        super().read()
        # load all info into memory
        # use hook_compressed to automatically handle gz files
        with hook_compressed(self.fname, mode='rt') as covars:
            covar_text = reader(covars, delimiter="\t")
            header = next(covar_text)
            # there should at least two columns
            assert (
                len(header) >= 2
            ), "The covariates TSV should have at least two columns."
            # the first column should be called "sample"
            assert header[0] == "sample", (
                "The first column of the covariates TSV should contain sample IDs and"
                " should be named 'sample' in the header line"
            )
            # convert to list and subset samples if need be
            if samples:
                samples = set(samples)
                covar_text = [covar for covar in covar_text if covar[0] in samples]
            else:
                covar_text = list(covar_text)
        # the second column should be castable to a float
        try:
            float(covar_text[0][1])
        except:
            raise AssertionError(
                "Every column in the covariates file (besides the sample column) must"
                " be numeric."
            )
        # fill out the samples and data properties
        data = list(zip(*covar_text))
        self.samples = data[0]
        self.names = tuple(header[1:])
        # coerce strings to floats
        self.data = np.transpose(np.array(data[1:], dtype="float64"))
