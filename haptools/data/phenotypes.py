from __future__ import annotations
from csv import reader
from pathlib import Path
from collections import namedtuple
from fileinput import hook_compressed

import numpy as np

from .data import Data


class Phenotypes(Data):
    """
    A class for processing phenotypes from a file

    Attributes
    ----------
    data : np.array
        The phenotypes in an n (samples) x 1 (phenotype value) array
    fname : Path
        The path to the read-only file containing the data
    samples : tuple
        The names of each of the n samples
    log: Logger
        A logging instance for recording debug statements.

    Examples
    --------
    >>> phenotypes = Phenotypes.load('tests/data/simple.tsv')
    """

    def __init__(self, fname: Path, log: Logger = None):
        super().__init__(fname, log)
        self.samples = tuple()

    @classmethod
    def load(cls: Phenotypes, fname: Path, samples: list[str] = None) -> Phenotypes:
        """
        Load phenotypes from a TSV file

        Read the file contents and standardize the phenotypes

        Parameters
        ----------
        fname
            See documentation for :py:attr:`~.Data.fname`
        samples : list[str], optional
            See documentation for :py:meth:`~.Data.Phenotypes.read`

        Returns
        -------
        phenotypes
            A Phenotypes object with the data loaded into its properties
        """
        phenotypes = cls(fname)
        phenotypes.read(samples)
        phenotypes.standardize()
        return phenotypes

    def read(self, samples: list[str] = None):
        """
        Read phenotypes from a TSV file into a numpy matrix stored in :py:attr:`~.Penotypes.data`

        Parameters
        ----------
        samples : list[str], optional
            A subset of the samples from which to extract phenotypes

            Defaults to loading phenotypes from all samples

        Raises
        ------
        AssertionError
            If the provided file doesn't follow the expected format
        """
        super().read()
        # load all info into memory
        # use hook_compressed to automatically handle gz files
        with hook_compressed(self.fname, mode="rt") as phens:
            phen_text = reader(phens, delimiter="\t")
            # convert to list and subset samples if need be
            if samples:
                samples = set(samples)
                phen_text = [phen for phen in phen_text if phen[0] in samples]
            else:
                phen_text = list(phen_text)
        # there should only be two columns
        if len(phen_text[0]) != 2:
            self.log.warning("The phenotype TSV should only have two columns.")
        # the second column should be castable to a float
        try:
            float(phen_text[0][1])
        except:
            self.log.error("The second column of the TSV file must numeric.")
        # fill out the samples and data properties
        self.samples, self.data = zip(*phen_text)
        # coerce strings to floats
        self.data = np.array(self.data, dtype="float64")

    def iterate(self, samples: list[str] = None) -> Iterator[namedtuple]:
        """
        Read phenotypes from a TSV line by line without storing anything

        Parameters
        ----------
        samples : list[str], optional
            A subset of the samples from which to extract phenotypes

            Defaults to loading phenotypes from all samples

        Yields
        ------
        Iterator[namedtuple]
            An iterator over each line in the file, where each line is encoded as a
            namedtuple containing each of the class properties
        """
        with hook_compressed(self.fname, mode="rt") as phens:
            phen_text = reader(phens, delimiter="\t")
            Record = namedtuple("Record", "data samples")
            for phen in phen_text:
                if samples is None or phen[0] in samples:
                    try:
                        yield Record(float(phen[1]), phen[0])
                    except:
                        self.log.error(
                            "The second column of the TSV file must numeric."
                        )

    def standardize(self):
        """
        Standardize phenotypes so they have a mean of 0 and a stdev of 1

        This function modifies :py:attr:`~.Genotypes.data` in-place
        """
        self.data = (self.data - np.mean(self.data)) / np.std(self.data)
