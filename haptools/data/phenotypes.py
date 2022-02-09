from __future__ import annotations
import numpy as np
from csv import reader
from pathlib import Path

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

    Examples
    --------
    >>> phenotypes = Phenotypes.load('tests/data/simple.tsv')
    """

    def __init__(self, fname: Path):
        super().__init__(fname)
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
        with open(self.fname) as phens:
            phen_text = reader(phens, delimiter="\t")
            # convert to list and subset samples if need be
            if samples:
                samples = set(samples)
                phen_text = [phen for phen in phen_text if phen[0] in samples]
            else:
                phen_text = list(phen_text)
        # there should only be two columns
        assert len(phen_text[0]) == 2, "The phenotype TSV should only have two columns."
        # the second column should be castable to a float
        try:
            float(phen_text[0][1])
        except:
            raise AssertionError("The second column of the TSV file must numeric.")
        # fill out the samples and data properties
        self.samples, self.data = zip(*phen_text)
        # coerce strings to floats
        self.data = np.array(self.data, dtype="float64")

    def standardize(self):
        """
        Standardize phenotypes so they have a mean of 0 and a stdev of 1

        This function modifies :py:attr:`~.Genotypes.data` in-place
        """
        self.data = (self.data - np.mean(self.data)) / np.std(self.data)
