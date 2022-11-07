from __future__ import annotations
from csv import reader
from pathlib import Path
from io import TextIOBase
from collections.abc import Iterable
from logging import getLogger, Logger
from collections import namedtuple, Counter

import numpy as np
import numpy.typing as npt

from .data import Data


class Phenotypes(Data):
    """
    A class for processing phenotypes from a file

    Attributes
    ----------
    data : np.array
        The phenotypes in an n (samples) x m (phenotypes) array
    fname : Path | str
        The path to the file containing the data
    samples : tuple
        The names of each of the n samples
    names : tuple[str]
        The names of the phenotypes
    log: Logger
        A logging instance for recording debug statements.

    Examples
    --------
    >>> phenotypes = Phenotypes.load('tests/data/simple.pheno')
    """

    def __init__(self, fname: Path | str, log: Logger = None):
        super().__init__(fname, log)
        self.samples = tuple()
        self.names = tuple()
        self._ext = "pheno"

    @classmethod
    def load(
        cls: Phenotypes, fname: Path | str, samples: list[str] = None
    ) -> Phenotypes:
        """
        Load phenotypes from a pheno file

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
        Read phenotypes from a pheno file into a numpy matrix stored in :py:attr:`~.Penotypes.data`

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
        # call self.__iter__() to obtain the data and samples
        data, self.samples = zip(*self.__iter__(samples))
        self.log.info(f"Loaded {len(self.samples)} samples from .{self._ext} file")
        # fill out the samples and data properties
        # collect data in a np array
        self.data = np.array(data)

    def _iterate(
        self, phens: TextIOBase, phen_text: Iterable, samples: set[str] = None
    ):
        """
        A generator over the lines of a pheno

        This is a helper function for :py:meth:`~.Phenotypes.__iter__`

        Parameters
        ----------
        phens: TextIOBase
            The file handler for the stream
        phen_text: Iterable
            The csv.reader object containing the lines of text from the file as lists
        samples : set[str], optional
            A subset of the samples from which to extract phenotypes

            Defaults to loading phenotypes from all samples

        Yields
        ------
        Iterator[namedtuple]
            An iterator over each line in the file, where each line is encoded as a
            namedtuple containing each of the class properties
        """
        self.log.info(f"Loading {len(self.names)} columns from .{self._ext} file")
        Record = namedtuple("Record", "data samples")
        for phen in phen_text:
            if samples is None or phen[0] in samples:
                try:
                    yield Record(np.array(phen[1:], dtype="float64"), phen[0])
                except:
                    self.log.error(
                        f"Every column in the .{self._ext} file (besides the sample"
                        " column) must be numeric."
                    )
        phens.close()

    def __iter__(self, samples: list[str] = None) -> Iterable[namedtuple]:
        """
        Read phenotypes from a pheno line by line without storing anything

        Parameters
        ----------
        samples : list[str], optional
            A subset of the samples from which to extract phenotypes

            Defaults to loading phenotypes from all samples

        Returns
        ------
        Iterable[namedtuple]
            See documentation for :py:meth:`~.Phenotypes._iterate`
        """
        phens = self.hook_compressed(self.fname, mode="r")
        phen_text = reader(phens, delimiter="\t")
        # ignore all of the comment lines
        while True:
            header = next(phen_text)
            if not header[0].startswith("#") or header[0].startswith("#IID"):
                break

        # there should be at least two columns
        if len(header) < 2:
            raise ValueError(f"The .{self._ext} file should have at least two columns.")
        # the first column should be called "#IID"
        if header[0] != "#IID":
            self.log.warning(
                f"The first column of the .{self._ext} file should contain sample IDs"
                " and should be named '#IID' in the header line"
            )
        self.names = tuple(header[1:])
        samples = set(samples) if samples else None
        # call another function to force the lines above to be run immediately
        # see https://stackoverflow.com/a/36726497
        return self._iterate(phens, phen_text, samples)

    def write(self):
        """
        Write the phenotypes in this class to a file at :py:attr:`~.Phenotypes.fname`

        Examples
        --------
        To write to a file, you must first initialize a Phenotypes object and then
        fill out the names, data, and samples properties:
        >>> phenotypes = Phenotypes('tests/data/simple.pheno')
        >>> phenotypes.names = ('height',)
        >>> phenotypes.data = np.array([1, 1, 2], dtype='float64')
        >>> phenotypes.samples = ('HG00096', 'HG00097', 'HG00099')
        >>> phenotypes.write()
        """
        # make sure the names are unique
        uniq_names = Counter()
        names = [None] * len(self.names)
        for idx, name in enumerate(self.names):
            suffix = ""
            if uniq_names[name]:
                suffix = f"-{uniq_names[name]}"
            names[idx] = name + suffix
            uniq_names[name] += 1
        # now we can finally write the file
        with self.hook_compressed(self.fname, mode="w") as phens:
            phens.write("#IID\t" + "\t".join(names) + "\n")
            formatter = {"float_kind": lambda x: "%.2f" % x}
            for samp, phen in zip(self.samples, self.data):
                line = np.array2string(
                    phen,
                    separator="\t",
                    formatter=formatter,
                    max_line_width=np.inf,
                    threshold=np.inf,
                    edgeitems=np.inf,
                )[1:-1]
                phens.write(f"{samp}\t" + line + "\n")

    # TODO: check_missing() (they'll be encoded as NA, nan, or -9)

    def standardize(self):
        """
        Standardize phenotypes so they have a mean of 0 and a stdev of 1

        This function modifies :py:attr:`~.Phenotypes.data` in-place
        """
        std = np.std(self.data, axis=0)
        self.data = (self.data - np.mean(self.data, axis=0)) / std
        # for phenotypes where the stdev is 0, just set all values to 0 instead of nan
        zero_elements = std == 0
        self.data[:, zero_elements] = np.zeros(
            (self.data.shape[0], np.sum(zero_elements))
        )

    def append(self, name: str, data: npt.NDArray):
        """
        Append a new set of phenotypes to the current set

        Parameters
        ----------
        name: str
            The name of the new phenotype
        data: npt.NDArray
            A 1D np array of the same length as :py:attr:`~.Phenotypes.samples`,
            containing the phenotype values for each sample. Must have the same dtype
            as :py:attr:`~.Phenotypes.data.`
        """
        if len(self.samples):
            if len(self.samples) != len(data):
                self.log.error(
                    "The data provided to the add() method is not of the appropriate"
                    "length"
                )
        else:
            self.log.warning(
                "Set the samples property of the Phenotypes instance before calling "
                "the add() method"
            )
        if self.unset():
            self.data = data[:, np.newaxis]
        else:
            self.data = np.concatenate((self.data, data[:, np.newaxis]), axis=1)
        self.names = self.names + (name,)
