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
        self._samp_idx = None
        self._name_idx = None

    @classmethod
    def load(
        cls: Phenotypes, fname: Path | str, samples: set[str] = None
    ) -> Phenotypes:
        """
        Load phenotypes from a pheno file

        Read the file contents and standardize the phenotypes

        Parameters
        ----------
        fname
            See documentation for :py:attr:`~.Data.fname`
        samples : set[str], optional
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

    def read(self, samples: set[str] = None):
        """
        Read phenotypes from a pheno file into a numpy matrix stored in :py:attr:`~.Penotypes.data`

        Parameters
        ----------
        samples : set[str], optional
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

    def __iter__(self, samples: set[str] = None) -> Iterable[namedtuple]:
        """
        Read phenotypes from a pheno line by line without storing anything

        Parameters
        ----------
        samples : set[str], optional
            A subset of the samples from which to extract phenotypes

            Defaults to loading phenotypes from all samples

        Returns
        -------
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
        >>> phenotypes.data = np.array([[1, 1, 2]], dtype='float64')
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
        # check that the phenotypes aren't a single vector; they must have a 2D shape!
        # otherwise, negative signs won't get written correctly to the output
        if len(self.data.shape) <= 1:
            raise ValueError("The data property must have a 2D shape.")
        # now we can finally write the file
        with self.hook_compressed(self.fname, mode="w") as phens:
            phens.write("#IID\t" + "\t".join(names) + "\n")
            for samp, phen in zip(self.samples, self.data):
                line = np.array2string(
                    phen,
                    sign="-",
                    legacy=False,
                    separator="\t",
                    threshold=np.inf,
                    edgeitems=np.inf,
                    floatmode="unique",
                    suppress_small=False,
                    max_line_width=np.inf,
                )[1:-1]
                phens.write(f"{samp}\t" + line + "\n")

    def check_missing(self, discard_also=False):
        """
        Check that each sample has a phenotype value

        Raises
        ------
        ValueError
            If any of the samples have missing phenotypes, represented by -9

        Parameters
        ----------
        discard_also : bool, optional
            If True, discard any samples that are missing phenotypes without raising a
            ValueError
        """
        if len(self.data.shape) <= 1:
            raise ValueError("The data property must have a 2D shape.")
        # check: are there any samples that have phenotypes values that are -9?
        mask = self.data == -9
        missing = np.any(mask, axis=1)
        if np.any(missing):
            samp_idx = np.nonzero(missing)[0]
            missing_phens = np.nonzero(np.any(mask, axis=0))[0]
            if discard_also:
                original_num_samples = len(self.samples)
                self.data = np.delete(self.data, samp_idx, axis=0)
                self.samples = tuple(np.delete(self.samples, samp_idx))
                self.log.warning(
                    "Ignoring missing phenotypes from "
                    f"{original_num_samples - len(self.samples)} samples"
                )
                self._samp_idx = None
            else:
                raise ValueError(
                    "Sample with ID {} for phenotype '{}' is missing".format(
                        self.samples[samp_idx[0]],
                        self.names[missing_phens[0]],
                    )
                )
        if discard_also and not self.data.shape[0]:
            self.log.warning(
                "All samples were discarded! Check that that none of your samples have"
                " missing phenotypes (a value of -9, NA, or na)."
            )

    def standardize(self):
        """
        Standardize phenotypes so they have a mean of 0 and a stdev of 1

        This function modifies :py:attr:`~.Phenotypes.data` in-place
        """
        if len(self.data.shape) <= 1:
            raise ValueError("The data property must have a 2D shape.")
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
        if len(data.shape) != 1:
            raise ValueError("The data argument must have a 1D shape.")
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
            if len(self.data.shape) <= 1:
                raise ValueError("The data property must have a 2D shape.")
            self.data = np.concatenate((self.data, data[:, np.newaxis]), axis=1)
        if self._name_idx is not None:
            self._name_idx[name] = len(self.names)
        self.names = self.names + (name,)

    def index(self, samples: bool = True, names: bool = True):
        """
        Call this function once to improve the amortized time-complexity of look-ups of
        samples and names by their ID. This is useful if you intend to later subset
        by a set of samples or name IDs.
        The time complexity of this function should be roughly O(n+p) if both
        parameters are True. Otherwise, it will be either O(n) or O(p).

        Parameters
        ----------
        samples: bool, optional
            Whether to index the samples for fast loop-up. Adds complexity O(n).
        names: bool, optional
            Whether to index the names for fast look-up. Adds complexity O(p).

        Raises
        ------
        ValueError
            If any samples or names appear more than once
        """
        if samples and self._samp_idx is None:
            self._samp_idx = dict(zip(self.samples, range(len(self.samples))))
            if len(self._samp_idx) < len(self.samples):
                duplicates = Counter(self.samples).items()
                duplicates = [samp_id for samp_id, count in duplicates if count > 1]
                a_few = 5 if len(duplicates) > 5 else len(duplicates)
                raise ValueError(f"Found duplicate sample IDs: {duplicates[:a_few]}")
        if names and self._name_idx is None:
            self._name_idx = dict(zip(self.names, range(len(self.names))))
            if len(self._name_idx) < len(self.names):
                duplicates = Counter(self.names).items()
                duplicates = [name_id for name_id, count in duplicates if count > 1]
                a_few = 5 if len(duplicates) > 5 else len(duplicates)
                raise ValueError(f"Found duplicate name IDs: {duplicates[:a_few]}")

    def subset(
        self,
        samples: tuple[str] = None,
        names: tuple[str] = None,
        inplace: bool = False,
    ):
        """
        Subset these phenotypes to a smaller set of samples or a smaller set of names

        The order of the samples and names in the subsetted instance will match
        the order in the provided tuple parameters.

        Parameters
        ----------
        samples: tuple[str]
            A subset of samples to keep
        names: tuple[str]
            A subset of phenotype IDs to keep
        inplace: bool, optional
            If False, return a new Phenotypes object; otherwise, alter the current one

        Returns
        -------
            A new Phenotypes object if inplace is set to False, else returns None
        """
        # First, initialize variables
        pts = self
        if not inplace:
            pts = self.__class__(self.fname, self.log)
        pts.samples = self.samples
        pts.names = self.names
        pts.data = self.data
        # Index the current set of samples and names so we can have fast look-up
        self.index(samples=(samples is not None), names=(names is not None))
        # Subset the samples
        if samples is not None:
            pts.samples = tuple(samp for samp in samples if samp in self._samp_idx)
            if len(pts.samples) < len(samples):
                diff = len(samples) - len(pts.samples)
                self.log.warning(
                    f"Saw {diff} fewer samples than requested. Proceeding with "
                    f"{len(pts.samples)} samples."
                )
            samp_idx = tuple(self._samp_idx[samp] for samp in pts.samples)
            if inplace:
                self._samp_idx = None
            pts.data = pts.data[samp_idx, :]
        # Subset the names
        if names is not None:
            pts.names = tuple(name for name in names if name in self._name_idx)
            if len(pts.names) < len(names):
                diff = len(names) - len(pts.names)
                self.log.warning(
                    f"Saw {diff} fewer names than requested. Proceeding with "
                    f"{len(pts.names)} names."
                )
            name_idx = tuple(self._name_idx[name] for name in pts.names)
            if inplace:
                self._name_idx = None
            pts.data = pts.data[:, name_idx]
        if not inplace:
            return pts
