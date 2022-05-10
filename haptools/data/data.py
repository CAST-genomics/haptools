from __future__ import annotations
from csv import reader
from pathlib import Path
from typing import Iterator
from collections import namedtuple
from abc import ABC, abstractmethod
from logging import getLogger, Logger

import numpy as np


class Data(ABC):
    """
    Abstract class for accessing read-only data files

    Attributes
    ----------
    fname : Path
        The path to the read-only file containing the data
    data : np.array
        The contents of the data file, once loaded
    log: Logger
        A logging instance for recording debug statements.
    """

    def __init__(self, fname: Path, log: Logger = None):
        self.fname = fname
        self.data = None
        self.log = log or getLogger(self.__class__.__name__)
        super().__init__()

    def __repr__(self):
        return str(self.fname)

    @classmethod
    @abstractmethod
    def load(cls: Data, fname: Path):
        """
        Read the file contents and perform any recommended pre-processing

        Parameters
        ----------
        fname : Path
            See documentation for :py:attr:`~.Data.fname`
        """
        pass

    def unset(self) -> bool:
        """
        Whether the data has been loaded into the object yet

        Returns
        -------
        bool
            True if :py:attr:`~.Data.data` is None else False
        """
        return self.data is None

    @abstractmethod
    def read(self):
        """
        Read the raw file contents into the class properties
        """
        if self.unset():
            self.log.warning("The data has already been loaded. Overriding.")

    @abstractmethod
    def __iter__(self) -> Iterator[namedtuple]:
        """
        Return an iterator over the raw file contents

        Yields
        ------
        Iterator[namedtuple]
            An iterator over each line in the file, where each line is encoded as a
            namedtuple containing each of the class properties
        """
        pass
