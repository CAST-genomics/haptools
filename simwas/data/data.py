from __future__ import annotations
import numpy as np
from csv import reader
from pathlib import Path
from abc import ABC, abstractmethod


class Data(ABC):
    """
    Abstract class for accessing read-only data files

    Attributes
    ----------
    fname : Path
        The path to the read-only file containing the data
    data : np.array
        The contents of the data file, once loaded
    """

    def __init__(self, fname: Path):
        self.fname = fname
        self.data = None
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

    @abstractmethod
    def read(self):
        """
        Read the raw file contents into the class properties
        """
        if self.data is not None:
            raise AssertionError("The data has already been loaded.")
