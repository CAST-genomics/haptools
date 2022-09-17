from __future__ import annotations
from csv import reader
from pathlib import Path
from typing import NewType
from collections import namedtuple
from collections.abc import Iterable
from logging import getLogger, Logger
from fileinput import hook_compressed

import numpy as np
import numpy.typing as npt

from .data import Data


# A haplotype block consists of
# 1) A population label (str), like 'YRI'
# 2) The end position of the block in bp (int), like 1001038
# 3) The end position of the block in cM (int), like 43.078
HapBlock = namedtuple("HapBlock", "pop bp cm")

# This dict maps chroms (as strings) to a tuple of two lists, one for each chromosome
# TODO: consider using IntervalTrees instead of simple lists:
#     https://github.com/chaimleib/intervaltree
# Let's define a type alias, "SampleBlocks", for future use...
SampleBlocks = NewType(
    "SampleBlocks", "dict[str, tuple[list[HapBlock], list[HapBlock]]]"
)

class Breakpoints(Data):
    """
    A class for processing breakpoints from a file

    Attributes
    ----------
    data : dict[str, SampleBlocks]
        The haplotype blocks for each chromosome in each sample
        This dict maps samples (as strings) to their haplotype blocks (as SampleBlocks)
    fname : Path | str
        The path to the file containing the data
    log: Logger
        A logging instance for recording debug statements.

    Examples
    --------
    >>> breakpoints = Breakpoints.load('tests/data/test.bp')
    """

    def __init__(self, fname: Path | str, log: Logger = None):
        super().__init__(fname, log)
        self._ext = "bp"

    @classmethod
    def load(
        cls: Breakpoints, fname: Path | str, samples: set[str] = None
    ) -> Breakpoints:
        """
        Load breakpoints from a TSV file

        Read the file contents and standardize the Breakpoints

        Parameters
        ----------
        fname
            See documentation for :py:attr:`~.Data.fname`
        samples : set[str], optional
            See documentation for :py:meth:`~.Breakpoints.read`

        Returns
        -------
        Breakpoints
            A Breakpoints object with the data loaded into its properties
        """
        breakpoints = cls(fname)
        breakpoints.read(samples)
        return breakpoints

    def read(self, samples: set[str] = None):
        """
        Read phenotypes from a TSV file into a numpy matrix stored in :py:attr:`~.Penotypes.data`

        Parameters
        ----------
        samples : set[str], optional
            A subset of the samples for which to extract breakpoints

            Defaults to loading breakpoints for all samples

        Raises
        ------
        AssertionError
            If the provided file doesn't follow the expected format
        """
        super().read()
        # call self.__iter__() to obtain the data
        self.data = dict(self.__iter__(samples))
        self.log.info(f"Loaded {len(self.data)} samples from .{self._ext} file")

    def __iter__(self, samples: set[str] = None) -> Iterable[str, SampleBlocks]:
        """
        Read phenotypes from a TSV line by line without storing more than a single
        sample at a time

        Parameters
        ----------
        samples : set[str], optional
            See documentation for :py:meth:`~.Breakpoints.read`

        Returns
        ------
        Iterable[str, SampleBlocks]
            An iterator over each sample in the file, where the sample if specified
            first as a string, and then followed by its SampleBlocks
        """
        bps = hook_compressed(self.fname, mode="rt")
        bp_text = reader(bps, delimiter="\t")
        samp = None
        blocks = {}
        for line in bp_text:
            # ignore all of the comment lines
            if line[0].startswith("#"):
                continue
            # check: is this a sample name or a hap block line?
            if len(line) == 1:
                line = line[0]
                # get the end of the sample name (ex: "_1" or "_2")
                strand_num = line.rsplit("_", 1)[-1]
                if not (strand_num in ("1", "2")):
                    self.log.warning(
                        f"Ignoring improperly formatted line in bp file: '{line}'"
                    )
                    continue
                strand_num = int(strand_num) - 1
                if strand_num:
                    if samp is None:
                        self.log.error(f"bp file does not start with first strand")
                    elif samp != line[:-2]:
                        self.log.error(
                            f"'{samp}_1' was not followed by '{samp}_2' "
                            "in the bp file"
                        )
                else:
                    if samp is not None and (samples is None or samp in samples):
                        # output the previous sample
                        yield samp, blocks
                    samp = line[:-2]
                    blocks = {}
            elif len(line) == 4:
                pop, chrom, pos, cm = line
                block = blocks.setdefault(chrom, ([], []))[strand_num]
                block.append(HapBlock(pop, int(pos), float(cm)))
            else:
                self.log.warning(f"Ignoring improperly formatted line in bp file: '{line}'")
        if samp is not None and (samples is None or samp in samples):
            # output the previous sample
            yield samp, blocks
        bps.close()

    def population_array(
        self,
        variants: np.array,
        samples: tuple[str] = None,
    ) -> tuple[dict[str, int], npt.NDArray[np.uint8]]:
        """
        Output an array denoting the population labels of each variant for each sample

        Parameters
        ----------
        variants : np.array
            Variant-level meta information in a mixed np array of dtypes:
            CHROM (str) and POS (int)
        samples : tuple[str], optional
            A subset of samples to include in the output, ordered by their given order

        Returns
        ------
        dict[str, int]
            A dict mapping population label strings to unique integers
        npt.NDArray[np.uint8]
            An array of shape: samples x variants x 2

            The array is composed of unique integers, where each integer encodes a
            population label in the returned dict
        """
        labels = {}
        label_ct = 0
        if samples is None:
            data = self.data
        else:
            data = {samp: self.data[samp] for samp in samples}
        arr = np.empty((len(data), len(variants), 2), dtype=np.uint8)
        # Note: Despite the fact that this code has four nested for-loops, it is still
        # 1-2 orders of magnitude faster than trying to load this array from a VCF
        for samp_idx, blocks in enumerate(data.values()):
            for var_idx, variant in enumerate(variants):
                chrom, pos = variant
                for strand_num in range(2):
                    # try to figure out the right pop label by iterating through them
                    # all in order from smallest bp to largest bp
                    pop = blocks[chrom][strand_num][0].pop
                    for block in blocks[chrom][strand_num]:
                        if block.bp > pos:
                            break
                        pop = block.pop
                    # obtain the proper pop label number
                    if pop not in labels:
                        labels[pop] = label_ct
                        label_ct += 1
                    arr[samp_idx, var_idx, strand_num] = labels[pop]
        return labels, arr

    def write(self):
        raise ValueError("Not Implemented")
