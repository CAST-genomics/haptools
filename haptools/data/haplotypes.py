from __future__ import annotations
import re
from pathlib import Path
from logging import getLogger, Logger
from dataclasses import dataclass, field, fields
from typing import Iterator, get_type_hints, Generator

import numpy as np

from .data import Data


# We declare this class to be a dataclass to automatically define __init__ and a few
# other methods.
@dataclass
class Variant:
    """
    A variant within the .hap format spec

    In order to use this class with the Haplotypes class, you should
    1) add properties to the class for each of extra fields
    2) override the extras() method
    3) change the default value of the fmt property to describe the line

    Attributes
    ----------
    start: int
        The chromosomal start position of the variant
    end: int
        The chromosomal end position of the variant

        In most cases this will be the same as the start position
    id: str
        The variant's unique ID
    allele: str
        The allele of this variant within the Haplotype
    fmt: str
        A format string describing the corresponding line from the .haps format spec

    Examples
    --------
    Let's extend this class and add an extra field called "score"

    >>> from dataclasses import dataclass, field
    >>> @dataclass
    >>> class CustomVariant(Variant):
    ...     score: float
    ...     fmt: str = field(default="V\\t{hap:s}\t{start:d}\t{end:d}\t{id:s}\t{allele:s}\t{score:.3f}", init=False)
    ...
    ...     @staticmethod
    ...     def extras() -> tuple:
    ...         return (
    ...             "#V\tscore\tfloat\tImportance of inclusion",
    ...         )
    """

    start: int
    end: int
    id: str
    allele: str
    fmt: str = field(
        default="V\t{hap:s}\t{start:d}\t{end:d}\t{id:s}\t{allele:s}", init=False
    )

    @property
    def ID(self):
        """
        Create an alias for the id property
        """
        return self.id

    @classmethod
    def from_hap_spec(cls: Variant, line: str) -> tuple[str, Variant]:
        """
        Convert a variant line into a Variant object in the .hap format spec

        Paramaters
        ----------
        line: str
            A variant (V) line from the .hap file

        Returns
        -------
        tuple[str, Variant]
            The haplotype ID and Variant object for the variant
        """
        line = line[2:].split("\t")
        hap_id = line[0]
        var_fields = {}
        idx = 1
        for name, val in get_type_hints(cls).items():
            if name != "fmt":
                var_fields[name] = val(line[idx])
                idx += 1
        return hap_id, cls(**var_fields)

    def to_hap_spec(self, hap_id: str) -> str:
        """
        Convert a Variant object into a variant line in the .hap format spec

        Parameters
        ----------
        hap_id: str
            The ID of the haplotype associated with this variant

        Returns
        -------
        str
            A valid variant line (V) in the .hap format spec
        """
        return fmt.format(**self.__dict__, hap=hap_id)

    @staticmethod
    def extras() -> tuple:
        """
        Return the header lines of the extra fields that are supported

        Returns
        -------
        tuple
            The header lines of the extra fields
        """
        return tuple()


# We declare this class to be a dataclass to automatically define __init__ and a few
# other methods.
@dataclass
class Haplotype:
    """
    A haplotype within the .hap format spec

    In order to use this class with the Haplotypes class, you should
    1) add properties to the class for each of extra fields
    2) override the extras() method
    3) change the default value of the fmt property to describe the line

    Attributes
    ----------
    chrom: str
        The contig to which this haplotype belongs
    start: int
        The chromosomal start position of the haplotype
    end: int
        The chromosomal end position of the haplotype
    id: str
        The haplotype's unique ID
    variants: list[Variant]
        A list of the variants in this haplotype
    fmt: str
        A format string describing the corresponding line from the .haps format spec

    Examples
    --------
    Let's extend this class and add an extra field called "ancestry"

    >>> from dataclasses import dataclass, field
    >>> @dataclass
    >>> class CustomHaplotype(Haplotype):
    ...     ancestry: str
    ...     fmt: str = field(default="H\\t{chrom:s}\\t{start:d}\\t{end:d}\\t{id:s}\\t{ancestry:s}", init=False)
    ...
    ...     @staticmethod
    ...     def extras() -> tuple:
    ...         return (
    ...             "#H\\tancestry\\tstr\\tLocal ancestry",
    ...         )
    """

    chrom: str
    start: int
    end: int
    id: str
    variants: tuple = field(default_factory=tuple, init=False)
    fmt: str = field(default="H\t{chrom:s}\t{start:d}\t{end:d}\t{id:s}", init=False)

    @property
    def ID(self):
        """
        Create an alias for the id property
        """
        return self.id

    @classmethod
    def from_hap_spec(
        cls: Haplotype, line: str, variants: tuple = tuple()
    ) -> Haplotype:
        """
        Convert a variant line into a Haplotype object in the .hap format spec

        Paramaters
        ----------
        line: str
            A variant (V) line from the .hap file

        Returns
        -------
        Haplotype
            The Haplotype object for the variant
        """
        line = line[2:].split("\t")
        hap_fields = {}
        idx = 0
        for name, val in get_type_hints(cls).items():
            if name not in ("fmt", "variants"):
                hap_fields[name] = val(line[idx])
                idx += 1
        hap = cls(**hap_fields)
        hap.variants = variants
        return hap

    def to_hap_spec(self) -> str:
        """
        Convert a Haplotype object into a haplotype line in the .hap format spec

        Returns
        -------
        str
            A valid haplotype line (H) in the .hap format spec
        """
        return self.fmt.format(**self.__dict__)

    @staticmethod
    def extras() -> tuple:
        """
        Return the header lines of the extra fields that are supported

        Returns
        -------
        tuple
            The header lines of the extra fields
        """
        return tuple()


class Haplotypes(Data):
    """
    A class for processing haplotypes from a file

    Attributes
    ----------
    fname: Path
        The path to the file containing the data
    data: dict[str, Haplotype]
        A dict of Haplotype objects keyed by their IDs
    types: dict
        A dict of class names keyed by the symbol denoting their line type

        Ex: {'H': Haplotype, 'V': Variant}
    version: str
        A string denoting the current file format version
    log: Logger
        A logging instance for recording debug statements.

    Examples
    --------
    >>> haplotypes = Haplotypes.load('tests/data/example.hap')

    >>> haplotypes = Haplotypes('tests/data/example.hap.gz', Haplotype, Variant)
    >>> haplotypes.read()
    >>> haps = haplotypes.data
    """

    def __init__(
        self,
        fname: Path,
        haplotype: type[Haplotype] = Haplotype,
        variant: type[Variant] = Variant,
        log: Logger = None,
    ):
        super().__init__(fname, log)
        self.data = {}
        self.types = {"H": haplotype, "V": variant}
        self.version = "0.0.1"

    @classmethod
    def load(
        cls: Haplotypes, fname: Path, region: str = None, haplotypes: set[str] = None
    ) -> Haplotypes:
        """
        Load haplotypes from a .hap file

        Read the file contents

        Parameters
        ----------
        fname: Path
            See documentation for :py:attr:`~.Data.fname`
        region: str, optional
            See documentation for :py:meth:`~.Haplotypes.read`
        haplotypes: list[str], optional
            See documentation for :py:meth:`~.Haplotypes.read`

        Returns
        -------
        Haplotypes
            A Haplotypes object with the data loaded into its properties
        """
        haps = cls(fname)
        haps.read(region, haplotypes)
        return haps

    def check_header(self, lines: list[str], check_version=False):
        """
        Check 1) that the version number matches and 2) that extra fields declared in
        # the .haps file can be handled by the the Variant and Haplotype classes
        # provided in __init__()

        Parameters
        ----------
        lines: list[str]
            Header lines from the .hap file
        check_version: bool = False
            Whether to also check the version of the file

        Raises
        ------
        ValueError
            If any of the header lines are not supported
        """
        if check_version:
            version_line = lines[0].split("\t")
            assert version_line[1] == "version"
            if version_line[2] != self.version:
                self.log.warning(
                    f"The version of the provided .hap file is {version_line} but this"
                    f" tool expected {self.version}"
                )
        for line in lines:
            if (
                line[1] in self.types.keys()
                and line not in self.types[line[1]].extras()
            ):
                # extract the name of the extra field
                name = line.split("\t", maxsplit=1)[1]
                raise ValueError(
                    f"The extra field '{name}' is declared in the header of the .hap"
                    " filebut is not accepted by this tool."
                )

    def _line_type(self, line: str) -> type:
        """
        Return the type of line that this line matches

        Parameters
        ----------
        line: str
            A line of the .hap file

        Returns
        -------
        type
            The name of the class corresponding with the type of this line
        """
        line_types = self.types.keys()
        if line[0] in line_types:
            return line[0]
        else:
            # if none of the lines matched, return None
            return None

    def read(self, region: str = None, haplotypes: set[str] = None):
        """
        Read haplotypes from a .hap file into a list stored in :py:attr:`~.Haplotypes.data`

        Parameters
        ----------
        region: str, optional
            The region from which to extract haplotypes; ex: 'chr1:1234-34566' or 'chr7'

            For this to work, the .hap file must be indexed and the seqname must match!

            Defaults to loading all haplotypes
        haplotypes: list[str], optional
            A list of haplotype IDs corresponding to a subset of the haplotypes to
            extract

            Defaults to loading haplotypes from all samples

            Note that finding the haplotype lines can be slow if the region param is
            not specified
        """
        super().read()
        self.data = {}
        # if the user requested a specific region or set of haplotypes, then we should
        # handle it using tabix
        # else, we use a regular text opener
        if region or haplotypes:
            haps_file = pysam.TabixFile(self.fname)
            self.check_header(list(haps_file.header))
            if region:
                # split the region string so each portion is an element
                region = re.split(":|-", region)
                if len(region) > 1:
                    region[1:] = [int(pos) for pos in region[1:] if pos]
                # fetch region
                # we already know that each line will start with an H, so we don't
                # need to check that
                for line in haps_file.fetch(region):
                    hap = self.types["H"].from_hap_spec(line)
                    if haplotypes is None or hap.id in haplotypes:
                        self.data[hap.id] = hap
            else:
                for line in haps_file.fetch():
                    # we only want lines that start with an H
                    line_type = self._line_type(line)
                    if line_type == "H":
                        hap = self.types["H"].from_hap_spec(line)
                        if hap.id in haplotypes:
                            self.data[hap.id] = hap
                    elif line_type > "H":
                        # if we've already passed all of the H's, we can just exit
                        # We assume the file has been sorted so that all of the H lines
                        # come before the V lines
                        break
            # query for the variants of each haplotype
            for hap_id in self.data:
                self.data[hap_id]["variants"] = {}
                # exclude variants outside the desired region
                hap_region = hap_id
                if region:
                    hap_region = hap_id + ":" + region.split(":", maxsplit=1)[1]
                # fetch region
                # we already know that each line will start with a V, so we don't
                # need to check that
                for variant in haps_file.fetch(*hap_region):
                    var = self.types["V"].from_hap_spec(line)[1]
                    self.data[hap_id].variants.append(var)
            haps_file.close()
        else:
            # the file is not indexed, so we can't assume it's sorted, either
            # use hook_compressed to automatically handle gz files
            with hook_compressed(fname, mode="rt") as haps:
                haps = {}
                for line in haps:
                    line_type = self._line_type(line)
                    header_lines = []
                    if line_type == "#":
                        # store header for later
                        try:
                            header_lines.append(line)
                        except AttributeError:
                            # this happens when we encounter a line beginning with a #
                            # after already having seen an H or V line
                            # in this case, it's usually just a comment, so we can ignore
                            pass
                    else:
                        if header_lines:
                            self.check_header(header_lines)
                            header_lines = None
                        if line_type == "H":
                            hap = self.types["H"].from_hap_spec(line)
                            self.data[hap.id] = hap
                        elif line_type == "V":
                            hap_id, var = self.types["V"].from_hap_spec(line)
                            haps.set_default(hap_id, []).append(var)
                        else:
                            self.log.warning(
                                f"Ignoring unsupported line type '{line[0]}'"
                            )
                for hap in haps:
                    self.data[hap].variants = tuple(haps[hap])

    def iterate(self) -> Iterator[Variant | Haplotype]:
        """
        Read haplotypes from a .hap file line by line without storing anything

        Yields
        ------
        Iterator[Variant|Haplotype]
            An iterator over each line in the file, where each line is encoded as a
            Variant or Haplotype containing each of the class properties
        """
        with hook_compressed(self.fname, mode="rt") as haps:
            # note that we do not assume the file is indexed or sorted here
            # it doesn't really matter anyway, since we aren't storing anything
            hap_text = reader(haps, delimiter="\t")
            for line in hap_text:
                line_type = self._line_type(line)
                header_lines = []
                if line_type == "#":
                    # store header for later
                    try:
                        header_lines.append(line)
                    except AttributeError:
                        # this happens when we encounter a line beginning with a #
                        # after already having seen an H or V line
                        # in this case, it's usually just a comment, so we can ignore
                        pass
                else:
                    if header_lines:
                        self.check_header(header_lines)
                        header_lines = None
                    if line_type == "H":
                        hap = self.types["H"].from_hap_spec(line)
                        yield hap
                    elif line_type == "V":
                        hap_id, var = self.types["V"].from_hap_spec(line)
                        # add the haplotype, since otherwise, the user won't know
                        # which haplotype this variant belongs to
                        var.hap = hap_id
                        yield var

    def to_str(self) -> Generator[str, None, None]:
        """
        Create a string representation of this Haplotype

        Yields
        ------
        Generator[str, None, None]
            A list of lines (strings) to include in the output
        """
        yield "# version " + self.version
        yield from Haplotype.extras()
        yield from Variant.extras()
        for hap in self.data:
            yield self.types["H"].to_hap_spec(hap)
        for hap in self.data:
            for var in hap.variants:
                yield self.types["V"].to_hap_spec(var, hap.id)

    def __repr__(self):
        return "\n".join(self.to_str())

    def write(self):
        """
        Write the contents of this Haplotypes object to the file given by fname

        Parameters
        ----------
        file: TextIO
            A file-like object to which this Haplotypes object should be written.
        """
        with hook_compressed(fname, mode="wt") as haps:
            for line in self.to_str():
                haps.write(line + "\n")
