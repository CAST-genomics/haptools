from __future__ import annotations
import re
from pathlib import Path
from typing import Iterator
from collections import namedtuple

import numpy as np

from .data import Data
from .data import Variant


class Haplotype:
    """
    A haplotype within the tree

    Attributes
    ----------
    nodes : tuple[tuple[Variant, int]]
        An ordered collection of pairs, where each pair is a node and its allele
    data : npt.NDArray[np.bool_]
        A np array (with shape n x 2, num_samples x num_chromosomes) denoting the
        presence of this haplotype in each chromosome of each sample
    """

    # TODO: consider using a named tuple?
    nodes: tuple[tuple[Variant, int]]
    data: npt.NDArray[np.bool_]

    def __init__(
        self,
        nodes: tuple[tuple[Variant, int]] = tuple(),
        data: npt.NDArray[np.bool_] = None,
        num_samples: int = None,
    ):
        """
        Initialize an empty haplotype

        Parameters
        ----------
        nodes : tuple[tuple[Variant, int]]
            An ordered collection of pairs, where each pair is a node and its allele
        data : npt.NDArray[np.bool_]
            A np array (with length n x 2, num_samples x num_chromosomes) denoting the
            presence of this haplotype in each chromosome of each sample
        num_samples : int
            The number of samples in this haplotype
        """
        self.nodes = nodes
        if num_samples and data is None:
            self.data = np.ones((num_samples, 2), dtype=np.bool_)
        elif num_samples is None:
            self.data = data
        else:
            raise ValueError(
                "The data and num_samples arguments are mutually exclusive. Provide"
                " either one or the other."
            )

    def __repr__(self):
        return str(self.nodes)

    @classmethod
    def from_node(
        cls, node: Variant, allele: int, variant_genotypes: npt.NDArray[np.bool_]
    ) -> Haplotype:
        """
        Create a new haplotype with a single node entry

        Parameters
        ----------
        node : Variant
            The initializing node for this haplotype
        allele : int
            The allele associated with node
        variant_genotypes : npt.NDArray[np.bool_]
            A np array (with length n x 2, num_samples x num_chromosomes) denoting the
            presence of this haplotype in each chromosome of each sample

        Returns
        -------
        Haplotype
            The newly created haplotype object containing node and allele
        """
        return cls(((node, allele),), variant_genotypes)

    def append(
        self, node: Variant, allele: int, variant_genotypes: npt.NDArray[np.bool_]
    ) -> Haplotype:
        """
        Append a new node (variant) to this haplotype

        Parameters
        ----------
        node : Variant
            The node to add to this haplotype
        allele : int
            The allele associated with this node
        variant_genotypes : npt.NDArray[np.bool_]
            A np array (with length n x 2, num_samples x num_chromosomes) denoting the
            presence of this haplotype in each chromosome of each sample

        Returns
        -------
        Haplotype
            A new haplotype object extended by the node and its allele
        """
        new_haplotype = self.nodes + ((node, allele),)
        new_haplotype_values = self.data & variant_genotypes
        return Haplotype(new_haplotype, new_haplotype_values)

    # TODO: use @cached_property from python3.8, instead
    @property
    def node_indices(self) -> tuple[int]:
        """
        Get the indices of the nodes in this haplotype

        Returns
        -------
        tuple
            The indices of the nodes
        """
        return tuple(node[0].idx for node in self.nodes)

    def transform(self, genotypes: Genotypes, allele: int) -> npt.NDArray[np.bool_]:
        """
        Transform a genotypes matrix via the current haplotype:

        Each entry in the returned matrix denotes the presence of the current haplotype
        extended by each of the variants in the genotype matrix

        Parameters
        ----------
        genotypes : Genotypes
            The genotypes which to transform using the current haplotype
        allele : int
            The allele (either 0 or 1) of the SNPs we're adding

        Returns
        -------
        npt.NDArray[np.bool_]
            A 3D haplotype matrix similar to the genotype matrix but with haplotypes
            instead of variants in the columns. It will have the same shape except that
            the number of columns (second dimension) will have decreased by the number
            of variants in this haplotype.
        """
        # first, remove any variants that are already in this haplotype using np.delete
        # TODO: consider moving this outside of this function
        gens = np.delete(genotypes.data, self.node_indices, axis=1)
        # add extra axes to match shape of gens
        hap_data = self.data[:, np.newaxis]
        # use np.logical_and to superimpose the current haplotype onto the GT matrix
        return np.logical_and(gens == allele, hap_data)


class Haplotypes(Data):
    """
    A class for processing haplotypes from a file

    Attributes
    ----------
    fname : Path
        The path to the file containing the data
    data : list[dict]
        A list of dict describing the composition of a series of haplotypes

        Each haplotype dictionary is composed of these items:
            1) id (str): A haplotype ID
            2) chrom (str): The chromosome that this haplotype belongs to
            3) start (int): The start position of this haplotype
            4) end (int): The end position of this haplotype
            5) info (dict): Other information belonging to this haplotype. Some examples are:
                - tree (int): A tree ID
                - beta (float): The effect size of the haplotype-phenotype association
                - pval (float): The p-value of the haplotype-phenotype association
                - pip (float): A PIP from running SuSiE or some other tool
            6) variants (list[dict]): A list of dictionaries, one for each variant...


        Each variants dictionary is composed of these items:
            1) id (str): A variant ID
            2) pos (int): The start position of this variant
            3) allele (bool): The allele for this variant
            4) info (dict): Other information belonging to this variant. Some examples are:
                - score (float): The score of this variant within its haplotype
    info : dict
        A dict containing information about the haplotype and variant info fields
    format : dict
        A dictionary describing the types of lines in the file format spec as well as
        their format and data types
    version : str
        A string denoting the current file format version
    log: Logger
        A logging instance for recording debug statements.

    Examples
    --------
    >>> haplotypes = Haplotypes.load('tests/data/simple.hapss')
    """

    def __init__(self, fname: Path, log: Logger = None):
        super().__init__(fname, log)
        self.format = {
            "meta": {
                "id": "M",
                "val": ["version"],
                "fmt": ["s"],
            },
            "hap": {
                "id": "H",
                "val": ["id", "chrom", "start", "end"],
                "fmt": ["s", "s", "d", "d"],
            },
            "var": {
                "id": "V",
                "val": ["id", "hap", "start", "end", "allele"],
                "fmt": ["s", "s", "d", "d", "d"],
            },
        }
        self.version = "0.0.1"
        for val in self.format.keys():
            self.format[val]["str"] = self._create_fmt_str(self.format[val])
            self.format[val]["rgx"] = re.compile(self.format[val]["str"])
        self.data = []

    def _create_fmt_str(self, fmts):
        return (
            fmts["id"]
            + "\t"
            + "\t".join(
                [
                    "{" + val + ":" + fmt + "}"
                    for val, fmt in zip(fmts["val"], fmts["fmt"])
                ]
            )
            + "\n"
        )

    def _line_type(self, line: str, validate=False):
        """
        Return the type of line that this line matches

        Parameters
        ----------
        line : str
            A line of the .haps file
        validate : bool
            Whether to validate that the line's contents follow the spec
        """
        if validate:
            checker = lambda line: bool(self.format[line_type]["rgx"].match(line))
        else:
            checker = lambda line: self.format[line_type]["id"] == line[0]
        for line_type in self.format.keys():
            if checker(line):
                return self.format[line_type]["id"]
        # if none of the lines matched, return None
        return None


    def read(self, region: str = None, haplotypes: set[str] = None):
        """
        Read haplotypes from a .haps file into a list stored in :py:attr:`~.Haplotypes.data`

        Parameters
        ----------
        region : str, optional
            The region from which to extract haplotypes; ex: 'chr1:1234-34566' or 'chr7'

            For this to work, the .haps file must be indexed and the seqname must match!

            Defaults to loading all haplotypes
        haplotypes : list[str], optional
            A list of haplotype IDs corresponding to a subset of the haplotypes to
            extract

            Defaults to loading haplotypes from all samples

            Note that finding the haplotype lines can be slow if the region param is
            not specified
        """
        super().read()
        self.data = []
        # if the user requested a specific region or set of haplotypes, then we should
        # handle it using tabix
        # else, we use a regular text opener
        if region or haplotypes:
            if region:
                # split the region string so each portion is an element
                region = re.split(':|-', region)
                if len(region) > 1:
                    region[1:] = [int(pos) for pos in region[1:] if pos]
                # fetch region
                # we already know that each line will start with an H, so we don't
                # need to check that
                self.data = [
                    line.split("\t") # TODO
                    for line in pysam.TabixFile(self.fname).fetch(*region)
                ]
            else:
                for line in pysam.tabix_iterator(self.fname):
                    # we only want lines that start with an H
                    if self._line_type(line) == 'H':
                        # TODO: store the line
            # TODO: query for the variants
        else:
            # use hook_compressed to automatically handle gz files
            with hook_compressed(fname, mode="rt") as haps:
                hap_text = reader(haps, delimiter="\t")
                for line in hap_text:
                    # TODO: store the line
        if not region and not haplotypes:
            haps.close()

    def iterate(self, region: str = None, haplotypes: set[str] = None) -> Iterator[namedtuple]:
        """
        Read haplotypes from a .haps file line by line without storing anything

        Parameters
        ----------
        region : str, optional
            The region from which to extract haplotypes; ex: 'chr1:1234-34566' or 'chr7'

            For this to work, the .haps file must be indexed and the seqname must match!

            Defaults to loading all haplotypes
        haplotypes : list[str], optional
            A list of haplotype IDs corresponding to a subset of the haplotypes to
            extract

            Defaults to loading haplotypes from all samples

        Yields
        ------
        Iterator[namedtuple]
            An iterator over each line in the file, where each line is encoded as a
            namedtuple containing each of the class properties
        """
        with hook_compressed(self.fname, mode="rt") as haps:
            hap_text = reader(haps, delimiter="\t")
            pass

    def to_str(self) -> Generator[str, None, None]:
        """
        Create a string representation of this Haplotype

        Yields
        ------
        Generator[str, None, None]
            A list of lines (strings) to include in the output
        """
        yield self.format["meta"]["str"].format(version=self.version)
        for hap in self.data:
            yield self.format["hap"]["str"].format(**hap)
            for var in hap["variants"]:
                yield self.format["var"]["str"].format(**var)

    def __repr__(self):
        return "\n".join(self.to_str())

    def write(self, file: TextIO):
        """
        Write the contents of this Haplotypes object to the file given by fname

        Parameters
        ----------
        file : TextIO
            A file-like object to which this Haplotypes object should be written.
        """
        for line in self.to_str():
            file.write(line)
