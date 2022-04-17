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
    data : dict[dict]
        A dict of dict describing the composition of a series of haplotypes

        Each haplotype dictionary is keyed by an ID and composed of these items:
            1) chrom (str): The chromosome that this haplotype belongs to
            2) start (int): The start position of this haplotype
            3) end (int): The end position of this haplotype
            4) id (str): A haplotype ID; this is also the element key
            5) info (dict): Other information belonging to this haplotype. Some examples are:
                - tree (int): A tree ID
                - beta (float): The effect size of the haplotype-phenotype association
                - pval (float): The p-value of the haplotype-phenotype association
                - pip (float): A PIP from running SuSiE or some other tool
                - ancestry (string): A population code denoting the haplotype ancestry
            6) variants (dict[dict]): A dict of dictionaries, one for each variant...


        Each variants dictionary is keyed by an ID and composed of these items:
            1) pos (int): The start position of this variant
            2) id (str): A variant ID; this is also the element key
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
                "val": ["chrom", "start", "end", "id"],
                "fmt": ["s", "d", "d", "s"],
            },
            "var": {
                "id": "V",
                "val": ["hap", "start", "end", "id", "allele"],
                "fmt": ["s", "d", "d", "s", "d"],
            },
        }
        self.version = "0.0.1"
        for val in self.format.keys():
            self.format[val]["str"] = self._create_fmt_str(self.format[val])
            self.format[val]["rgx"] = re.compile(self.format[val]["str"])
        self.data = {}

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

    def _line_type(self, line: str, validate: bool = False):
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

    def _extract_haplotype(self, line: str) -> dict:
        """
        Convert a haplotype line into a dictionary for inclusion in self.data

        Paramaters
        ----------
        line : str
            A haplotype (H) line from the .haps file

        Returns
        -------
        A dictionary that can be included as a haplotype in self.data
        """
        hap = dict(zip(
            ('chrom', 'start', 'end', 'id', 'info'),
            line[2:].split("\t", maxsplit=4)
        ))
        # TODO: find a way to do this by using self.format?
        hap['start'] = int(hap['start'])
        hap['end'] = int(hap['end'])
        hap['info'] = dict(zip(self.info['hap'].keys(), hap['info'].split("\t")))
        return hap['id'], hap

    def _extract_variant(self, line: strs) -> dict:
        """
        Convert a variant line into a dictionary for inclusion in self.data

        Paramaters
        ----------
        line : str
            A variant (V) line from the .haps file

        Returns
        -------
        A dictionary that can be included as a variant in self.data
        """
        var = dict(zip(
            ('hap_id', 'pos', 'end', 'id', 'allele', 'info'),
            line[2:].split("\t", maxsplit=5)
        ))
        # TODO: find a way to do this by using self.format?
        hap_id = var['hap']
        var['pos'] = int(var['pos'])
        del var['hap'], var['end']
        var['allele'] = int(var['allele'])
        var['info'] = dict(zip(self.info['var'].keys(), var['info'].split("\t")))
        return var['id'], hap_id, var

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
        self.data = {}
        # if the user requested a specific region or set of haplotypes, then we should
        # handle it using tabix
        # else, we use a regular text opener
        if region or haplotypes:
            haps_file = pysam.TabixFile(self.fname)
            if region:
                # split the region string so each portion is an element
                region = re.split(':|-', region)
                if len(region) > 1:
                    region[1:] = [int(pos) for pos in region[1:] if pos]
                # fetch region
                # we already know that each line will start with an H, so we don't
                # need to check that
                for line in haps_file.fetch(region):
                    hap_id, hap = self._extract_haplotype(line)
                    if haplotypes is None or hap_id in haplotypes:
                        self.data[hap_id] = hap
            else:
                for line in haps_file.fetch():
                    # we only want lines that start with an H
                    line_type = self._line_type(line)
                    if line_type == 'H':
                        hap_id, hap = self._extract_haplotype(line)
                        if hap_id in haplotypes:
                            self.data[hap_id] = hap
                    elif line_type > 'H':
                        # if we've already passed all of the H's, we can just exit
                        # We assume the file has been sorted so that all of the H lines
                        # come before the V lines
                        break
            # query for the variants of each haplotype
            for hap_id in self.data:
                self.data[hap_id]['variants'] = {}
                # exclude variants outside the desired region
                hap_region = hap_id
                if region:
                    hap_region = hap_id + ":" + region.split(':', maxsplit=1)[1]
                # fetch region
                # we already know that each line will start with a V, so we don't
                # need to check that
                for variant in haps_file.fetch(*hap_region):
                    var_id, _, var = self._extract_variant(line)
                    self.data[hap_id]['variants'][var_id] = var
            haps_file.close()
        else:
            # use hook_compressed to automatically handle gz files
            with hook_compressed(fname, mode="rt") as haps:
                hap_text = reader(haps, delimiter="\t")
                for line in hap_text:
                    line_type = self._line_type(line)
                    if line_type == "H":
                        hap_id, hap = self._extract_haplotype(line)
                        self.data[hap_id].update(hap)
                    elif line_type == "V":
                        var_id, hap_id, var = self._extract_variant(line)
                        self.data.setdefault(hap_id, {})['variants'][var_id] = var

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
            HaplotypeRecord = namedtuple("HaplotypeRecord", "data info")
            VariantRecord = namedtuple("VariantRecord", "data info")
            for line in hap_text:
                line_type = self._line_type(line)
                if line_type == "H":
                    hap_id, hap = self._extract_haplotype(line)
                    yield HaplotypeRecord(hap, self.info)
                elif line_type == "V":
                    var_id, hap_id, var = self._extract_variant(line)
                    # add the haplotype back in, since otherwise, the user won't know
                    # which haplotype this variant belongs to
                    var['hap_id'] = hap_id
                    yield VariantRecord(var, self.info)

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
