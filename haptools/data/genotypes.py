from __future__ import annotations
import re
from csv import reader
from pathlib import Path
from typing import Iterator
from collections import namedtuple
from logging import getLogger, Logger
from fileinput import hook_compressed

import numpy as np
import numpy.typing as npt
from cyvcf2 import VCF, Variant
from pysam import VariantFile, TabixFile

from .data import Data


class Genotypes(Data):
    """
    A class for processing genotypes from a file

    Attributes
    ----------
    data : npt.NDArray
        The genotypes in an n (samples) x p (variants) x 2 (strands) array
    fname : Path
        The path to the read-only file containing the data
    samples : tuple[str]
        The names of each of the n samples
    variants : np.array
        Variant-level meta information:
            1. ID
            2. CHROM
            3. POS
            4. AAF: allele freq of alternate allele (or MAF if to_MAC() is called)
    log: Logger
        A logging instance for recording debug statements
    _prephased : bool
        If True, assume that the genotypes are phased. Otherwise, extract their phase
        when reading from the VCF.
    _samp_idx : dict[str, int]
        Sample index; maps samples to indices in self.samples
    _var_idx : dict[str, int]
        Variant index; maps variant IDs to indices in self.variants

    Examples
    --------
    >>> genotypes = Genotypes.load('tests/data/simple.vcf')
    >>> # directly access the loaded variants, samples, and genotypes (in data)
    >>> genotypes.variants
    >>> genotypes.samples
    >>> genotypes.data
    """

    def __init__(self, fname: Path, log: Logger = None):
        super().__init__(fname, log)
        self.samples = tuple()
        self.variants = np.array(
            [],
            dtype=[
                ("id", "U50"),
                ("chrom", "U10"),
                ("pos", np.uint32),
                ("aaf", np.float64),
            ],
        )
        self._prephased = False
        self._samp_idx = None
        self._var_idx = None

    @classmethod
    def load(
        cls: Genotypes,
        fname: Path,
        region: str = None,
        samples: list[str] = None,
        variants: set[str] = None,
    ) -> Genotypes:
        """
        Load genotypes from a VCF file

        Read the file contents, check the genotype phase, and create the MAC matrix

        Parameters
        ----------
        fname
            See documentation for :py:attr:`~.Data.fname`
        region : str, optional
            See documentation for :py:meth:`~.Genotypes.read`
        samples : list[str], optional
            See documentation for :py:meth:`~.Genotypes.read`
        variants : set[str], optional
            See documentation for :py:meth:`~.Genotypes.read`

        Returns
        -------
        genotypes
            A Genotypes object with the data loaded into its properties
        """
        genotypes = cls(fname)
        genotypes.read(region, samples, variants)
        genotypes.check_missing()
        genotypes.check_biallelic()
        genotypes.check_phase()
        # genotypes.to_MAC()
        return genotypes

    def read(
        self,
        region: str = None,
        samples: list[str] = None,
        variants: set[str] = None,
        max_variants: int = None,
    ):
        """
        Read genotypes from a VCF into a numpy matrix stored in :py:attr:`~.Genotypes.data`

        Raises
        ------
        ValueError
            If the genotypes array is empty

        Parameters
        ----------
        region : str, optional
            The region from which to extract genotypes; ex: 'chr1:1234-34566' or 'chr7'

            For this to work, the VCF must be indexed and the seqname must match!

            Defaults to loading all genotypes
        samples : list[str], optional
            A subset of the samples from which to extract genotypes

            Defaults to loading genotypes from all samples
        variants : set[str], optional
            A set of variant IDs for which to extract genotypes

            All other variants will be ignored. This may be useful if you're running
            out of memory.
        max_variants : int, optional
            The maximum mumber of variants to load from the file. Setting this value
            helps preallocate the arrays, making the process faster and less memory
            intensive. You should use this option if your processes are frequently
            "Killed" from memory overuse.

            If you don't know how many variants there are, set this to a large number
            greater than what you would except. The np array will be resized
            appropriately. You can also use the bcftools "counts" plugin to obtain the
            number of expected sites within a region.

            Note that this value is ignored if the variants argument is provided.
        """
        super().read()
        records = self.__iter__(region=region, samples=samples, variants=variants)
        if variants is not None:
            max_variants = len(variants)
        # check whether we can preallocate memory instead of making copies
        if max_variants is None:
            self.log.warning(
                "The max_variants parameter was not specified. We have no choice but to"
                " append to an ever-growing array, which can lead to memory overuse!"
            )
            variants_arr = []
            data_arr = []
            for rec in records:
                variants_arr.append(rec.variants)
                data_arr.append(rec.data)
            self.log.info(f"Copying {len(variants_arr)} variants into np arrays.")
            # convert to np array for speedy operations later on
            self.variants = np.array(variants_arr, dtype=self.variants.dtype)
            self.data = np.array(data_arr, dtype=np.uint8)
        else:
            # preallocate arrays! this will save us lots of memory and speed b/c
            # appends can sometimes make copies
            self.variants = np.empty((max_variants,), dtype=self.variants.dtype)
            # in order to check_phase() later, we must store the phase info, as well
            self.data = np.empty(
                (max_variants, len(self.samples), (2 + (not self._prephased))),
                dtype=np.uint8,
            )
            num_seen = 0
            for rec in records:
                self.variants[num_seen] = rec.variants
                self.data[num_seen] = rec.data
                num_seen += 1
            if max_variants > num_seen:
                self.log.info(
                    f"Removing {max_variants-num_seen} unneeded variant records that "
                    "were preallocated b/c max_variants was specified."
                )
                self.variants = self.variants[:num_seen]
                self.data = self.data[:num_seen]
        if self.data.shape == (0, 0, 0):
            self.log.warning(
                "Failed to load genotypes. If you specified a region, check that the"
                " contig name matches! For example, double-check the 'chr' prefix."
            )
        # transpose the GT matrix so that samples are rows and variants are columns
        self.log.info("Transposing genotype matrix.")
        self.data = self.data.transpose((1, 0, 2))

    def _variant_arr(self, record: Variant):
        """
        Construct a np array from the metadata in a line of the VCF

        This is a helper function for :py:meth:`~.Genotypes._iterate`. It's separate
        so that it can easily be overridden in any child classes.

        Parameters
        ----------
        record: Variant
            A cyvcf2.Variant object from which to fetch metadata

        Returns
        -------
        npt.NDArray
            A row from the :py:attr:`~.Genotypes.variants` array
        """
        return np.array(
            (record.ID, record.CHROM, record.POS, record.aaf),
            dtype=self.variants.dtype,
        )

    def _iterate(self, vcf: VCF, region: str = None, variants: set[str] = None):
        """
        A generator over the lines of a VCF

        This is a helper function for :py:meth:`~.Genotypes.__iter__`

        Parameters
        ----------
        vcf: VCF
            The cyvcf2.VCF object from which to fetch variant records
        region : str, optional
            See documentation for :py:meth:`~.Genotypes.read`
        variants : set[str], optional
            See documentation for :py:meth:`~.Genotypes.read`

        Yields
        ------
        Iterator[namedtuple]
            An iterator over each line in the file, where each line is encoded as a
            namedtuple containing each of the class properties
        """
        self.log.info(f"Loading genotypes from {len(self.samples)} samples")
        Record = namedtuple("Record", "data variants")
        # iterate over each line in the VCF
        # note, this can take a lot of time if there are many samples
        for variant in vcf(region):
            if variants is not None and variant.ID not in variants:
                continue
            # save meta information about each variant
            variant_arr = self._variant_arr(variant)
            # extract the genotypes to a matrix of size 1 x p x 3
            # the last dimension has three items:
            # 1) presence of REF in strand one
            # 2) presence of REF in strand two
            # 3) whether the genotype is phased (if self._prephased is False)
            data = np.array(variant.genotypes, dtype=np.uint8)
            data = data[:, : (2 + (not self._prephased))]
            yield Record(data, variant_arr)
        vcf.close()

    def __iter__(
        self, region: str = None, samples: list[str] = None, variants: set[str] = None
    ) -> Iterator[namedtuple]:
        """
        Read genotypes from a VCF line by line without storing anything

        Parameters
        ----------
        region : str, optional
            See documentation for :py:meth:`~.Genotypes.read`
        samples : list[str], optional
            See documentation for :py:meth:`~.Genotypes.read`
        variants : set[str], optional
            See documentation for :py:meth:`~.Genotypes.read`

        Returns
        -------
        Iterator[namedtuple]
            See documentation for :py:meth:`~.Genotypes._iterate`
        """
        vcf = VCF(str(self.fname), samples=samples, lazy=True)
        self.samples = tuple(vcf.samples)
        # call another function to force the lines above to be run immediately
        # see https://stackoverflow.com/a/36726497
        return self._iterate(vcf, region, variants)

    def index(self, samples: bool = True, variants: bool = True):
        """
        Call this function once to improve the amortized time-complexity of look-ups of
        samples and variants by their ID. This is useful if you intend to later subset
        by a set of samples or variant IDs.
        The time complexity of this function should be roughly O(n+m) if both
        parameters are True. Otherwise, it will be either O(n) or O(m).

        Parameters
        ----------
        samples: bool, optional
            Whether to index the samples for fast loop-up. Adds complexity O(n).
        variants: bool, optional
            Whether to index the variants for fast look-up. Adds complexity O(m).
        """
        if samples and self._samp_idx is None:
            self._samp_idx = dict(zip(self.samples, range(len(self.samples))))
        if variants and self._var_idx is None:
            self._var_idx = dict(zip(self.variants["id"], range(len(self.variants))))

    def subset(
        self,
        samples: tuple[str] = None,
        variants: tuple[str] = None,
        inplace: bool = False,
    ):
        """
        Subset these genotypes to a smaller set of samples or a smaller set of variants

        The order of the samples and variants in the subsetted instance will match
        the order in the provided tuple parameters.

        Parameters
        ----------
        samples: tuple[str]
            A subset of samples to keep
        variants: tuple[str]
            A subset of variant IDs to keep
        inplace: bool, optional
            If False, return a new Genotypes object; otherwise, alter the current one

        Returns
        -------
            A new Genotypes object if inplace is set to False, else returns None
        """
        # First, initialize variables
        gts = self
        if not inplace:
            gts = self.__class__(self.fname, self.log)
        gts.samples = self.samples
        gts.variants = self.variants
        gts.data = self.data
        # Index the current set of samples and variants so we can have fast look-up
        self.index(samples=(samples is not None), variants=(variants is not None))
        # Subset the samples
        if samples is not None:
            gts.samples = tuple(samp for samp in samples if samp in self._samp_idx)
            if len(gts.samples) < len(samples):
                diff = len(samples) - len(gts.samples)
                self.log.warning(
                    f"Saw {diff} fewer samples than requested. Proceeding with "
                    f"{len(gts.samples)} samples."
                )
            samp_idx = tuple(self._samp_idx[samp] for samp in gts.samples)
            if inplace:
                self._samp_idx = None
            gts.data = gts.data[samp_idx, :]
        # Subset the variants
        if variants is not None:
            var_idx = [self._var_idx[var] for var in variants if var in self._var_idx]
            if len(var_idx) < len(variants):
                diff = len(variants) - len(var_idx)
                self.log.warning(
                    f"Saw {diff} fewer variants than requested. Proceeding with "
                    f"{len(var_idx)} variants."
                )
            gts.variants = self.variants[var_idx]
            if inplace:
                self._var_idx = None
            gts.data = gts.data[:, var_idx]
        if not inplace:
            return gts

    def check_missing(self, discard_also=False):
        """
        Check that each sample is properly genotyped

        Raises
        ------
        ValueError
            If any of the samples have missing genotypes 'GT: .|.'

        Parameters
        ----------
        discard_also : bool, optional
            If True, discard any samples that are missing genotypes without raising a
            ValueError
        """
        # check: are there any samples that have genotype values that are empty?
        # A genotype value equal to the max for uint8 indicates the value was missing
        missing = np.any(self.data[:, :, :2] == np.iinfo(np.uint8).max, axis=2)
        if np.any(missing):
            samp_idx, variant_idx = np.nonzero(missing)
            if discard_also:
                self.log.info(
                    f"Ignoring missing genotypes from {len(samp_idx)} samples"
                )
                self.data = np.delete(self.data, samp_idx, axis=0)
                self.samples = tuple(np.delete(self.samples, samp_idx))
                self._samp_idx = None
            else:
                raise ValueError(
                    "Genotype with ID {} at POS {}:{} is missing for sample {}".format(
                        *tuple(self.variants[variant_idx[0]])[:3],
                        self.samples[samp_idx[0]],
                    )
                )
        if discard_also and not self.data.shape[0]:
            self.log.warning(
                "All samples were discarded! Check that that none of your variants are"
                " missing genotypes (GT: '.|.')."
            )

    def check_biallelic(self, discard_also=False):
        """
        Check that each genotype is composed of only two alleles

        This function modifies the dtype of :py:attr:`~.Genotypes.data` from uint8 to bool

        Raises
        ------
        ValueError
            If any of the genotypes have more than two alleles

        Parameters
        ----------
        discard_also : bool, optional
            If True, discard any multiallelic variants without raising a ValueError
        """
        if self.data.dtype == np.bool_:
            self.log.warning("All genotypes are already biallelic")
            return
        # check: are there any variants that have genotype values above 1?
        # A genotype value above 1 would imply the variant has more than one ALT allele
        multiallelic = np.any(self.data[:, :, :2] > 1, axis=2)
        if np.any(multiallelic):
            samp_idx, variant_idx = np.nonzero(multiallelic)
            if discard_also:
                self.log.info(f"Ignoring {len(variant_idx)} multiallelic variants")
                self.data = np.delete(self.data, variant_idx, axis=1)
                self.variants = np.delete(self.variants, variant_idx)
                self._var_idx = None
            else:
                raise ValueError(
                    "Variant with ID {} at POS {}:{} is multiallelic for sample {}"
                    .format(
                        *tuple(self.variants[variant_idx[0]])[:3],
                        self.samples[samp_idx[0]],
                    )
                )
        if discard_also and not self.data.shape[1]:
            self.log.warning(
                "All variants were discarded! Check that there are biallelic variants "
                "in your dataset."
            )
        self.data = self.data.astype(np.bool_)

    def check_phase(self):
        """
        Check that the genotypes are phased then remove the phasing info from the data

        This function modifies :py:attr:`~.Genotypes.data` in-place

        Raises
        ------
        ValueError
            If any heterozgyous genotpyes are unphased
        """
        if self._prephased or self.data.shape[2] < 3:
            self.log.warning("Phase information has already been removed from the data")
            return
        # check: are there any variants that are heterozygous and unphased?
        unphased = (self.data[:, :, 0] ^ self.data[:, :, 1]) & (~self.data[:, :, 2])
        if np.any(unphased):
            samp_idx, variant_idx = np.nonzero(unphased)
            raise ValueError(
                "Variant with ID {} at POS {}:{} is unphased for sample {}".format(
                    *tuple(self.variants[variant_idx[0]])[:3], self.samples[samp_idx[0]]
                )
            )
        # remove the last dimension that contains the phase info
        self.data = self.data[:, :, :2]

    def to_MAC(self):
        """
        Convert the ALT count GT matrix into a matrix of minor allele counts

        This function modifies :py:attr:`~.Genotypes.data` in-place

        It also changes the 'aaf' record in :py:attr:`~.Genotypes.variants` to 'maf'

        Raises
        ------
        AssertionError
            If the matrix has already been converted
        """
        if self.variants.dtype.names[3] == "maf":
            self.log.warning(
                "The matrix already counts instances of the minor allele rather than "
                "the ALT allele."
            )
            return
        need_conversion = self.variants["aaf"] > 0.5
        # flip the count on the variants that have an alternate allele frequency
        # above 0.5
        self.data[:, need_conversion, :2] = ~self.data[:, need_conversion, :2]
        # also encode an MAF instead of an AAF in self.variants
        self.variants["aaf"][need_conversion] = (
            1 - self.variants["aaf"][need_conversion]
        )
        # replace 'aaf' with 'maf' in the matrix
        self.variants.dtype.names = [
            (x, "maf")[x == "aaf"] for x in self.variants.dtype.names
        ]


class GenotypesRefAlt(Genotypes):
    """
    A class for processing genotypes from a file
    Unlike the base Genotypes class, this class also includes REF and ALT alleles in
    the variants array

    Attributes
    ----------
    data : np.array
        See documentation for :py:attr:`~.Genotypes.data`
    fname : Path
        See documentation for :py:attr:`~.Genotypes.fname`
    samples : tuple[str]
        See documentation for :py:attr:`~.Genotypes.samples`
    variants : np.array
        Variant-level meta information:
            1. ID
            2. CHROM
            3. POS
            4. AAF: allele freq of alternate allele (or MAF if to_MAC() is called)
            5. REF
            6. ALT
    log: Logger
        See documentation for :py:attr:`~.Genotypes.log`
    """

    def __init__(self, fname: Path, log: Logger = None):
        super(Genotypes, self).__init__(fname, log)
        self.samples = tuple()
        self.variants = np.array(
            [],
            dtype=[
                ("id", "U50"),
                ("chrom", "U10"),
                ("pos", np.uint32),
                ("aaf", np.float64),
                ("ref", "U100"),
                ("alt", "U100"),
            ],
        )
        self._prephased = False
        self._samp_idx = None
        self._var_idx = None

    def _variant_arr(self, record: Variant):
        """
        See documentation for :py:meth:`~.Genotypes._variant_arr`
        """
        return np.array(
            (
                record.ID,
                record.CHROM,
                record.POS,
                record.aaf,
                record.REF,
                record.ALT[0],
            ),
            dtype=self.variants.dtype,
        )

    def write(self):
        """
        Write the variants in this class to a VCF at :py:attr:`~.GenotypesRefAlt.fname`
        """
        vcf = VariantFile(str(self.fname), mode="w")
        # make sure the header is properly structured
        for contig in set(self.variants["chrom"]):
            vcf.header.contigs.add(contig)
        vcf.header.add_meta(
            "FORMAT",
            items=[
                ("ID", "GT"),
                ("Number", 1),
                ("Type", "String"),
                ("Description", "Genotype"),
            ],
        )
        try:
            vcf.header.add_samples(self.samples)
        except AttributeError:
            self.log.warning(
                "Upgrade to pysam >=0.19.1 to reduce the time required to create "
                "VCFs. See https://github.com/pysam-developers/pysam/issues/1104"
            )
            for sample in self.samples:
                vcf.header.add_sample(sample)
        self.log.info("Writing VCF records")
        for var_idx, var in enumerate(self.variants):
            rec = {
                "contig": var["chrom"],
                "start": var["pos"],
                "stop": var["pos"] + 1,
                "qual": None,
                "alleles": tuple(var[["ref", "alt"]]),
                "id": var["id"],
                "filter": None,
            }
            # handle pysam increasing the start site by 1
            rec["start"] -= 1
            # parse the record into a pysam.VariantRecord
            record = vcf.new_record(**rec)
            for samp_idx, sample in enumerate(self.samples):
                record.samples[sample]["GT"] = tuple(self.data[samp_idx, var_idx])
                record.samples[sample].phased = True
            # write the record to a file
            vcf.write(record)
        vcf.close()


class GenotypesPLINK(GenotypesRefAlt):
    """
    A class for processing genotypes from a PLINK .pgen file

    NOTE: this class is still under-development and NOT fit for usage

    Attributes
    ----------
    data : np.array
        See documentation for :py:attr:`~.GenotypesRefAlt.data`
    samples : tuple
        See documentation for :py:attr:`~.GenotypesRefAlt.data`
    variants : np.array
        See documentation for :py:attr:`~.GenotypesRefAlt.data`
    log: Logger
        See documentation for :py:attr:`~.GenotypesRefAlt.data`
    _prephased: bool
        See documentation for :py:attr:`~.GenotypesRefAlt.data`

    Examples
    --------
    >>> genotypes = GenotypesPLINK.load('tests/data/simple.pgen')
    """

    def read_samples(self, samples: list[str] = None):
        """
        Read sample IDs from a PSAM file into a list stored in
        :py:attr:`~.GenotypesPLINK.samples`

        One of either variants or max_variants MUST be specified!

        Parameters
        ----------
        samples : list[str], optional
            See documentation for :py:attr:`~.GenotypesRefAlt.read`

        Returns
        -------
        npt.NDArray[np.uint32]
            The indices of each of the samples within the PSAM file
        """
        if len(self.samples) != 0:
            self.log.warning("Sample data has already been loaded. Overriding.")
        if samples is not None and not isinstance(samples, set):
            samples = set(samples)
        with hook_compressed(self.fname.with_suffix(".psam"), mode="rt") as psam:
            psamples = reader(psam, delimiter="\t")
            # find the line that declares the header
            for header in psamples:
                if header[0].startswith("#FID") or header[0].startswith("#IID"):
                    break
            # we need the header to contain the IID column
            if header[0][0] != "#":
                raise ValueError("Your PSAM file is missing a header!")
                col_idx = 1
                self.log.info(
                    "Your PVAR file lacks a proper header! Assuming first five columns"
                    " are [CHROM, POS, ID, REF, ALT] in that order..."
                )
                # TODO: add header back in to psamples using itertools.chain?
            else:
                header[0] = header[0][1:]
                try:
                    col_idx = header.index("IID")
                except ValueError:
                    raise ValueError("Your PSAM file must have an IID column.")
            self.samples = {
                ct: samp[col_idx]
                for ct, samp in enumerate(psamples)
                if (samples is None) or (samp[col_idx] in samples)
            }
            indices = np.array(list(self.samples.keys()), dtype=np.uint32)
            self.samples = list(self.samples.values())
            return indices

    def _check_region(
        self, pos: tuple, chrom: str, start: int = 0, end: int = float("inf")
    ):
        """
        Check that pos lies within the provided chrom, start, and end coordinates

        This is a helper function for :py:meth:`~.GenotypesPLINK.read_variants`

        Parameters
        ----------
        pos : tuple[str, int]
            A tuple of two elements: contig (str) and chromosomal position (int)
        chrom: str
            The contig of the region to check
        start: int, optional
            The start position of the region
        end: int, optional
            The end position of the region
        """
        return (pos[0] == chrom) and (start <= pos[1]) and (end >= pos[1])

    def read_variants(
        self,
        region: str = None,
        variants: set[str] = None,
        max_variants: int = None,
    ):
        """
        Read variants from a PVAR file into a numpy array stored in
        :py:attr:`~.GenotypesPLINK.variants`

        One of either variants or max_variants MUST be specified!

        Parameters
        ----------
        region : str, optional
            See documentation for :py:attr:`~.GenotypesRefAlt.read`
        variants : set[str], optional
            See documentation for :py:attr:`~.GenotypesRefAlt.read`
        max_variants : int, optional
            See documentation for :py:attr:`~.GenotypesRefAlt.read`

        Returns
        -------
        npt.NDArray[np.uint32]
            The indices of each of the variants within the PVAR file
        """
        if len(self.variants) != 0:
            self.log.warning("Variant data has already been loaded. Overriding.")
        if variants is not None:
            max_variants = len(variants)
        if max_variants is None:
            raise ValueError("Provide either the variants or max_variants parameter!")
        # split the region string so each portion is an element
        if region is not None:
            region = re.split(":|-", region)
            if len(region) > 1:
                region[1:] = [int(pos) for pos in region[1:] if pos]
        with hook_compressed(self.fname.with_suffix(".pvar"), mode="rt") as pvar:
            pvariants = reader(pvar, delimiter="\t")
            # first, preallocate the array
            self.variants = np.empty((max_variants,), dtype=self.variants.dtype)
            # find the line that declares the header
            for header in pvariants:
                if not header[0].startswith("##"):
                    break
            # there should be at least five columns
            if len(header) < 5:
                raise ValueError("Your PVAR file should have at least five columns.")
            if header[0][0] == "#":
                header[0] = header[0][1:]
            else:
                raise ValueError("Your PVAR file is missing a header!")
                header = ["CHROM", "POS", "ID", "REF", "ALT"]
                self.log.info(
                    "Your PVAR file lacks a proper header! Assuming first five columns"
                    " are [CHROM, POS, ID, REF, ALT] in that order..."
                )
                # TODO: add header back in to pvariants using itertools.chain?
            cid = {
                item: header.index(item.upper())
                for item in self.variants.dtype.names
                if item != "aaf"
            }
            indices = np.empty((max_variants,), dtype=np.uint32)
            num_seen = 0
            for ct, rec in enumerate(pvariants):
                if region and not self._check_region(
                    (rec[cid["chrom"]], int(rec[cid["pos"]])), *region
                ):
                    continue
                if variants is not None:
                    if rec[cid["id"]] not in variants:
                        continue
                indices[num_seen] = ct
                self.variants[num_seen] = np.array(
                    (
                        rec[cid["id"]],
                        rec[cid["chrom"]],
                        rec[cid["pos"]],
                        0,
                        rec[cid["ref"]],
                        rec[cid["alt"]],
                    ),
                    dtype=self.variants.dtype,
                )
                num_seen += 1
        if max_variants > num_seen:
            self.log.info(
                f"Removing {max_variants-num_seen} unneeded variant records that "
                "were preallocated b/c max_variants was specified."
            )
            indices = indices[:num_seen]
            self.variants = self.variants[:num_seen]
        return indices

    def read(
        self,
        region: str = None,
        samples: list[str] = None,
        variants: set[str] = None,
        max_variants: int = None,
    ):
        """
        Read genotypes from a PGEN file into a numpy matrix stored in
        :py:attr:`~.GenotypesPLINK.data`

        Parameters
        ----------
        region : str, optional
            See documentation for :py:attr:`~.GenotypesRefAlt.read`
        samples : list[str], optional
            See documentation for :py:attr:`~.GenotypesRefAlt.read`
        variants : set[str], optional
            See documentation for :py:attr:`~.GenotypesRefAlt.read`
        max_variants : int, optional
            See documentation for :py:attr:`~.GenotypesRefAlt.read`
        """
        super(Genotypes, self).read()
        # TODO: figure out how to install this package
        from pgenlib import PgenReader

        sample_idxs = self.read_samples(samples)
        with PgenReader(
            bytes(str(self.fname), "utf8"), sample_subset=sample_idxs
        ) as pgen:

            if variants is not None:
                max_variants = len(variants)
            if max_variants is None:
                # use the pgen file to figure out how many variants there are
                max_variants = pgen.get_variant_ct()
            indices = self.read_variants(region, variants, max_variants)

            # the genotypes start out as a simple 2D array with twice the number of samples
            if not self._prephased:
                # raise an error message b/c this is untested code that doesn't really
                # seem to work properly
                # for ex, the phasing info is not boolean for some reason?
                raise ValueError("Not implemented yet!")
                # ...so each column is a different chromosomal strand
                self.data = np.empty(
                    (len(indices), len(sample_idxs) * 2), dtype=np.int32
                )
                phasing = np.empty((len(indices), len(sample_idxs)), dtype=np.uint8)
                # The haplotype-major mode of read_alleles_and_phasepresent_list has
                # not been implemented yet, so we need to read the genotypes in sample-
                # major mode and then transpose them
                pgen.read_alleles_and_phasepresent_list(indices, self.data, phasing)
                # missing alleles will have a value of -9
                # let's make them be -1 to be consistent with cyvcf2
                self.data[self.data == -9] = -1
                self.data = np.dstack((self.data[:, ::2], self.data[:, 1::2])).astype(
                    np.uint8
                )
                self.data = np.concatenate(
                    (self.data, phasing[:, :, np.newaxis]), axis=2
                )
                # transpose the GT matrix so that samples are rows and variants are columns
                self.data = self.data.transpose((1, 0, 2))
            else:
                # ...so each row is a different chromosomal strand
                self.data = np.empty(
                    (len(sample_idxs) * 2, len(indices)), dtype=np.int32
                )
                pgen.read_alleles_list(indices, self.data, hap_maj=True)
                # missing alleles will have a value of -9
                # let's make them be -1 to be consistent with cyvcf2
                self.data[self.data == -9] = -1
                self.data = np.dstack((self.data[::2, :], self.data[1::2, :])).astype(
                    np.uint8
                )
