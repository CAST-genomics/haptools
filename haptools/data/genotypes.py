from __future__ import annotations
import re
import gc
from csv import reader
from pathlib import Path
from typing import Iterator
from logging import getLogger, Logger
from collections import namedtuple, Counter

import numpy as np
import numpy.typing as npt
from cyvcf2 import VCF, Variant
from pysam import VariantFile, TabixFile

try:
    import trtools.utils.tr_harmonizer as trh
except ModuleNotFoundError:
    from . import tr_harmonizer as trh

from .data import Data


class Genotypes(Data):
    """
    A class for processing genotypes from a file

    Attributes
    ----------
    data : npt.NDArray
        The genotypes in an n (samples) x p (variants) x 2 (strands) array
    fname : Path | str
        The path to the read-only file containing the data
    samples : tuple[str]
        The names of each of the n samples
    variants : np.array
        Variant-level meta information:
            1. ID
            2. CHROM
            3. POS
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

    def __init__(self, fname: Path | str, log: Logger = None):
        super().__init__(fname, log)
        self.samples = tuple()
        self.variants = np.array(
            [],
            dtype=[
                ("id", "U50"),
                ("chrom", "U10"),
                ("pos", np.uint32),
            ],
        )
        self._prephased = False
        self._samp_idx = None
        self._var_idx = None

    @classmethod
    def load(
        cls: Genotypes,
        fname: Path | str,
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
        Genotypes
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
                if num_seen >= max_variants:
                    break
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
        if 0 in self.data.shape:
            self.log.warning(
                "Failed to load genotypes. If you specified a region, check that the"
                " contig name matches! For example, double-check the 'chr' prefix."
            )
        # transpose the GT matrix so that samples are rows and variants are columns
        self.log.info(f"Transposing genotype matrix of size {self.data.shape}")
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
            (record.ID, record.CHROM, record.POS),
            dtype=self.variants.dtype,
        )

    def _vcf_iter(self, vcf: cyvcf2.VCF, region: str):
        """
        Yield all variants within a region in the VCF file.

        Parameters
        ----------
        vcf: VCF
            The cyvcf2.VCF object from which to fetch variant records
        region : str, optional
            See documentation for :py:meth:`~.Genotypes.read`

        Returns
        -------
        vcffile : cyvcf2.VCF
            Iterable cyvcf2 instance.
        """
        return vcf(region)

    def _return_data(self, variant: Variant):
        """
        Collect genotypes from current variant

        Parameters
        ----------
        variant: cyvcf2.Variant
            A cyvcf2.Variant object from which to fetch genotypes

        Returns
        -------
        data: npt.NDArray[np.uint8]
            Numpy array storing all genotypes
        """
        return np.array(variant.genotypes, dtype=np.uint8)

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
        num_seen = 0
        # iterate over each line in the VCF
        # note, this can take a lot of time if there are many samples
        for variant in self._vcf_iter(vcf, region):
            if variants is not None and variant.ID not in variants:
                if num_seen >= len(variants):
                    # exit early if we've already found all the variants
                    break
                continue
            # save meta information about each variant
            variant_arr = self._variant_arr(variant)
            # extract the genotypes to a matrix of size n x 3
            # the last dimension has three items:
            # 1) presence of REF in strand one
            # 2) presence of REF in strand two
            # 3) whether the genotype is phased (if self._prephased is False)
            data = self._return_data(variant)
            data = data[:, : (2 + (not self._prephased))]
            yield Record(data, variant_arr)
            num_seen += 1
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
        The time complexity of this function should be roughly O(n+p) if both
        parameters are True. Otherwise, it will be either O(n) or O(p).

        Parameters
        ----------
        samples: bool, optional
            Whether to index the samples for fast loop-up. Adds complexity O(n).
        variants: bool, optional
            Whether to index the variants for fast look-up. Adds complexity O(p).

        Raises
        ------
        ValueError
            If any samples or variants appear more than once
        """
        if samples and self._samp_idx is None:
            self._samp_idx = dict(zip(self.samples, range(len(self.samples))))
            if len(self._samp_idx) < len(self.samples):
                duplicates = Counter(self.samples).items()
                duplicates = [samp_id for samp_id, count in duplicates if count > 1]
                a_few = 5 if len(duplicates) > 5 else len(duplicates)
                raise ValueError(f"Found duplicate sample IDs: {duplicates[:a_few]}")
        if variants and self._var_idx is None:
            self._var_idx = dict(zip(self.variants["id"], range(len(self.variants))))
            if len(self._var_idx) < len(self.variants):
                duplicates = Counter(self.variants["id"]).items()
                duplicates = [var_id for var_id, count in duplicates if count > 1]
                a_few = 5 if len(duplicates) > 5 else len(duplicates)
                raise ValueError(f"Found duplicate variant IDs: {duplicates[:a_few]}")

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
        # A genotype value equal to the max or one less than max for uint8 indicates
        #   the value was missing
        missing = np.any(self.data[:, :, :2] >= np.iinfo(np.uint8).max - 1, axis=2)
        if np.any(missing):
            samp_idx, variant_idx = np.nonzero(missing)
            if discard_also:
                original_num_samples = len(self.samples)
                self.data = np.delete(self.data, samp_idx, axis=0)
                self.samples = tuple(np.delete(self.samples, samp_idx))
                self.log.warning(
                    "Ignoring missing genotypes from "
                    f"{original_num_samples - len(self.samples)} samples"
                )
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
        data = self.data
        if data.dtype != np.bool_:
            data = self.data.astype(np.bool_)
        unphased = (data[:, :, 0] ^ data[:, :, 1]) & (~data[:, :, 2])
        if np.any(unphased):
            samp_idx, variant_idx = np.nonzero(unphased)
            raise ValueError(
                "Variant with ID {} at POS {}:{} is unphased for sample {}".format(
                    *tuple(self.variants[variant_idx[0]])[:3], self.samples[samp_idx[0]]
                )
            )
        # remove the last dimension that contains the phase info
        self.data = self.data[:, :, :2]

    def check_maf(
        self,
        threshold: float = None,
        discard_also: bool = False,
        warn_only: bool = False,
    ) -> npt.NDArray[np.float64]:
        """
        Check the minor allele frequency of each variant

        Raise a ValueError if any variant's MAF doesn't satisfy the threshold, if
        one is provided

        .. note::
            You should call :py:meth:`~.Genotypes.check_missing` and
            :py:meth:`~.Genotypes.check_biallelic` before executing this method, for
            best results. Otherwise, the frequencies may be computed incorrectly.

        Parameters
        ----------
        threshold: float, optional
            If a variant has a minor allele frequency (MAF) rarer than this threshold,
            raise a ValueError
        discard_also : bool, optional
            If True, discard any variants that would otherwise cause a ValueError

            This parameter will be ignored if a threshold is not specified
        warn_only: bool, optional
            Just raise a warning instead of a ValueError

        Raises
        ------
        ValueError
            If any variant does not meet the provided threshold minor allele frequency

        Returns
        -------
            The minor allele frequency of each variant
        """
        num_strands = 2 * self.data.shape[0]
        # TODO: make this work for multi-allelic variants, too?
        ref_af = self.data[:, :, :2].astype(np.bool_).sum(axis=(0, 2)) / num_strands
        maf = np.array([ref_af, 1 - ref_af]).min(axis=0)
        if threshold is None:
            return maf
        rare_variants = maf < threshold
        if np.any(rare_variants):
            idx = np.nonzero(rare_variants)[0]
            if discard_also:
                original_num_variants = len(self.variants)
                self.data = np.delete(self.data, idx, axis=1)
                self.variants = np.delete(self.variants, idx)
                maf = np.delete(maf, idx)
                self.log.info(
                    "Ignoring missing genotypes from "
                    f"{original_num_variants - len(self.variants)} samples"
                )
                self._var_idx = None
            else:
                vals = tuple(self.variants[idx[0]])[:3] + (maf[idx[0]], threshold)
                msg = "Variant with ID {} at POS {}:{} has MAF {} < {}".format(*vals)
                if warn_only:
                    self.log.warning(msg)
                else:
                    # raise error if the minor allele frequency of a variant does not
                    # meet the threshold
                    raise ValueError(msg)
        return maf

    def check_sorted(self):
        """
        Check that the variant coordinates are sorted

        Raise a ValueError if any of variants in any of the chromosomes are unsorted

        Raises
        ------
        ValueError
            If any variant position is less than a position that comes before it within
            the same chromosome
        """
        chroms = set(self.variants["chrom"])
        for chrom in chroms:
            positions = self.variants["pos"][self.variants["chrom"] == chrom]
            if not np.all(positions[:-1] <= positions[1:]):
                raise ValueError(
                    f"The variants in chromosome '{chrom}' are not sorted by position"
                )

    @classmethod
    def merge_variants(
        cls, objs: tuple[Genotypes], check_samples: bool = True, **kwargs
    ) -> Genotypes:
        """
        Merge genotypes objects with different sets of variants together

        .. note::
            The input genotypes objects are not expected to have any overlapping sets
            of variants. Also, all samples in the input genotypes must be the same.

        Parameters
        ----------
        objs: tuple[Genotypes]
            The objects that should be merged together
        check_samples: bool, optional
            Whether to check that the set of provided samples is *exactly* the same
            for all genotypes. This can take a while so you may want to avoid it
        **kwargs
            Any parameters to pass to :py:meth:`~.Genotypes._init__`

        Raises
        ------
        ValueError
            If the set of samples in each input object is not the same

        Returns
        -------
        Genotypes
            A new object containing merged versions of the properties in each object
        """
        gts = cls(**kwargs)
        if check_samples:
            for obj in objs[1:]:
                if objs[0].samples != obj.samples:
                    raise ValueError("Samples must be shared among all Genotypes")
        else:
            num_samps = [len(obj.samples) for obj in objs]
            if all(num_samps[0] == num_samps[1:]):
                gts.samples = tuple(samp for obj in objs for samp in obj.samples)
            else:
                raise ValueError("Samples must be shared among all Genotypes")
        gts.samples = objs[0].samples
        dtypes = list(gts.variants.dtype.names)
        gts.variants = np.concatenate(tuple(obj.variants[dtypes] for obj in objs))
        unphased = [obj.data.shape[2] == 3 for obj in objs]
        # check: do we have a mix of phased and unphased objects?
        if any(unphased) and not all(unphased):
            data = (
                obj.data if phase else np.insert(obj.data, 2, 1, axis=2)
                for phase, obj in zip(unphased, objs)
            )
        else:
            data = (obj.data for obj in objs)
        # TODO: fix Genotypes.check_biallelic so it always keeps data as np.uint8 and then adjust this code accordingly
        dtype = (
            np.bool_ if all(obj.data.dtype == np.bool_ for obj in objs) else np.uint8
        )
        gts.data = np.concatenate(tuple(data), axis=1, dtype=dtype)
        return gts


class GenotypesVCF(Genotypes):
    """
    A class for processing genotypes from a file
    Unlike the base Genotypes class, this class also includes REF and ALT alleles as
    a list of alleles in the variants array

    Attributes
    ----------
    data : npt.NDArray
        See documentation for :py:attr:`~.Genotypes.data`
    fname : Path | str
        See documentation for :py:attr:`~.Genotypes.fname`
    samples : tuple[str]
        See documentation for :py:attr:`~.Genotypes.samples`
    variants : np.array
        Variant-level meta information:
            1. ID
            2. CHROM
            3. POS
            4. [REF, ALT1, ALT2, ...]
    log: Logger
        See documentation for :py:attr:`~.Genotypes.log`
    """

    def __init__(self, fname: Path | str, log: Logger = None):
        super().__init__(fname, log)
        dtype = {k: v[0] for k, v in self.variants.dtype.fields.items()}
        self.variants = np.array([], dtype=list(dtype.items()) + [("alleles", object)])

    def _variant_arr(self, record: Variant):
        """
        See documentation for :py:meth:`~.Genotypes._variant_arr`
        """
        return np.array(
            (
                record.ID,
                record.CHROM,
                record.POS,
                (record.REF, *record.ALT),
            ),
            dtype=self.variants.dtype,
        )

    def write(self):
        """
        Write the variants in this class to a VCF at :py:attr:`~.GenotypesVCF.fname`
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
        phased = self._prephased or (self.data.shape[2] < 3)
        missing_val = np.iinfo(np.uint8).max
        for var_idx, var in enumerate(self.variants):
            rec = {
                "contig": var["chrom"],
                "start": var["pos"],
                "stop": var["pos"] + len(var["alleles"][0]) - 1,
                "qual": None,
                "alleles": var["alleles"],
                "id": var["id"],
                "filter": None,
            }
            # handle pysam increasing the start site by 1
            rec["start"] -= 1
            # parse the record into a pysam.VariantRecord
            record = vcf.new_record(**rec)
            for samp_idx, sample in enumerate(self.samples):
                record.samples[sample]["GT"] = tuple(
                    None if val == missing_val else val
                    for val in self.data[samp_idx, var_idx, :2]
                )
                # add proper phasing info
                if phased:
                    record.samples[sample].phased = True
                else:
                    record.samples[sample].phased = self.data[samp_idx, var_idx, 2]
            # write the record to a file
            vcf.write(record)
        vcf.close()


class GenotypesTR(Genotypes):
    """
    A class for processing TR genotypes from a file
    Unlike the base Genotypes class, this class genotypes will be repeat number
    in the variants array

    Attributes
    ----------
    data : npt.NDArray
        See documentation for :py:attr:`~.Genotypes.data`
    fname : Path | str
        See documentation for :py:attr:`~.Genotypes.fname`
    samples : tuple[str]
        See documentation for :py:attr:`~.Genotypes.samples`
    variants : np.array
        Variant-level meta information:
            1. ID
            2. CHROM
            3. POS
    log: Logger
        See documentation for :py:attr:`~.Genotypes.log`
    vcftype: str
        TR vcf type currently being read.
        {'auto', 'gangstr', 'advntr', 'hipstr', 'eh', 'popstr'}
    """

    def __init__(self, fname: Path | str, log: Logger = None, vcftype: str = "auto"):
        super().__init__(fname, log)
        self.vcftype = vcftype

    @classmethod
    def load(
        cls: GenotypesTR,
        fname: Path | str,
        region: str = None,
        samples: list[str] = None,
        variants: set[str] = None,
        vcftype: str = "auto",
    ) -> Genotypes:
        """
        Load STR genotypes from a VCF file

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
        Genotypes
            A Genotypes object with the data loaded into its properties
        """
        genotypes = cls(fname, vcftype=vcftype)
        genotypes.read(region, samples, variants)
        genotypes.check_phase()
        return genotypes

    def _vcf_iter(self, vcf: cyvcf2.VCF, region: str = None):
        """
        Collect GTs (trh.TRRecord objects) to iterate over

        Parameters
        ----------
        vcf: VCF
            The cyvcf2.VCF object from which to fetch variant records
        region : str, optional
            See documentation for :py:meth:`~.Genotypes.read`

        Returns
        -------
        tr_records: trh.TRRecord
            TRRecord objects yielded from TRRecordHarmonizer
        """
        for record in trh.TRRecordHarmonizer(
            vcffile=vcf, vcfiter=vcf(region), region=region, vcftype=self.vcftype
        ):
            record.ID = record.record_id
            record.CHROM = record.chrom
            record.POS = record.pos
            yield record

    def _return_data(self, variant: trh.TRRecord):
        """
        Collect Genotypes, transform them to copy number, and return them.

        Parameters
        ----------
        variant: trh.TRRecord
            A trh.TRRecord object from which to collect copy number genotypes using the
            GetLengthGenotypes() function

        Returns
        -------
        data: npt.NDArray[np.uint8]
            Numpy array storing all genotypes
        """
        # Grab GT Lengths and round to lowest integer
        gts = np.rint(variant.GetLengthGenotypes())

        # If only one GT present fill rest with empty gts (-1)
        if gts.shape[1] == 2:
            self.log.warning(
                "The current variant in the VCF only has one allele per sample."
            )
            data = []
            # Only one GT so phase will always be 0
            zeros = np.zeros((gts.shape[0], 1))
            missing = -1 * np.ones((gts.shape[0],))
            gts = np.concatenate((gts, zeros), axis=1)
            gts[:, 1] = missing

        return gts.astype(np.uint8)

    def check_biallelic(self):
        """
        See documentation for :py:meth:`~.Genotypes.check_biallelic`
        """
        raise NotImplementedError

    def check_maf(self):
        """
        See documentation for :py:meth:`~.Genotypes.check_maf`
        """
        raise NotImplementedError


class GenotypesPLINK(GenotypesVCF):
    """
    A class for processing genotypes from a PLINK ``.pgen`` file

    Attributes
    ----------
    data : npt.NDArray
        See documentation for :py:attr:`~.GenotypesVCF.data`
    samples : tuple
        See documentation for :py:attr:`~.GenotypesVCF.data`
    variants : np.array
        See documentation for :py:attr:`~.GenotypesVCF.data`
    log: Logger
        See documentation for :py:attr:`~.GenotypesVCF.data`
    chunk_size: int, optional
        The max number of variants to fetch from and write to the PGEN file at any
        given time

        If this value is provided, variants from the PGEN file will be loaded in
        chunks so as to use less memory
    _prephased: bool
        See documentation for :py:attr:`~.GenotypesVCF.data`

    Examples
    --------
    >>> genotypes = GenotypesPLINK.load('tests/data/simple.pgen')
    """

    def __init__(self, fname: Path | str, log: Logger = None, chunk_size: int = None):
        super().__init__(fname, log)
        self.chunk_size = chunk_size
        try:
            global pgenlib
            import pgenlib
        except ImportError:
            raise ImportError(
                f"We cannot read PGEN files without the pgenlib library. Please "
                f"reinstall haptools with the 'files' extra requirement via\n"
                f"pip install haptools[files]"
            )

    def read_samples(self, samples: list[str] = None):
        """
        Read sample IDs from a PSAM file into a list stored in
        :py:attr:`~.GenotypesPLINK.samples`

        This method is called automatically by :py:meth:`~.GenotypesPLINK.read`

        Parameters
        ----------
        samples : list[str], optional
            See documentation for :py:attr:`~.GenotypesVCF.read`

        Returns
        -------
        npt.NDArray[np.uint32]
            The indices of each of the samples within the PSAM file
        """
        if len(self.samples) != 0:
            self.log.warning("Sample data has already been loaded. Overriding.")
        if samples is not None and not isinstance(samples, set):
            samples = set(samples)
        with self.hook_compressed(self.fname.with_suffix(".psam"), mode="r") as psam:
            psamples = reader(psam, delimiter="\t")
            # find the line that declares the header
            for header in psamples:
                if header[0].startswith("#FID") or header[0].startswith("#IID"):
                    break
            # we need the header to contain the IID column
            if header[0][0] != "#":
                raise ValueError("Your PSAM file is missing a header!")
                # TODO: add proper support for missing headers
                # The PSAM spec says "If no header lines are present, the file is
                # treated as if it had a '#FID IID PAT MAT SEX PHENO1' header line if
                # it has 6 or more columns, and '#FID IID PAT MAT SEX' if there are
                # exactly five"
                header = ["FID", "IID", "PAT", "MAT", "SEX"]
                self.log.info(
                    "Your PSAM file lacks a proper header! Assuming first five columns"
                    " are [FID, IID, PAT, MAT, SEX] in that order..."
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
            self.samples = tuple(self.samples.values())
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

    def _variant_arr(
        self,
        record: list[str],
        cid: dict[str, int] = dict(zip(["CHROM", "POS", "ID", "REF", "ALT"], range(5))),
    ):
        """
        Construct a np array from the metadata in a line of the PVAR file

        This is a helper function for :py:meth:`~.GenotypesPLINK._iterate_variants`.
        It's separate so that it can easily be overridden in any child classes.

        Parameters
        ----------
        record: list[str]
            A list of the tab-separated fields in a line from the PVAR file
        cid: dict[str, int]
            A dictionary mapping each column of the PVAR file to its index, as declared
            in the header

            This helps cover cases where the order of the fields might be different

        Returns
        -------
        npt.NDArray
            A row from the :py:attr:`~.GenotypesPLINK.variants` array
        """
        return np.array(
            (
                record[cid["ID"]],
                record[cid["CHROM"]],
                record[cid["POS"]],
                (record[cid["REF"]], record[cid["ALT"]]),
            ),
            dtype=self.variants.dtype,
        )

    def _iterate_variants(
        self,
        region: str = None,
        variants: set[str] = None,
    ):
        """
        A generator over the lines of a PVAR file

        This is a helper function for :py:meth:`~.GenotypesPLINK._iterate` and
        :py:meth:`~.GenotypesPLINK.read_variants`

        Parameters
        ----------
        region : str, optional
            See documentation for :py:attr:`~.GenotypesVCF.read`
        variants : set[str], optional
            See documentation for :py:attr:`~.GenotypesVCF.read`

        Yields
        ------
        Iterator[tuple[int, npt.NDArray]]
            An iterator of tuples over each line in the file

            The first value is the index of the variant and the second is a line from
            the file, encoded as a numpy mixed array type
        """
        # split the region string so each portion is an element
        if region is not None:
            region = re.split(":|-", region)
            if len(region) > 1:
                region[1:] = [int(pos) for pos in region[1:] if pos]
        with self.hook_compressed(self.fname.with_suffix(".pvar"), mode="r") as pvar:
            pvariants = reader(pvar, delimiter="\t")
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
                # TODO: add proper support for missing headers
                # The PVAR spec says "If no header lines are present, the file is
                # treated as if it had a '#CHROM ID CM POS ALT REF' header line if it
                # has 6 or more columns, and '#CHROM ID POS ALT REF' if there are
                # exactly five"
                header = ["CHROM", "POS", "ID", "REF", "ALT"]
                self.log.info(
                    "Your PVAR file lacks a proper header! Assuming first five columns"
                    " are [CHROM, POS, ID, REF, ALT] in that order..."
                )
                # TODO: add header back in to pvariants using itertools.chain?
            # create a dictionary that translates between the variant dtypes and thee
            # columns of the PVAR file
            cid = {item: header.index(item.upper()) for item in ("chrom", "pos", "id")}
            num_seen = 0
            for ct, rec in enumerate(pvariants):
                if region and not self._check_region(
                    (rec[cid["chrom"]], int(rec[cid["pos"]])), *region
                ):
                    continue
                if variants is not None:
                    if rec[cid["id"]] not in variants:
                        if num_seen >= len(variants):
                            # exit early if we've already found all the variants
                            break
                        continue
                yield ct, self._variant_arr(rec)

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

        This method is called automatically by :py:meth:`~.GenotypesPLINK.read`

        Parameters
        ----------
        region : str, optional
            See documentation for :py:attr:`~.GenotypesVCF.read`
        variants : set[str], optional
            See documentation for :py:attr:`~.GenotypesVCF.read`
        max_variants : int, optional
            See documentation for :py:attr:`~.GenotypesVCF.read`

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
        # first, preallocate the array and the indices of each variant
        self.variants = np.empty((max_variants,), dtype=self.variants.dtype)
        indices = np.empty((max_variants,), dtype=np.uint32)
        num_seen = 0
        # iterate through each variant and save it
        for ct, rec in self._iterate_variants(region, variants):
            if num_seen >= max_variants:
                break
            indices[num_seen] = ct
            self.variants[num_seen] = rec
            num_seen += 1
        # did we see more variants than was asked for?
        if max_variants > num_seen:
            self.log.info(
                f"Removing {max_variants-num_seen} unneeded variant records that "
                "were preallocated b/c max_variants was specified."
            )
            indices = indices[:num_seen]
            self.variants = self.variants[:num_seen]
        if not len(indices):
            self.log.warning(
                "Failed to load any variants. If you specified a region, check that "
                "the contig name matches! For example, double-check the 'chr' prefix."
            )
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
            See documentation for :py:attr:`~.GenotypesVCF.read`
        samples : list[str], optional
            See documentation for :py:attr:`~.GenotypesVCF.read`
        variants : set[str], optional
            See documentation for :py:attr:`~.GenotypesVCF.read`
        max_variants : int, optional
            See documentation for :py:attr:`~.GenotypesVCF.read`
        """
        super(Genotypes, self).read()
        import pgenlib

        sample_idxs = self.read_samples(samples)
        with pgenlib.PgenReader(
            bytes(str(self.fname), "utf8"), sample_subset=sample_idxs
        ) as pgen:
            # how many variants to load?
            if variants is not None:
                max_variants = len(variants)
            # use the pgen file to figure out how many variants there are
            if max_variants is None:
                max_variants = pgen.get_variant_ct()
            else:
                max_variants = min(max_variants, pgen.get_variant_ct())
            indices = self.read_variants(region, variants, max_variants)
            mat_shape = (len(sample_idxs), len(indices), (2 + (not self._prephased)))
            self.log.debug(
                f"Allocating memory for genotype matrix of shape {mat_shape} and "
                "dtype np.uint8"
            )
            # initialize the data array
            self.data = np.empty(mat_shape, dtype=np.uint8)
            # how many variants should we load at once?
            chunks = self.chunk_size
            if chunks is None or chunks > len(indices):
                chunks = len(indices)
            self.log.info(
                f"Reading genotypes from {len(self.samples)} samples and "
                f"{len(indices)} variants in chunks of size {chunks} variants"
            )
            # iterate through chunks of variants
            for start in range(0, len(indices), chunks):
                end = start + chunks
                if end > len(indices):
                    end = len(indices)
                size = end - start
                self.log.debug(f"Loading from variant #{start} to variant #{end}")
                # the genotypes start out as a simple 2D array with twice the number
                # of samples
                if not self._prephased:
                    # ...each column is a different chromosomal strand
                    try:
                        data = np.empty((size, len(sample_idxs) * 2), dtype=np.int32)
                        phasing = np.zeros((size, len(sample_idxs)), dtype=np.uint8)
                    except np.core._exceptions._ArrayMemoryError as e:
                        raise ValueError(
                            "You don't have enough memory to load these genotypes! Try"
                            " specifying a value to the chunk_size parameter, instead"
                        ) from e
                    # The haplotype-major mode of read_alleles_and_phasepresent_list
                    # has not been implemented yet, so we need to read the genotypes
                    # in sample-major mode and then transpose them
                    pgen.read_alleles_and_phasepresent_list(
                        indices[start:end], data, phasing
                    )
                    # missing alleles will have a value of -9
                    # let's make them be -1 to be consistent with cyvcf2
                    data[data == -9] = -1
                    # add phase info, then transpose the GT matrix so that samples are
                    # rows and variants are columns
                    self.data[:, start:end, :2] = data.reshape(
                        (size, mat_shape[0], 2)
                    ).transpose((1, 0, 2))
                    self.data[:, start:end, 2] = phasing.transpose()
                else:
                    # ...each row is a different chromosomal strand
                    data = np.empty((size, len(sample_idxs) * 2), dtype=np.int32)
                    pgen.read_alleles_list(indices[start:end], data)
                    # missing alleles will have a value of -9
                    # let's make them be -1 to be consistent with cyvcf2
                    data[data == -9] = -1
                    self.data[:, start:end] = data.reshape(
                        (size, mat_shape[0], 2)
                    ).transpose((1, 0, 2))
                del data
                gc.collect()

    def _iterate(
        self,
        pgen: pgenlib.PgenReader,
        region: str = None,
        variants: set[str] = None,
    ):
        """
        A generator over the lines of a PGEN-PVAR file pair

        This is a helper function for :py:meth:`~.GenotypesPLINK.__iter__`

        Parameters
        ----------
        pgen: pgenlib.PgenReader
            The pgenlib.PgenReader object from which to fetch variant records
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

        # iterate over each line in the PVAR file
        for idx, variant_arr in self._iterate_variants(region, variants):
            # the genotypes start out as a simple 2D array with twice the number of samples
            data = np.empty(len(self.samples) * 2, dtype=np.int32)
            if not self._prephased:
                phasing = np.empty(len(self.samples), dtype=np.uint8)
                # The haplotype-major mode of read_alleles_and_phasepresent_list has
                # not been implemented yet, so we need to read the genotypes in sample-
                # major mode and then transpose them
                pgen.read_alleles_and_phasepresent(idx, data, phasing)
            else:
                pgen.read_alleles(idx, data)
            # missing alleles will have a value of -9
            # let's make them be -1 to be consistent with cyvcf2
            data[data == -9] = -1
            # strand 1 is at even indices and strand 2 is at odd indices
            data = np.dstack((data[::2], data[1::2]))[0].astype(np.uint8)
            # concatenate phasing info into the data
            if not self._prephased:
                data = np.concatenate((data, phasing[:, np.newaxis]), axis=1)
            # we extracted the genotypes to a matrix of size p x 3
            # the last dimension has three items:
            # 1) presence of REF in strand one
            # 2) presence of REF in strand two
            # 3) whether the genotype is phased (if self._prephased is False)
            yield Record(data, variant_arr)
        pgen.close()

    def __iter__(
        self, region: str = None, samples: list[str] = None, variants: set[str] = None
    ) -> Iterator[namedtuple]:
        """
        Read genotypes from a PGEN line by line without storing anything

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
            See documentation for :py:meth:`~.GenotypesPLINK._iterate`
        """
        super(Genotypes, self).read()
        import pgenlib

        sample_idxs = self.read_samples(samples)
        pgen = pgenlib.PgenReader(
            bytes(str(self.fname), "utf8"), sample_subset=sample_idxs
        )
        # call another function to force the lines above to be run immediately
        # see https://stackoverflow.com/a/36726497
        return self._iterate(pgen, region, variants)

    def write_samples(self):
        """
        Write sample IDs to a PSAM file from a list stored in
        :py:attr:`~.GenotypesPLINK.samples`

        This method is called automatically by :py:meth:`~.GenotypesPLINK.write`
        """
        with self.hook_compressed(self.fname.with_suffix(".psam"), mode="w") as psam:
            psam.write("#IID\n")
            psam.write("\n".join(self.samples))
            psam.write("\n")

    def write_variants(self):
        """
        Write variant IDs to a PVAR file from the numpy array stored in
        :py:attr:`~.GenotypesPLINK.variants`

        This method is called automatically by :py:meth:`~.GenotypesPLINK.write`
        """
        with VariantFile(str(self.fname.with_suffix(".pvar")), mode="w") as vcf:
            # make sure the header is properly structured
            for contig in set(self.variants["chrom"]):
                vcf.header.contigs.add(contig)
            self.log.info("Writing PVAR records")
            for var_idx, var in enumerate(self.variants):
                rec = {
                    "contig": var["chrom"],
                    "start": var["pos"],
                    "stop": var["pos"] + len(var["alleles"][0]) - 1,
                    "qual": None,
                    "alleles": var["alleles"],
                    "id": var["id"],
                    "filter": None,
                }
                # handle pysam increasing the start site by 1
                rec["start"] -= 1
                # parse the record into a pysam.VariantRecord
                record = vcf.new_record(**rec)
                # write the record to a file
                vcf.write(record)

    def write(self):
        """
        Write the variants in this class to PLINK2 files at
        :py:attr:`~.GenotypesPLINK.fname`
        """
        import pgenlib

        # write the psam and pvar files
        self.write_samples()
        self.write_variants()
        self.log.debug(f"Transposing genotype matrix of size {self.data.shape}")
        # transpose the data b/c pgenwriter expects things in "variant-major" order
        # (ie where variants are rows instead of samples)
        data = self.data.transpose((1, 0, 2))[:, :, :2]
        # how many variants should we write at once?
        chunks = self.chunk_size
        if chunks is None or chunks > len(self.variants):
            chunks = len(self.variants)

        # write the pgen file
        with pgenlib.PgenWriter(
            filename=bytes(str(self.fname), "utf8"),
            sample_ct=len(self.samples),
            variant_ct=len(self.variants),
            nonref_flags=False,
            hardcall_phase_present=True,
        ) as pgen:
            self.log.info(
                f"Writing genotypes from {len(self.samples)} samples and "
                f"{len(self.variants)} variants in chunks of size {chunks} variants"
            )
            # iterate through chunks of variants
            for start in range(0, len(self.variants), chunks):
                end = start + chunks
                if end > len(self.variants):
                    end = len(self.variants)
                size = end - start
                try:
                    missing = np.ascontiguousarray(
                        data[start:end] == np.iinfo(np.uint8).max
                    )
                    subset_data = np.ascontiguousarray(data[start:end], dtype=np.int32)
                except np.core._exceptions._ArrayMemoryError as e:
                    raise ValueError(
                        "You don't have enough memory to write these genotypes! Try"
                        " specifying a value to the chunk_size parameter, instead"
                    ) from e
                subset_data.resize((len(self.variants), len(self.samples) * 2))
                missing.resize((len(self.variants), len(self.samples) * 2))
                # convert any missing genotypes to -9
                subset_data[missing] = -9
                if self._prephased or self.data.shape[2] < 3:
                    pgen.append_alleles_batch(subset_data, all_phased=True)
                else:
                    # TODO: figure out why this sometimes leads to a corrupted file?
                    subset_phase = self.data[:, start:end, 2].T.copy(order="C")
                    pgen.append_partially_phased_batch(subset_data, subset_phase)
                    del subset_phase
            del subset_data
            del missing
            gc.collect()
