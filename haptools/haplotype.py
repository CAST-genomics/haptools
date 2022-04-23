from __future__ import annotations
from dataclasses import dataclass, field

import numpy as np

from .data import Haplotype


@dataclass
class HaptoolsHaplotype(Haplotype):
    """
    A haplotype

    Properties and functions are shared with the base Haplotype object
    """
    ancestry: str
    beta: float
    _fmt: str = field(default="H\t{chrom:s}\t{start:d}\t{end:d}\t{id:s}\t{ancestry:s}\t{beta:.3f}", init=False, repr=False)

    @staticmethod
    def extras() -> tuple:
        return (
            "#H\tancestry\tstr\tLocal ancestry",
            "#H\tbeta\tfloat\tEffect size in linear model",
        )

    def transform(self, genotypes: Genotypes) -> npt.NDArray[np.bool_]:
        """
        Transform a genotypes matrix via the current haplotype:

        Each entry in the returned matrix denotes the presence of the current haplotype
        in each chromosome of each sample in the Genotypes object

        Parameters
        ----------
        genotypes : Genotypes
            The genotypes which to transform using the current haplotype

        Returns
        -------
        npt.NDArray[np.bool_]
            A 3D haplotype matrix similar to the genotype matrix but with haplotypes
            instead of variants in the columns. It will have the same shape except that
            the number of columns (second dimension) will have decreased by the number
            of variants in this haplotype.
        """
        # TODO: implement this
        pass
