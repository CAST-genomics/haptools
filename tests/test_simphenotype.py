import os
from pathlib import Path

import pytest
import numpy as np
import numpy.lib.recfunctions as rfn

from haptools.sim_phenotype import Haplotype
from haptools.data import (
    Genotypes,
    GenotypesRefAlt,
    Phenotypes,
    Haplotypes,
    Variant,
    Haplotype as HaplotypeBase,
)


DATADIR = Path(__file__).parent.joinpath("data")

def test_simphenotype(self):
    pass

