from __future__ import annotations
from dataclasses import dataclass, field

import numpy as np

from .data import Extra, Haplotype


@dataclass
class HaptoolsHaplotype(Haplotype):
    """
    A haplotype with sufficient fields for simphenotype

    Properties and functions are shared with the base Haplotype object
    """

    ancestry: str
    beta: float
    _extras: tuple = field(
        repr=False,
        init=False,
        default=(
            Extra("ancestry", "s", "Local ancestry"),
            Extra("beta", ".2f", "Effect size in linear model"),
        ),
    )
