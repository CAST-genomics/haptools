import os
from pathlib import Path

import pytest
import numpy as np
import numpy.lib.recfunctions as rfn

from haptools.sim_phenotype import Haplotype, PhenoSimulator
from haptools.data import (
    Genotypes,
    Phenotypes,
    Haplotypes,
)


DATADIR = Path(__file__).parent.joinpath("data")


class TestSimPhenotype:
    def _get_fake_gens(self):
        gts = Genotypes(fname=None)
        gts.data = np.array(
            [
                [[1, 0], [0, 1]],
                [[0, 1], [0, 1]],
                [[0, 0], [1, 1]],
                [[1, 1], [1, 1]],
                [[0, 0], [0, 0]],
            ],
            dtype=np.uint8,
        )
        gts.variants = np.array(
            [
                ("1:10114:T:C", "1", 10114, 0.7),
                ("1:10116:A:G", "1", 10116, 0.6),
            ],
            dtype=[
                ("id", "U50"),
                ("chrom", "U10"),
                ("pos", np.uint32),
                ("aaf", np.float64),
            ],
        )
        gts.samples = ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")
        return gts

    def _get_fake_haps(self):
        return [
            Haplotype("1", 10114, 10115, "1:10114:T:C", 0.25),
            Haplotype("1", 10116, 10117, "1:10116:A:G", 0.75),
        ]

    def _get_expected_phens(self):
        pts = Phenotypes(fname=None)
        pts.data = np.array(
            [
                [-0.13363062, False],
                [-0.13363062, False],
                [0.53452248, True],
                [1.20267559, True],
                [-1.46993683, False],
            ],
            dtype=np.float64,
        )
        pts.samples = ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")
        pts.names = ("1:10114:T:C", "1:10116:A:G")
        return pts

    def test_one_hap_zero_noise(self):
        gts = self._get_fake_gens()
        hps = self._get_fake_haps()
        expected = self._get_expected_phens()

        pt_sim = PhenoSimulator(gts, seed=42)
        data = pt_sim.run([hps[0]], heritability=1)
        data = data[:, np.newaxis]
        phens = pt_sim.phens

        # check the data and the generated phenotype object
        assert phens.data.shape == (5, 1)
        np.testing.assert_allclose(phens.data, data)
        assert phens.data[0, 0] == phens.data[1, 0]
        assert phens.data[2, 0] == phens.data[4, 0]
        assert phens.data[3, 0] > phens.data[0, 0]
        assert phens.data[2, 0] < phens.data[0, 0]
        assert phens.samples == expected.samples
        assert phens.names[0] == expected.names[0]

    def test_one_hap_zero_noise_all_same(self):
        gts = self._get_fake_gens()
        hps = self._get_fake_haps()
        expected = self._get_expected_phens()

        gts_shape = list(gts.data.shape)
        gts_shape[1] = 1
        # set the genotypes
        gts.variants = gts.variants[:1]
        gts.data = np.zeros(tuple(gts_shape), dtype=gts.data.dtype) + 1
        # set the expected phenotypes
        expected.names = expected.names[:1]
        expected.data = np.zeros((gts_shape[0], 1), dtype=expected.data.dtype)

        pt_sim = PhenoSimulator(gts, seed=42)
        data = pt_sim.run([hps[0]], heritability=1)
        data = data[:, np.newaxis]
        phens = pt_sim.phens

        # check the data and the generated phenotype object
        assert phens.data.shape == (5, 1)
        np.testing.assert_allclose(phens.data, expected.data)
        assert phens.samples == expected.samples
        assert phens.names[0] == expected.names[0]

    def test_one_hap_zero_noise_neg_beta(self):
        """
        the same test as test_one_phen_zero_noise but with a negative beta this time
        """
        gts = self._get_fake_gens()
        hps = self._get_fake_haps()
        # make the beta value negative
        hps[0].beta = -hps[0].beta
        expected = self._get_expected_phens()

        pt_sim = PhenoSimulator(gts, seed=42)
        data = pt_sim.run([hps[0]], heritability=1)
        data = data[:, np.newaxis]
        phens = pt_sim.phens

        # check the data and the generated phenotype object
        assert phens.data.shape == (5, 1)
        np.testing.assert_allclose(phens.data, data)
        assert phens.data[0, 0] == phens.data[1, 0]
        assert phens.data[2, 0] == phens.data[4, 0]
        assert phens.data[3, 0] < phens.data[0, 0]
        assert phens.data[2, 0] > phens.data[0, 0]
        assert phens.samples == expected.samples
        assert phens.names[0] == expected.names[0]

    def test_two_haps_zero_noise(self):
        gts = self._get_fake_gens()
        hps = self._get_fake_haps()
        expected = self._get_expected_phens()

        pt_sim = PhenoSimulator(gts, seed=42)
        pt_sim.run([hps[0]], heritability=1)
        pt_sim.run([hps[1]], heritability=1)
        phens = pt_sim.phens

        # check the data and the generated phenotype object
        assert phens.data.shape == (5, 2)
        assert phens.data[0, 0] == phens.data[1, 0]
        assert phens.data[2, 0] == phens.data[4, 0]
        assert phens.data[3, 0] > phens.data[0, 0]
        assert phens.data[2, 0] < phens.data[0, 0]

        assert phens.data[0, 1] == phens.data[1, 1]
        assert phens.data[2, 1] == phens.data[3, 1]
        assert phens.data[3, 1] > phens.data[0, 1]
        assert phens.data[4, 1] < phens.data[0, 1]

        assert phens.samples == expected.samples
        assert phens.names == expected.names

    def test_combined_haps_zero_noise(self):
        gts = self._get_fake_gens()
        hps = self._get_fake_haps()
        expected = self._get_expected_phens()

        pt_sim = PhenoSimulator(gts, seed=42)
        pt_sim.run(hps, heritability=1)
        phens = pt_sim.phens

        # check the data and the generated phenotype object
        assert phens.data.shape == (5, 1)
        np.testing.assert_allclose(phens.data[:, 0], expected.data[:, 0])

        assert phens.samples == expected.samples
        assert phens.names == ("-".join(expected.names),)

    def test_noise(self):
        gts = self._get_fake_gens()
        hps = self._get_fake_haps()
        expected = self._get_expected_phens()

        pt_sim = PhenoSimulator(gts, seed=42)
        pt_sim.run(hps, heritability=1)
        pt_sim.run(hps, heritability=0.96)
        pt_sim.run(hps, heritability=0.5)
        pt_sim.run(hps, heritability=0.2)
        phens = pt_sim.phens

        # check the data and the generated phenotype object
        assert phens.data.shape == (5, 4)
        np.testing.assert_allclose(phens.data[:, 0], expected.data[:, 0])
        diff1 = np.abs(phens.data[:, 1] - phens.data[:, 0]).sum()
        assert diff1 > 0
        diff2 = np.abs(phens.data[:, 2] - phens.data[:, 0]).sum()
        assert diff2 > diff1
        diff3 = np.abs(phens.data[:, 3] - phens.data[:, 0]).sum()
        assert diff3 > diff2

    def test_case_control(self):
        gts = self._get_fake_gens()
        hps = self._get_fake_haps()
        expected = self._get_expected_phens()
        all_false = np.zeros(expected.data.shape, dtype=np.float64)
        some_true = expected.data[:, 1]
        all_true = (~(all_false.astype(np.bool_))).astype(np.float64)

        pt_sim = PhenoSimulator(gts, seed=42)
        pt_sim.run(hps, heritability=1, prevalence=0.4)
        pt_sim.run(hps, heritability=1, prevalence=1)
        pt_sim.run(hps, heritability=1, prevalence=0)
        pt_sim.run(hps, heritability=0.5, prevalence=0.4)
        phens = pt_sim.phens
        assert phens.data.shape == (5, 4)
        np.testing.assert_allclose(phens.data[:, 0], some_true)
        np.testing.assert_allclose(phens.data[:, 1], all_true[:, 1])
        np.testing.assert_allclose(phens.data[:, 2], all_false[:, 1])
        diff1 = (phens.data[:, 3] == phens.data[:, 0]).sum()
        assert diff1 > 0
