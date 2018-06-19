import unittest
from simtools.simtools import tools
from simtools import simtools as si
import numpy as np

class TestTools(unittest.TestCase):

    def setUp(self):
        self.ff = "/home/robert/Documents/projects/risk_prediction/data/subset_10k"
        self.plink_loc = "/home/robert/software/plink"

        self.sim = si.Simtools(self.ff)
        self.pheno = self.sim.simple_phenotype(0.2, 0.6, n=100)
        self.plink = tools.Plink(self.ff, self.plink_loc)

    def test_plink(self):
        results = self.plink.gwas_plink(self.pheno,
                subjects=self.sim.last_random_subjects)

        self.assertGreater(results.shape[0], 1)
        self.assertGreater(results.shape[1], 1)
        self.assertEqual(np.max(results.NMISS), 100)

    def test_prs(self):
        results = self.plink.gwas_plink(self.pheno,
                subjects=self.sim.last_random_subjects)
        # format prs ready
        prs_computation = results[['SNP', 'A1', 'BETA']]
        prs = self.plink.prs(prs_computation,
                subjects=self.sim.last_random_subjects)
        self.assertEqual(prs.shape[0], 100)

    def test_clumping(self):
        results = self.plink.gwas_plink(self.pheno,
                subjects=self.sim.last_random_subjects)

        clumping = self.plink.clumping(results)
        self.assertGreater(clumping.shape[0], 1)

    def test_qqplot(self):
        results = self.plink.gwas_plink(self.pheno,
                subjects=self.sim.last_random_subjects)

        tools.qqplot(results.P.values)
