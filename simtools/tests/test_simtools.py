import unittest
from simtools.simtools2 import Simtools
import numpy as np
import os


class TestSimtools(unittest.TestCase):

    """Docstring for TestSimtools. """

    @classmethod
    def setUpClass(self):
        self.plink = "/home/robert/Documents/gits/simtools/data/1kg_phase1_chr22"
        self.vcf = "/home/robert/Documents/gits/simtools/data/example.vcf.gz"
        self.sim_plink = Simtools(self.plink)
        self.sim_vcf = Simtools(self.vcf)

    def test_simple_sim(self):
        n = 100
        for cl in [self.sim_vcf, self.sim_plink]:
            print("currently running:", cl.type)
            pheno = cl.simple_phenotype(0.05, 0.4, n=n)
            var = pheno.Pheno.var()
            mean = pheno.Pheno.mean()
            np.testing.assert_almost_equal(var, 1.000, decimal=0.2)
            np.testing.assert_almost_equal(mean, 0.000, decimal=0.2)

    def test_liability_model(self):
        n_a = n_u = 100
        liab_para = (0.4, n_a, n_u)
        for cl in [self.sim_vcf, self.sim_plink]:
            pheno = cl.simple_phenotype(0.1, 0.3, liab_para)
            pheno = pheno.Pheno.values
            self.assertEqual(len(pheno)-np.sum(pheno), n_a)
            self.assertEqual(sum(pheno), n_a)

    def test_multi_pheno(self):
        B = np.zeros((3, 3))
        lamb = np.zeros((3, 3))
        num_causal = [3, 3, 3]
        n = 200
        for cl in [self.sim_vcf, self.sim_plink]:
            pheno = cl.multi_phenotype(lamb, B, num_causal, n=n)
            self.assertEqual(pheno.shape[0], n)
            self.assertEqual(pheno.shape[1], 3)
            var = pheno.var().values
            print(var)
            self.assertAlmostEqual(var[0], 1, delta=0.2)
            self.assertAlmostEqual(var[1], 1, delta=0.2)
            self.assertAlmostEqual(var[2], 1, delta=0.2)


if __name__ == '__main__':
    unittest.main()
