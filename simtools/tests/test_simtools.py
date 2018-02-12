import unittest
from simtools import simtools as si
import numpy as np


class TestSimtools(unittest.TestCase):

    """Docstring for TestSimtools. """

    def setUp(self):
        self.n = 100
        self.p = 1000
        self.ff = "/home/robert/Documents/projects/risk_prediction/data/subset_10k"
        self.n_cases = 50
        self.n_controls = 50
        self.sim = si.Simtools(self.ff)

    def test_simple_sim(self):
        pheno = self.sim.simple_phenotype(0.05, 0.4, n=self.n)
        self.assertAlmostEqual(np.var(pheno), 1.000, delta=0.2)
        self.assertAlmostEqual(np.mean(pheno), 0.000, delta=0.2)

    def test_liability_model(self):
        pheno = self.sim.simple_phenotype(0.1, 0.3, (0.4, self.n_cases, self.n_controls))
        self.assertEqual(len(pheno)-np.sum(pheno), self.n_cases)
        self.assertEqual(sum(pheno), self.n_cases)

    def test_multi_pheno(self):
        B = np.zeros((3, 3))
        lamb = np.zeros((3, 3))
        num_causal = [3, 3, 3]
        pheno = self.sim.multi_phenotype(lamb, B, num_causal, n=self.n)
        self.assertEqual(pheno.shape[1], self.n)
        self.assertEqual(pheno.shape[0], 3)
        var = np.var(pheno, axis=1)
        self.assertAlmostEqual(var[0], 1, delta=0.2)
        self.assertAlmostEqual(var[1], 1, delta=0.2)
        self.assertAlmostEqual(var[2], 1, delta=0.2)


if __name__ == '__main__':
    unittest.main()
    exit()
