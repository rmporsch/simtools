import unittest
from simtools import genotypes as gp
from simtools import simtools as si
import numpy as np
import pandas as pd

class TestSimtools(unittest.TestCase):

    """Docstring for TestSimtools. """

    def setUp(self):
        self.n = 100
        self.p = 1000
        self.n_cases = 50
        self.n_controls = 50
        self.genotypematrix = gp.simple_genotype_matrix(self.n, self.p, 0.05, 0.5)
        self.sim = si.Simtools(self.genotypematrix)

    
    def test_causal_sim(self):
        vec_causal = self.sim.define_causal(0.1)
        self.assertEqual(len(vec_causal), self.p, 'wrong number of causal variants')
        self.assertAlmostEqual(np.mean(vec_causal), 0.1, 1)
        self.assertEqual(np.max(vec_causal), 1)
        self.assertEqual(np.min(vec_causal), 0)

    def test_causal_check(self):
        causal = np.random.binomial(1, 0.1, self.p)
        causal_check = self.sim.define_causal(causal)
        self.assertTrue(np.all(causal == causal_check))

    def test_simple_sim(self):
        pheno = self.sim.simple_phenotype(0.2, 0.4)
        self.assertAlmostEqual(np.var(pheno), 1.000, 1)
        self.assertAlmostEqual(np.mean(pheno), 0.000, 1)

    def test_liabilitymodel(self):
        pheno = self.sim.simple_phenotype(0.1, 0.3, (0.4, self.n_cases, self.n_controls))
        self.assertEqual(len(pheno)-np.sum(pheno), self.n_cases)
        self.assertEqual(sum(pheno), self.n_cases)

    def test_gwas(self):
        pheno = self.sim.simple_phenotype(0.1, 0.3, (0.4, self.n_cases, self.n_controls))
        # single thread
        output = self.sim.gwas(pheno)
        # multi thread
        output = self.sim.gwas(pheno, num_threads=4)
        self.assertEqual(output.shape[0], self.p)
        self.assertEqual(output.shape[1], 3)


if __name__ == '__main__':
    unittest.main()
