import unittest
from simtools import genotypes as gp
from simtools import simtools as si
import numpy as np
import pandas as pd
import pymp

class TestSimtools(unittest.TestCase):

    """Docstring for TestSimtools. """

    def setUp(self):
        self.n = 100
        self.p = 1000
        self.n_cases = 50
        self.n_controls = 50
        self.genotypematrix = gp.simple_genotype_matrix(self.n, self.p, 0.05, 0.5)
        self.sim = si.Simtools(self.genotypematrix)

    def test_simple_sim(self):
        pheno = self.sim.simple_phenotype(0.2, 0.4)
        self.assertAlmostEqual(np.var(pheno), 1.000, delta=0.2)
        self.assertAlmostEqual(np.mean(pheno), 0.000, delta=0.2)

    def test_liabilitymodel(self):
        pheno = self.sim.simple_phenotype(0.1, 0.3, (0.4, self.n_cases, self.n_controls))
        self.assertEqual(len(pheno)-np.sum(pheno), self.n_cases)
        self.assertEqual(sum(pheno), self.n_cases)

#    def test_gwas(self):
#        #pheno = self.sim.simple_phenotype(0.1, 0.3, (0.4, self.n_cases, self.n_controls))
#        # single thread
#        pheno = np.random.rand(100)
#        self.genotypematrix = np.random.rand(100, 1000)
#        output = si.gwas(pheno, self.genotypematrix)
#        # multi thread
#        #output = si.gwas(pheno, self.genotypematrix, num_threads=4)
#        self.assertEqual(output.shape[0], self.p)
#        self.assertEqual(output.shape[1], 3)


if __name__ == '__main__':
    unittest.main()
    exit()
