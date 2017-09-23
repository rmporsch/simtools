import unittest
import genotypes as gp
import simtools as si
import numpy as np
import pandas as pd

class TestSimtools(unittest.TestCase):

    """Docstring for TestSimtools. """

    def setUp(self):
        self.genotypematrix = gp.simple_genotype_matrix(1000, 1000, 0.05, 0.5)
        self.sim = si.Simtools(self.genotypematrix)

    
    def test_causal_sim(self):
        """TODO: Docstring for test_causal.
        :returns: TODO

        """
        vec_causal = self.sim.define_causal(0.1)
        self.assertEqual(len(vec_causal), 1000, 'wrong number of causal variants')
        self.assertAlmostEqual(np.mean(vec_causal), 0.1, 1)
        self.assertEqual(np.max(vec_causal), 1)
        self.assertEqual(np.min(vec_causal), 0)

    def test_causal_check(self):
        causal = np.random.binomial(1, 0.1, 1000)
        causal_check = self.sim.define_causal(causal)
        self.assertTrue(np.all(causal == causal_check))

    def test_simple_sim(self):
        pheno = self.sim.simple_phenotype(0.2, 0.4)
        self.assertAlmostEqual(np.var(pheno), 1.0, 1)
        self.assertAlmostEqual(np.mean(pheno), 0.0, 1)

    def test_liabilitymodel(self):
        pheno = self.sim.simple_phenotype(0.1, 0.3, (0.4, 100, 100))
        self.assertEqual(len(pheno)-np.sum(pheno), 100)
        self.assertEqual(sum(pheno), 100)

    def test_gwas(self):
        pheno = self.sim.simple_phenotype(0.1, 0.3, (0.4, 100, 100))
        output = self.sim.gwas(pheno)


if __name__ == '__main__':
    unittest.main()
