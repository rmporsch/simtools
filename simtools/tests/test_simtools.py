import unittest
from simtools.simtools2 import Simtools
import numpy as np


class TestSimtools(unittest.TestCase):

    """Docstring for TestSimtools. """

    def setUp(self):
        self.plink = "data/1kg_phase1_chr22"
        self.vcf = "data/example.vcf.gz"
        self.sim_plink = Simtools(self.plink)
        self.sim_vcf = Simtools(self.vcf)

    def test_correct_read(self):
        self.assertTrue(self.sim_vcf.type == 'vcf')
        self.assertTrue(self.sim_plink.type == 'plink')
        for cl in [self.sim_vcf, self.sim_plink]:
            self.assertGreater(len(cl.variants), 1)
            self.assertGreater(len(cl.subject), 1)
            self.assertIsInstance(cl.variants[0], str)
            self.assertIsInstance(cl.subject[0], str)
            self.assertIsInstance(cl.N, int)
            self.assertIsInstance(cl.P, int)

    def test_genotypematrix(self):
        n = 100
        p = 1000
        for cl in [self.sim_vcf, self.sim_plink]:
            temp = cl.get_gentoypematrix(cl.variants[0:p], cl.subject[0:n])
            self.assertEqual(temp.shape[0], n)
            self.assertEqual(temp.shape[1], p)

    def test_simple_sim(self):
        n = 100
        for cl in [self.sim_vcf, self.sim_plink]:
            pheno = cl.simple_phenotype(0.05, 0.4, n=n)
            self.assertAlmostEqual(np.var(pheno), 1.000, delta=0.2)
            self.assertAlmostEqual(np.mean(pheno), 0.000, delta=0.2)

    def test_liability_model(self):
        n_a = n_u = 100
        liab_para = (0.4, n_a, n_u)
        for cl in [self.sim_vcf, self.sim_plink]:
            pheno = cl.simple_phenotype(0.1, 0.3, liab_para)
            self.assertEqual(len(pheno)-np.sum(pheno), n_a)
            self.assertEqual(sum(pheno), n_a)

    def test_multi_pheno(self):
        B = np.zeros((3, 3))
        lamb = np.zeros((3, 3))
        num_causal = [3, 3, 3]
        n = 100
        for cl in [self.sim_vcf, self.sim_plink]:
            pheno = cl.multi_phenotype(lamb, B, num_causal, n=n)
            self.assertEqual(pheno.shape[1], n)
            self.assertEqual(pheno.shape[0], 3)
            var = np.var(pheno, axis=1)
            self.assertAlmostEqual(var[0], 1, delta=0.2)
            self.assertAlmostEqual(var[1], 1, delta=0.2)
            self.assertAlmostEqual(var[2], 1, delta=0.2)


if __name__ == '__main__':
    unittest.main()
    exit()
