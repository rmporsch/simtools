import unittest
from simtools import genotypes as gp
import numpy as np
import pandas as pd
import os

DATA = os.path.join(os.path.dirname(__file__), 'data')

class TestGenotypes(unittest.TestCase):

    def setUp(self):
        self.n = 1000
        self.p = 100
        self.testread = gp.ReadVCF('data/example.vcf.gz')

    def test_simple_genotype(self):
        temp = gp.simple_genotype_matrix(self.n, self.p)
        self.assertEqual(temp.shape[0], self.n, 'incorrect row numbers')
        self.assertEqual(temp.shape[1], self.p, 'incorrect col numbers')
        sum_genotypes = np.sum(temp)
        self.assertLess(sum_genotypes, self.n*self.p, 'sum of all variants is >(n*p)')

    def test_vcf_read(self):
        self.testread.read_vcf(index_file=True)
        temp = pd.read_table('.temp/variant.list', header=None, sep=' ')
        self.assertEqual(temp.shape[0], 970, 'variant file seems incorrect')
        self.assertEqual(temp.shape[1], 4, 'variant file seems incorrect')

    def test_sampling(self):
        self.testread.read_vcf(index_file=False)
        temp = self.testread.sample(0.1, self.n, self.p)
        self.assertEqual(temp.shape[0], self.n, 'incorrect row numbers')
        maf = temp.mean(axis=0)
        self.assertEqual(len(maf), self.p, 'incorrect col numbers')
        self.assertFalse(np.any(maf < 0.1), 'incorrect allele frequencies')

if __name__ == '__main__':
    unittest.main()
