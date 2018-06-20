import unittest
from simtools.read_vcf import ReadVCF
import numpy as np
import vcf
import pandas as pd


class TestVCF(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.n = 200
        self.p = 100
        self.f = 'data/example.vcf.gz'
        self.testread = ReadVCF('data/example.vcf.gz')

    def test_get_samples(self):
        self.assertGreater(len(self.testread.subjects), 1)

    def test_get_genotypes(self):
        samples = ['HG00096', 'HG00097', 'HG00099']
        reader = vcf.Reader(filename=self.f)
        record = next(reader)
        gt = self.testread._get_genotypes(samples, record, False)
        self.assertTrue(np.sum(gt) == 2)
        self.assertTrue(len(gt), 3)

    def test_get_variants_info(self):
        self.testread._get_allele_freq('output.txt')
        dat = pd.read_table('output.txt', delim_whitespace=True)
        self.assertTrue(dat.shape[0] > 1)
        self.assertTrue(dat.shape[1] == 4)

    def test_genotypematrix(self):
        n = 100
        p = 200
        temp = self.testread.get_gentoypematrix(
            self.testread.variants[0:p],
            self.testread.subjects[0:n])
        self.assertEqual(temp.shape[0], n, 'incorrect row numbers')
        self.assertEqual(temp.shape[1], p, 'incorrect col numbers')
