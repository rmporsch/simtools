import unittest
from simtools import genotypes as gp
import numpy as np
import vcf
import pandas as pd
import os

DATA = os.path.join(os.path.dirname(__file__), 'data')

class TestVCF(unittest.TestCase):

    def setUp(self):
        self.n = 1000
        self.p = 100
        self.f = 'data/example.vcf.gz'
        self.testread = gp.ReadVCF('data/example.vcf.gz')

    def test_simple_genotype(self):
        temp = gp.simple_genotype_matrix(self.n, self.p)
        self.assertEqual(temp.shape[0], self.n, 'incorrect row numbers')
        self.assertEqual(temp.shape[1], self.p, 'incorrect col numbers')
        sum_genotypes = np.sum(temp)
        self.assertLess(sum_genotypes, self.n*self.p, 'sum of all variants is >(n*p)')

    def test_get_samples(self):
        self.assertGreater(len(self.testread._samples), 1)

    def test_get_genotypes(self):
        samples = ['HG00096', 'HG00097', 'HG00099']
        reader = vcf.Reader(filename=self.f)
        record = next(reader)
        gt = self.testread._get_genotypes(samples, record, False)
        self.assertTrue(np.sum(gt)==2)
        self.assertTrue(len(gt), 3)

    def test_get_variants_info(self):
        self.testread.get_allele_freq('output.txt')
        dat = pd.read_table('output.txt', delim_whitespace=True)
        self.assertTrue(dat.shape[0] > 1)
        self.assertTrue(dat.shape[1] == 4)

    def test_sampling_variants(self):
        out = self.testread._sample_variants(10)
        dat = pd.read_csv(out, delim_whitespace=True, header=None)
        print(dat)
        self.assertTrue(dat.shape[0] == 10)
        self.assertTrue(dat.shape[1] == 4)

    def test_sampling_subjects(self):
        out = self.testread._sample_subjects(10)
        is_in = np.isin(out, self.testread._samples)
        self.assertTrue(len(out) == 10)
        self.assertTrue(np.all(is_in))

    def test_sample(self):
        self.testread.maf = 0.5
        gp = self.testread.sample(10,10)
        self.assertTrue(np.sum(np.sum(gp)) > 1)
        self.assertTrue(gp.shape[0] == 10)

    def test_genotype_loading(self):
        self.testread.maf = 0.5
        self.testread.get_allele_freq('output.txt')
        out = self.testread._sample_subjects(10)
        gp = self.testread.read_vcf(out, 'output.txt')
        self.assertTrue(gp.shape[0] == 10)
        self.assertTrue(gp.shape[1] > 10)
        self.assertTrue(np.sum(np.sum(gp)) > 10)

    def test_burden(self):
        def burden(cases, controls):
            return np.abs(np.sum(np.sum(cases)) - np.sum(np.sum(controls)))

        tests = {'burden': burden}
        self.testread.maf = 0.5
        # self.testread.get_allele_freq('output.txt')
        cases = self.testread._sample_subjects(10)
        controls = self.testread._sample_subjects(10)
        gp = self.testread.binary_test(cases, controls, tests, 'output.txt', 10)
        self.assertTrue(len(gp) > 1)
        self.assertTrue(gp['burden_p'] > 0)
        self.assertTrue(gp['burden_p'] < 1)
