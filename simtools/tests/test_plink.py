import unittest
import numpy as np
import pandas as pd
import os
from pandas_plink import example_file_prefix
from simtools import genotypes as gp

DATA = os.path.join(os.path.dirname(__file__), 'data')

class TestPlink(unittest.TestCase):

    def setUp(self):
        self.n = 1000
        self.p = 100
        self.ff = "/home/robert/Documents/projects/prs/data/genotypes/subset_10k"
        self.testread = gp.ReadPlink(self.ff)

    def test_sampling(self):
        temp = self.testread.sample(self.n, self.p)
        self.assertEqual(temp.shape[0], self.n, 'incorrect row numbers')
        self.assertEqual(temp.shape[1], self.p, 'incorrect col numbers')


if __name__ == '__main__':
    unittest.main()

