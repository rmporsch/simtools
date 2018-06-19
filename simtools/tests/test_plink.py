import unittest
from simtools.read_plink import ReadPlink


class TestPlink(unittest.TestCase):

    def setUp(self):
        self.ff = "data/1kg_phase1_chr22"
        self.testread = ReadPlink(self.ff)

    def test_genotypematrix(self):
        n = 100
        p = 1000
        temp = self.testread.get_gentoypematrix(
            self.testread.variants[0:p],
            self.testread.subject[0:n])
        self.assertEqual(temp.shape[0], n, 'incorrect row numbers')
        self.assertEqual(temp.shape[1], p, 'incorrect col numbers')


if __name__ == '__main__':
    unittest.main()
