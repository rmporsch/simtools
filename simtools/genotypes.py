#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import re
import vcf
import os
from pandas_plink import read_plink

def simple_genotype_matrix(n, p, min_maf=0.05, max_maf=0.5):
    """Generates a simple matrix containing either 0 or 1 of size nxp

    :n: number of samples
    :p: number of genotypes
    :min_maf: min frequency
    :max_maf: max frequency
    :returns: a numpy matrix of size nxp

    """
    genotypes = np.zeros(shape=(n, p))
    for item in range(0, p):
        genotypes[:,item] = np.random.binomial(1, np.random.uniform(0.1, 0.5, 1), n)

    return genotypes


class ReadVCF(object):

    """Reads a given vcf file"""

    def __init__(self, vcffile):
        """
        :vcffile: path to the vcf file

        """
        self._vcffile = vcffile
        
        
    def read_vcf(self, index_file=False):
        """Reads the vcf file and gets and index

        :index_file: location of the index file

        """
        vcf_reader = vcf.Reader(filename=self._vcffile)
        self._samples = vcf_reader.samples

        temp_folder = '.temp'
        temp_file = '.temp/variant.list'

        if not os.path.exists(temp_folder):
            os.makedirs(temp_folder)

        if not os.path.isfile(temp_file) or index_file:
            with open(temp_file, 'w') as f:
                for record in vcf_reader:
                    f.write('%s %i %s %f\n' %
                            (record.CHROM, record.POS, record.ID, record.INFO['AF']))


    def __sample_variants(self, maf, p, write_disk=False):
        """Sample a random set of variants

        :maf: minor allele frequency cutoff
        :p: number of variants to sample

        """
        self._variants = pd.read_table('./.temp/variant.list',
                header=None, sep=' ',
                names=['Chrom', 'Pos', 'ID', 'AF'])

        self._variants = self._variants.query('AF >= @maf and AF <= (1- @maf)')
        self._sampled = self._variants.sample(n=p)
        

        if write_disk:
            self._sampled.to_csv('.temp/subsample_variant.list')

    def __sample_subjects(self, n):
        """Sample a random set of subjects

        :n: number of subjects to sample

        """
        self._randomset = np.random.choice(self._samples, n, replace=True)

    def sample(self, n, p, maf=0.05, write_disk=False):
        """Random sample a set of variants and subjects

        :maf: minor allele frequency cutoff
        :n: number of subjects to sample
        :p: number of variants to sample
        :write_disk: bool, write to disk a list of variants
        :returns: a numpy matrix of size n*p

        """
        self.__sample_variants(maf, p, write_disk)
        self.__sample_subjects(n)

        genotypematrix = np.zeros((n,p))
        reader = vcf.Reader(filename=self._vcffile)
        i = 0
        for index, row in self._sampled.iterrows():
            record = reader.fetch(str(row[0]), row[1]-1, row[1])
            for rr in record:
                flip = True if rr.INFO['AF'] > 0.5 else False
                u = 0
                for sample in rr.samples:
                    if sample.sample in self._randomset:
                        genotypematrix[u,i] = self.__count_genotypes(sample['GT'], flip)
                        u += 1
            i += 1

        return genotypematrix

    def __count_genotypes(self, rr, flip):
        """counts the genotypes at a given position
        !!! missing is counted as 0 !!!

        :rr: genotype string from vcf object
        :returns: 0, 1 or 2

        """
        genotype = np.array(list(map(int, re.findall(r'(\d+)', rr))))
        if flip:
            genotypes = abs(genotype - 1)
        return sum(genotype)



class ReadPlink(object):

    """Reads plink files"""

    def __init__(self, plinkstem):
        """
        :plinkstem: plink stem file path

        """
        self._plinkstem = plinkstem
        (self.bim, self.fam, self.G) = read_plink(plinkstem, verbose=False)


    def sample(self, n, p, write_disk=False):
        """Samples from a plink file with random SNPs and subjects
        Currently pandas_plink does not support fancy indexing, hence
        sample will load the genotypes of all subjects before randomly sample
        subjects IDs.

        :maf: minor allele frequency cutoff
        :n: number of subjects to sample
        :p: number of variants to sample
        :write_disk: bool, write to disk a list of variants
        :returns: a numpy matrix of size n*p

        """
        self.__sample_subjects = np.random.choice(self.fam.i.values, n, replace=True)
        self.__sample_variants = np.random.choice(self.bim.i.values, p)

        if write_disk:
            self.bim.iloc[self.__sample_variants].to_csv('sampled_variants.csv')

        genotypematrix =  self.G[self.__sample_variants, :].compute()
        genotypematrix = genotypematrix[:, self.__sample_subjects]
        np.nan_to_num(genotypematrix, copy=False)

        return genotypematrix.T
