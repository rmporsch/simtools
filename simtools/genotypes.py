#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import re
import vcf
import os
from pyplink import PyPlink


def simple_genotype_matrix(n, p):
    """Generates a simple matrix containing either 0 or 1 of size nxp

    :n: number of samples
    :p: number of genotypes
    :min_maf: min frequency
    :max_maf: max frequency
    :returns: a numpy matrix of size nxp

    """
    genotypes = np.zeros(shape=(n, p))
    for item in range(0, p):
        genotypes[:, item] = np.random.binomial(1, np.random.uniform(0.1, 0.5, 1), n)

    return genotypes


class ReadVCF(object):

    """Reads a given vcf file"""

    def __init__(self, vcffile):
        """
        :vcffile: path to the vcf file

        """
        self._vcf_file = vcffile
        self._samples = None
        self._randomset = None
        self._get_samples()

    def _get_samples(self):
        reader = vcf.Reader(filename=self._vcf_file)
        self._samples = reader.samples

    def _get_genotypes(self, samples, records, switch):
        variant = np.zeros(len(samples))
        for idx, sample in enumerate(samples):
            try:
                gt = records.genotype(sample)['GT']
            except IndexError:
                print("something went wrong with:")
                print('sample:', sample, 'variant:', records, '-- set value to missing')
                gt = '.'
            if gt == '.':
                gt = 0
            else:
                gt = re.split('\||/', gt)
                gt = list(map(int, gt))
            variant[idx] = np.sum(gt)
        if switch:
            variant = np.abs(variant - 2)
        return variant

    def criteria(self, record):
        process = True
        switch = False
        if record.call_rate <= 0.9:
            process = False
        if not record.is_snp:
            process = False
        if len(record.aaf) > 1:
            process = False

        if record.aaf[0] > 0.5:
            if (1-record.aaf[0]) > 0.01:
                process = False
            switch = True
        else:
            if record.aaf[0] > 0.01:
                process = False
        return process, switch

    def _check_samples(self, samples):
        check = [k not in self._samples for k in samples]
        num_not_in_vcf = np.sum(check)
        if num_not_in_vcf > 0:
            print(num_not_in_vcf, 'were not in vcf file and were removed')
        return samples[~np.array(check)]

    def load_genotype_matrix(self, subjects, variants):
        if not os.path.isfile(variants):
            raise NameError('File does not exist')
        subjects = self._check_samples(subjects)
        variants = pd.read_table(variants, header=None, sep=' ')
        if variants.shape[1] < 2:
            raise NameError('Variant file does not seem to have the right amount of columns')
        matrix = np.zeros((len(subjects), variants.shape[0]))
        skiped_variants = []

        vcf_reader = vcf.Reader(filename=self._vcf_file)
        for index, row in variants.iterrows():
            record = vcf_reader.fetch(str(row[0]), row[1], row[1])
            process, switch = self.criteria(record)
            if not process:
                skiped_variants.append(index)
                continue
            matrix[:, index] = self._get_genotypes(samples, record, switch)

        # remove skiped variants
        matrix = np.delete(matrix, skiped_variants, 1)
        return matrix

    def get_allele_freq(self, output_path):
        """Reads the vcf file and gets and index

        :index_file: location of the index file

        """
        vcf_reader = vcf.Reader(filename=self._vcf_file)
        with open(output_path, 'w') as f:
            for record in vcf_reader:
                f.write('%s %i %s %f\n' %
                        (record.CHROM, record.POS, record.ID, record.aaf))

    def __sample_variants(self, p, write_disk=False):
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

    def sample(self, n, p, write_disk=False):
        """Random sample a set of variants and subjects

        :maf: minor allele frequency cutoff
        :n: number of subjects to sample
        :p: number of variants to sample
        :write_disk: bool, write to disk a list of variants
        :returns: a numpy matrix of size n*p

        """
        self.__sample_variants(p, write_disk)
        self.__sample_subjects(n)

        genotypematrix = np.zeros((n,p))
        reader = vcf.Reader(filename=self._vcf_file)
        i = 0
        for index, row in self._sampled.iterrows():
            record = reader.fetch(str(row[0]), row[1], row[1])
            for rr in record:
                flip = True if rr.aaf > 0.5 else False
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

    """Reads plink files and allows random sampling"""

    def __init__(self, plinkstem):
        """
        :plinkstem: plink stem file path

        """
        self._plinkstem = plinkstem
        self._bim_path = os.path.basename(self._plinkstem)+'.bim'
        self._bed_path = os.path.basename(self._plinkstem)+'.bed'
        self._fam_path = os.path.basename(self._plinkstem)+'.fam'

        self.plinkfile = PyPlink(self._plinkstem)
        self.fam = self.plinkfile.get_fam()
        self.bim = self.plinkfile.get_bim()
        self.N = self.fam.shape[0]
        self.P = self.bim.shape[0]


    def sample(self, n, p, write_disk=False):
        """Samples from a plink file with random SNPs and subjects
        Currently pandas_plink does not support fancy indexing, hence
        sample will load the genotypes of all subjects before randomly sample
        subjects IDs.

        :n: number of subjects to sample
        :p: number of variants to sample
        :write_disk: bool, write to disk a list of variants
        :returns: a numpy matrix of size n*p

        """
        self.__sample_subjects = np.random.choice(self.fam.index.values, n, replace=True)
        self.__sample_variants = np.random.choice(self.bim.index.values, p)

        if write_disk:
            self.bim.iloc[self.__sample_variants].to_csv('sampled_variants.csv')
            self.bim.iloc[self.__sample_subjects].to_csv('sampled_subjects.csv')

        genotypematrix =  self.read_bed(self.__sample_variants,
                self.__sample_subjects)

        return genotypematrix

    def read_bed(self, marker=None, subjects=None):
        """read bed file

        :marker: list of SNPs
        :subjects: list of subjects
        :returns: genotypematrix of size subjects*marker

        """
        if marker is None:
            P_size = self.P
            marker = self.bim.index.values
        else:
            P_size = len(marker)

        if subjects is None:
            N_size = self.N
            subjects = self.fam.index.values
        else:
            N_size = len(subjects)

        genotypematrix = np.zeros((N_size, P_size), dtype=np.int8)

        j = 0
        for m, g in self.plinkfile.iter_geno_marker(marker):
            genotypematrix[:,j] = g[subjects]
            j += 1

        genotypematrix[genotypematrix < 0] = 0

        return genotypematrix
