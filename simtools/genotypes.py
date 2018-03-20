#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import re
import vcf
import os
from pyplink import PyPlink
from tqdm import tqdm



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
        self.maf = 0.01

    def _get_samples(self):
        reader = vcf.Reader(filename=self._vcf_file)
        self._samples = reader.samples

    def _get_genotypes(self, samples, records, switch):
        """
        Gets the genotypes from records

        :param samples: list of subject IDs
        :param records: record object
        :param switch:  switch minor major allele *bool*
        :return: list of genotypes
        """

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
        """
        Criteria to filter out record

        :param record: record objects
        :return: two values *bool*; process (pass filter), switch (switch major/minor allele)
        """

        process = True
        switch = False
        if record.aaf[0] > 0.5:
            switch = True
        return process, switch

    def _check_samples(self, samples):
        """
        Check if samples are in vcf file

        :param samples: list of samples
        :return: samples in vcf file
        """
        check = [k not in self._samples for k in samples]
        num_not_in_vcf = np.sum(check)
        if num_not_in_vcf > 0:
            with open('/tmp/excluded.subjects', 'w') as f:
                for item in samples[np.array(check)]:
                    f.write("%s\n" % item)
            print(num_not_in_vcf, 'were not in vcf file and were removed; see /tmp/excluded.subjects')
        return samples[~np.array(check)]

    def load_genotype_matrix(self, subjects, variants):
        """
        Load genotype matrix

        :param subjects: list of subjects
        :param variants: path to variant file *str*; requires at least 2 columns [CHR, BP]; tab separated
        :return: matrix of genotypes
        """
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
            for record in vcf_reader.fetch(str(row[0]), row[1]-1, row[1]):
                process, switch = self.criteria(record)
                if not process:
                    skiped_variants.append(index)
                    continue
                matrix[:, index] = self._get_genotypes(subjects, record, switch)

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
                        (record.CHROM, record.POS, record.ID, record.aaf[0]))

    def _sample_variants(self, p, file_path=None):
        """Sample a random set of variants

        :maf: minor allele frequency cutoff
        :p: number of variants to sample

        """
        output_path = '/tmp/subsample_variant.list'
        if file_path is None:
            file_path = '/tmp/variant.list'
            if os.path.isfile(file_path):
                self._variants = pd.read_table(file_path,
                                               header=None, sep=' ',
                                               names=['Chrom', 'Pos', 'ID', 'AF'])
            else:
                self.get_allele_freq(file_path)
                self._variants = pd.read_table(file_path,
                                               header=None, sep=' ',
                                               names=['Chrom', 'Pos', 'ID', 'AF'])
        else:
            self._variants = pd.read_table(file_path, header=None, sep=' ')

        self._sampled = self._variants.sample(n=p)
        self._sampled.to_csv(output_path, sep=' ', index=False, header=False)
        return output_path

    def _sample_subjects(self, n):
        """Sample a random set of subjects

        :n: number of subjects to sample

        """
        return np.random.choice(self._samples, n, replace=True)

    def sample(self, n, p):
        """Random sample a set of variants and subjects

        :maf: minor allele frequency cutoff
        :n: number of subjects to sample
        :p: number of variants to sample
        :returns: a numpy matrix of size n*p

        """
        variant_file = self._sample_variants(p)
        subjects = self._sample_subjects(n)

        genotypematrix = self.load_genotype_matrix(subjects, variant_file)
        return genotypematrix

    def binary_test(self, cases, controls, tests, variants, iteration):
        """
        Computes p-values for a binary phenotype for various different tests

        :param cases: list of case names
        :param controls: list of control names
        :param tests: dict of tests
        :param variants: path to variant files
        :param iteration: number of iterations to compute p-values
        :return: dict with the results
        """

        cases = self.load_genotype_matrix(cases, variants)
        controls = self.load_genotype_matrix(controls, variants)

        n_cases = len(cases)
        n_controls = len(controls)

        output = {'n_cases': n_cases,
                  'n_controls': n_controls,
                  'n_var': cases.shape[1]}
        null = {}

        for i, (func_name, func) in enumerate(tests.items()):
            output[func_name] = func(cases, controls)
            null[func_name] = []

        genotype_matrix = np.vstack((cases, controls))
        del cases
        del controls
        n_subjects = genotype_matrix.shape[0]

        for i in tqdm(range(iteration)):
            temp = genotype_matrix
            np.random.shuffle(temp)
            for u, (func_name, func) in enumerate(tests.items()):
                null[func_name].append(func(temp[0:n_cases, :], temp[n_cases:n_subjects, :]))

        for i, (func_name, func) in enumerate(tests.items()):
            output[func_name+'_p'] = np.mean(null[func_name] >= output[func_name])

        return output


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
