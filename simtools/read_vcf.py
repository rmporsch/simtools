import pandas as pd
import numpy as np
import re
import vcf
import os


class ReadVCF(object):
    """Reads VCF File."""

    def __init__(self, vcffile):
        """Init function for vcf file.

        :param vcffile: path to the vcf file
        """
        self.subjects = None
        self.variants = np.array([])
        self.P = None
        self.N = None
        self.fam = None
        self._vcf_file = vcffile
        self._get_samples()
        self._get_allele_freq()
        self.N = len(self.subjects)

    def _get_samples(self):
        reader = vcf.Reader(filename=self._vcf_file)
        self.subjects = np.array(reader.samples)
        self.fam = pd.DataFrame({'iid': self.subjects})

    def _get_genotypes(self, samples, records, switch):
        """Get the genotypes from records.

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
                print(
                    'sample:', sample,
                    'variant:', records,
                    '-- set value to missing'
                    )
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

    def _criteria(self, record):
        """Criteria to filter out record.

        :param record: record objects
        :return: process (pass filter), switch (switch major/minor allele)
        """
        process = True
        switch = False
        if record.aaf[0] > 0.5:
            switch = True
        return process, switch

    def _check_samples(self, samples):
        """Check if samples are in vcf file.

        :param samples: list of samples
        :return: samples in vcf file
        """
        samples = np.array(samples)
        check = np.isin(samples, self.subjects, invert=True)
        num_not_in_vcf = np.sum(check)
        exclude_file_path = '.excluded.subjects'
        if num_not_in_vcf > 0:
            with open(exclude_file_path, 'w') as f:
                for item in samples[np.array(check)]:
                    f.write("%s\n" % item)
            print(
                num_not_in_vcf,
                'not in file and were removed; see %s' % exclude_file_path)
            return samples[~check]
        elif num_not_in_vcf >= len(samples):
            raise RuntimeError('No subjects is present in vcf file')
        else:
            return samples

    def get_gentoypematrix(self, variants=None, subjects=None):
        """Load genotype matrix.

        :param subjects: list of subjects
        :param variants: list of variants (format: chr1:234234234)
        :return: matrix of genotypes
        """
        if variants is None:
            p_size = self.P
            variants = np.array(self.variants)
        else:
            p_size = len(variants)
            variants = np.array(variants)

        if subjects is None:
            n_size = self.N
            subjects = self.subjects
        else:
            n_size = len(subjects)
            if not isinstance(subjects[0], str):
                subjects = self.subjects[subjects]

        subjects = self._check_samples(subjects)
        variants = pd.DataFrame([k.split(':') for k in variants])
        variants.iloc[:, 1] = pd.to_numeric(variants.iloc[:, 1])
        matrix = np.zeros((n_size, p_size))
        skiped_variants = []

        vcf_reader = vcf.Reader(filename=self._vcf_file)
        for index, row in variants.iterrows():
            for record in vcf_reader.fetch(str(row[0]), row[1]-1, row[1]):
                process, switch = self._criteria(record)
                if not process:
                    skiped_variants.append(index)
                    continue
                matrix[:, index] = self._get_genotypes(
                    subjects, record, switch)

        # remove skiped variants
        matrix = np.delete(matrix, skiped_variants, 1)
        return matrix

    def _get_allele_freq(self, output_path='.vcf_variants.bim'):
        """Read the vcf file and gets and index.

        :param index_file: location of the index file
        """
        counter = 0
        if os.path.isfile(output_path):
            with open(output_path, 'r') as f:
                for line in f:
                    line = line.split()
                    self.variants = np.append(self.variants,
                                              (str(line[0]+':'+str(line[1]))))

                    counter += 1
            self.P = counter
        else:
            print("Scanning vcf file and writing stats to %s" % output_path)
            vcf_reader = vcf.Reader(filename=self._vcf_file)
            with open(output_path, 'w') as f:
                for record in vcf_reader:
                    f.write('%s %i %s %f\n' %
                            (record.CHROM, record.POS, record.ID, record.aaf[0]))
                    self.variants = np.append(self.variants,
                                              str(record.CHROM)+':'+str(record.POS))
                    counter += 1
            self.P = counter
