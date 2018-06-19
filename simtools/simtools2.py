"""Simtool Class.

The Simtool class enables to simulate based on data given from a vcf file.
or plink file.
"""
import numpy as np
import pandas as pd
import re
from read_vcf import ReadVCF
from read_plink import ReadPlink
from scipy.optimize import root
import os


class Simtools(object):
    """docstring for Simtools."""

    def __init__(self, genotype_path, chunk_size):
        """Init for Simtools."""
        assert os.path.isfile(genotype_path)
        self.genotype_path = genotype_path
        assert isinstance(genotype_path, str)
        self.type = None
        if bool(re.search('vcf', genotype_path)):
            print('Assume file in vcf format')
            self.type = 'vcf'
            self._reader = ReadVCF(genotype_path)
            self.get_genotypematrix = self._reader.get_genotypematrix
        else:
            print('Assume file in plink format')
            self.type = 'plink'
            self._reader = ReadPlink(genotype_path)
            self.get_genotypematrix = self._reader.get_genotypematrix
        self.subject = self._reader.subject
        self.variants = self._reader.variants
        self.N = self._reader.N
        self.P = self._reader.P
        self.chunk_size = chunk_size

    def sample_genotypematrix(self, n, p):
        """Sample from genotype matrix.

        Samples from a plink file with random SNPs and subjects
        Currently pandas_plink does not support fancy indexing, hence
        sample will load the genotypes of all subjects before randomly sample
        subjects IDs.

        :param n: number of subjects to sample
        :param p: number of variants to sample
        :returns: a numpy matrix of size n*p

        """
        assert isinstance(n, int)
        assert isinstance(p, int)
        self._sample_subjects = np.random.choice(self.subject, n, replace=True)
        self._sample_variants = np.random.choice(self.variants, p)

        genotypematrix = self.get_gentoypematrix(
            self._sample_variants, self._sample_subjects)

        return genotypematrix

    def _causal_SNPs(self, causal, weights=None):
        """Define causal SNPs.

        :param causal: number, proportion or list of causal snps
        :param weights: how to weight causal SNPs (default is 1)
        :returns: list of causal SNPs and weights

        """
        causal_snps = []
        if isinstance(causal, int):
            causal_snps = np.random.choice(range(self.P), causal)
            causal_snps = self.variants[causal_snps]
        if isinstance(causal, list):
            causal_snps = self.variants[np.array(causal)]
        if isinstance(causal, float):
            n_causal = int(np.floor(self.P * causal))
            causal_snps = np.random.choice(range(self.P), n_causal)
            causal_snps = self.variants[causal_snps]

        if weights is None:
            weights = np.ones(len(causal_snps))
        else:
            if len(weights) != len(causal_snps):
                ValueError('Number of weights and number of SNPs do not agree')

        self.causal_snps = causal_snps
        return causal_snps, weights

    def _chunks(self, l, n):
        """Yield successive n-sized chunks from l.

        :param l: list of things
        :param n: number of chunks
        """
        for i in range(0, len(l), n):
            yield l[i:i + n]

    def _scale(self, x):
        """Scales a matrix or vector.

        :param x: numpy matrix
        :returns: scaled numpy matrix

        """
        return (x - np.mean(x, axis=0)) / np.std(x, axis=0)

    def _compute_geffect(self, causal_snps, weights, subjects):
        """Compute the genetic effect.

        Computes genetic effect from a list of causal snps
        with weights for specific subjects.

        :param causal_snps: array of causal SNPs
        :param weights: array of weights for the SNPs
        :param subjects: list of subjects
        :returns: vector of size n containing the genetic effect
        """
        snp_chunks = self._chunks(causal_snps, self.chunk_size)
        effect_chunks = self._chunks(weights, self.chunk_size)
        t = weights.shape
        n = len(subjects)

        # check number of phenotypes
        if len(t) > 1:
            t = t[1]
            geffect = np.zeros((n, t))
        else:
            t = 1
            geffect = np.zeros((n,))

        for snps, effect in zip(snp_chunks, effect_chunks):
            temp_matrix = self.get_gentoypematrix(
                marker=snps, subjects=subjects)
            geffect += np.dot(temp_matrix, effect)

        return self._scale(geffect)

    def simple_phenotype(self, causal, hera, liability=None, n=None):
        """Simulate a phenotypes (continues or binary).

        If a liability threshold is used, method generates a binary phenotype
        :param causal: Number of causal SNPs
        :param hera: Heritability
        :param liability: Liability Threshold
        :param n: number of subjects to sample
        :returns: Vector of the phenotype

        """
        assert isinstance(hera, float)
        assert hera <= 1.0
        assert n <= self.N
        subjects = np.arange(0, self.N)
        subject_names = self.subject
        if n is not None:
            subjects = np.random.choice(np.arange(0, self.N), n)
            subject_names = self.subjects[subjects]

        causal_snps, weights = self._causal_SNPs(causal)
        geffect = self._compute_geffect(causal_snps, weights, subjects)

        if liability is None:
            pheno = np.sqrt(hera) * geffect + np.sqrt(1 - hera) *\
                 np.random.normal(0, np.sqrt(1), len(subjects))
            return pd.DataFrame({'IID': subject_names, 'Pheno': pheno})

        else:
            assert isinstance(liability, tuple)
            assert len(liability) == 3
            threshold = liability[0]
            ncases = liability[1]
            ncontrols = liability[2]
            cases, controls = self._liability_model(ncases, ncontrols,
                                                    threshold, hera, geffect)
            pheno = np.append(
                np.repeat(1, len(cases)), np.repeat(0, len(controls)))
            id_cases_controls = np.append(cases, controls)
            subject_names = self.subjects[id_cases_controls]
            return pd.DataFrame({'IID': subject_names, 'Pheno': pheno})

    def _liability_model(self, num_cases, num_controls,
                         threshold, hera, geffect, max_iter=10000):
        """Simulate cases and controls.

        :param num_cases: number of cases
        :param num_controls: number of controls
        :param threshold: liability threshold
        :returns: subject IDs of cases and controls

        """
        container_cases = []
        container_controls = []

        for item in range(max_iter):
            pheno = np.sqrt(hera) *\
                 geffect + np.sqrt(1 - hera) *\
                 np.random.normal(0, np.sqrt(1), self.n)
            if len(container_cases) < num_cases:
                container_cases.append(np.argwhere(pheno >= threshold))

            if len(container_controls) < num_controls:
                container_controls.append(np.argwhere(pheno < threshold))

            if (len(container_cases) >= num_cases) and (
                len(container_controls) >= num_controls):
                break

        # unlist
        container_cases = [item for sublist in container_cases for item in sublist]
        container_controls = [item for sublist in container_controls for item in sublist]

        # remove overhanging samples
        container_cases = container_cases[0:num_cases]
        container_controls = container_controls[0:num_controls]

        return container_cases, container_controls

    def _estimateF(self, x, G, lamb, B):
        """Estimator function to estiamte scalar matrices.

        :param x: value to find
        :param G: genetic effect matrix
        :param lamb: adjacency matrix
        :param B: genetic effect matrix
        :returns: variance matrix of phenotypes

        """
        t = lamb.shape[0]
        invert = np.linalg.inv(np.identity(t) - lamb)
        D = np.diag(x)
        C = invert.dot(B)
        F = invert.dot(D)

        outcov = C.dot(np.cov(G.T)).dot(C.T) + F.dot(np.identity(t)).dot(F.T)
        return outcov.diagonal() - np.ones(t)

    def _compute_multi_pheno(self, x, G, lamb, B):
        """Compute multiple phenotype from estimates.

        :param x: optimal scaling matrix
        :param G: genetic effect
        :param lamb: adjacency matrix
        :param B: genetic cross effect matrix
        :return: phenotypes

        """
        t = lamb.shape[0]
        error = np.zeros((t, G.shape[0]))
        for i in range(t):
            error[i, :] = np.random.normal(scale=np.sqrt(1), size=G.shape[0])
        invert = np.linalg.inv(np.identity(t) - lamb)
        D = np.diag(x)
        temp = (B.dot(G.T) + D.dot(error))
        return invert.dot(temp)

    def _multiple_effect_vectors(self, t, num_causal_snps):
        """Generate matrix of effects.

        :param t: number of phenotypes
        :param num_causal_snps: number of causal snps for each phenotype
        :return: matrix of causal effects

        """
        Beta = np.zeros([self.p, t])
        causal_index = [np.random.randint(low=0, high=self.p, size=k)
                        for k in num_causal_snps]
        # simualte effect sizes for each variant
        for i in range(Beta.shape[1]):
            Beta[causal_index[i], i] = 1

        return Beta

    def multi_phenotype(self, lamb, B, num_causal, n=None):
        """Generate multiple inter-related phenotypes.

        :param lamb: adjacency matrix
        :param B: genetic effect matrix
        :param num_causal: number of causal variants (randomly choosen)
        :param n: number of samples to randomly chose
        :returns: matrix of size n*t

        """
        assert isinstance(n, int)
        if np.all(lamb.shape != B.shape):
            raise NameError("""dimensions of transition matrix and
                    genetic effect is not the same""")

        subjects = np.arange(0, self.N)

        if n is not None:
            subjects = np.random.choice(self.subjects, n)
            subject_names = self.subjects[subjects]

        t = lamb.shape[0]

        Beta = self._multiple_effect_vectors(t, num_causal)
        overall_snp_effect = np.sum(Beta, 1)
        causal_snps = np.where(overall_snp_effect > 0)
        Beta = Beta[causal_snps]
        causal_snps = self.variants[causal_snps[0]]

        G = self._compute_geffect(causal_snps, Beta, subjects)

        sol = root(lambda k: self._estimateF(k, G=G, lamb=lamb, B=B),
                   np.array([0.3, 0.3, 0.3]), method='krylov')
        phenotypes = self._compute_multi_pheno(sol.x, G, lamb, B)
        out = pd.DataFrame({'IID': subject_names})
        out = pd.concat([out, pd.DataFrame(phenotypes)])
        return out
