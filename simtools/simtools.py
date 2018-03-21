import sys
import numpy as np
from scipy.optimize import root
from simtools.genotypes import ReadPlink


class Simtools(object):
    """Initiates a simtools object to generate various different phenotypes"""

    def __init__(self, plink_stem, subjects=None, p=None):
        """
        :param plink_stem: plink stem file path
        :param subjects: iid of subjects to use
        :param p: number of variants to use
        """
        self._plink_stem = plink_stem
        self._plink = ReadPlink(self._plink_stem)

        if subjects is None:
            self.subjects = self._plink.fam.iid.values
            self._id_subjects = self._plink.fam.index.values
            self.n = self._plink.fam.shape[0]
        else:
            self.n = len(subjects)
            self.subjects = subjects
            self._id_subjects = self._plink.fam.iid.get_loc(self.subjects)

        if p is None:
            self.p = self._plink.P
        else:
            self.p = p

        self.chunk_size = 100
        self.causal_snps = None
        self.liability_cases_controls = None
        self.last_random_subjects = None

    def _causal_SNPs(self, causal, weights=None):
        """Define causal SNPs

        :param causal: number, proportion or list of causal snps
        :param weights: how to weight causal SNPs (default is 1)
        :returns: list of causal SNPs and weights

        """
        causal_snps = []
        if isinstance(causal, int):
            causal_snps = np.random.choice(range(self.p), causal)
            causal_snps = self._plink.bim.index.values[causal_snps]
        if isinstance(causal, list):
            causal_snps = self._plink.bim.index.values[np.array(causal)]
        if isinstance(causal, float):
            n_causal = int(np.floor(self.p * causal))
            causal_snps = np.random.choice(range(self.p), n_causal)
            causal_snps = self._plink.bim.index.values[causal_snps]

        if weights == None:
            weights = np.ones(len(causal_snps))
        else:
            if len(weights) != len(self.causal_snps):
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
        """scales a matrix or vector.

        :param x: numpy matrix
        :returns: scaled numpy matrix

        """
        return (x - np.mean(x, axis=0)) / np.std(x, axis=0)

    def _compute_geffect(self, causal_snps, weights, subjects):
        """Computes genetic effect from a list of causal snps with weights for specific subjects

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
            temp_matrix = self._plink.read_bed(marker=snps,
                                               subjects=subjects)
            geffect += np.dot(temp_matrix, effect)

        return self._scale(geffect)

    def simple_phenotype(self, causal, hera, liability=None, n=None):
        """simulates a phenotypes (continues or binary)
        If a liability threshold is used, the method generates a binary phenotype

        :param causal: Number of causal SNPs
        :param hera: Heritability
        :param liability: Liability Threshold
        :param n: number of subjects to sample
        :returns: Vector of the phenotype

        """
        subjects = self._id_subjects
        if n is not None:
            subjects = np.random.choice(self._id_subjects, n)
            self.last_random_subjects = self.subjects[subjects]

        causal_snps, weights = self._causal_SNPs(causal)
        geffect = self._compute_geffect(causal_snps, weights, subjects)

        if liability is None:
            pheno = np.sqrt(hera) * geffect + np.sqrt(1 - hera) * np.random.normal(0,
                                                                                   np.sqrt(1), len(subjects))
            return pheno

        elif isinstance(liability, tuple) and len(liability) == 3:
            threshold = liability[0]
            ncases = liability[1]
            ncontrols = liability[2]

            cases, controls = self._liability_model(ncases, ncontrols,
                                                    threshold, hera, geffect)
            pheno = np.append(np.repeat(1, len(cases)), np.repeat(0, len(controls)))
            self.liability_cases_controls = np.append(cases, controls)
            self.last_random_subjects = self.subjects[self.liability_cases_controls]
            return pheno

        else:
            sys.exit('values are missing')

    def _liability_model(self, num_cases, num_controls,
                         threshold, hera, geffect, max_iter=10000):
        """Simulates cases and controls

        :param num_cases: number of cases
        :param num_controls: number of controls
        :param threshold: liability threshold
        :returns: subject IDs of cases and controls

        """
        container_cases = []
        container_controls = []

        for item in range(max_iter):
            pheno = np.sqrt(hera) * geffect + np.sqrt(1 - hera) * np.random.normal(0,
                                                                                   np.sqrt(1), self.n)
            if len(container_cases) < num_cases:
                container_cases.append(np.argwhere(pheno >= threshold))

            if len(container_controls) < num_controls:
                container_controls.append(np.argwhere(pheno < threshold))

            if (len(container_cases) >= num_cases and
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
        """Estimator function to estiamte scalar matrices

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
        """Computes multiple phenotype from estimates

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
        """Generates matrix of effects

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
        """Generates multiple inter-related phenotypes

        :param lamb: adjacency matrix
        :param B: genetic effect matrix
        :param num_causal: number of causal variants (randomly choosen)
        :param n: number of samples to randomly chose
        :returns: matrix of size n*t

        """
        if np.all(lamb.shape != B.shape):
            raise NameError("""dimensions of transition matrix and
                    genetic effect is not the same""")

        subjects = self._id_subjects

        if n is not None:
            subjects = np.random.choice(self._id_subjects, n)
            self.last_random_subjects = self.subjects[subjects]

        t = lamb.shape[0]

        Beta = self._multiple_effect_vectors(t, num_causal)
        overall_snp_effect = np.sum(Beta, 1)
        causal_snps = np.where(overall_snp_effect > 0)
        Beta = Beta[causal_snps]
        causal_snps = self._plink.bim.index.values[causal_snps[0]]

        G = self._compute_geffect(causal_snps, Beta, subjects)

        sol = root(lambda k: self._estimateF(k, G=G, lamb=lamb, B=B),
                   np.array([0.3, 0.3, 0.3]), method='krylov')
        phenotypes = self._compute_multi_pheno(sol.x, G, lamb, B)
        return phenotypes
