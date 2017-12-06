import pandas as pd
import numpy as np
import sys
from sklearn.preprocessing import scale
from scipy.optimize import root
import pymp
import statsmodels.api as sm


class Simtools(object):

    """Containing multiple methods to simulate phenotypes"""

    def __init__(self, matrix):
        """Initiate a simtools object

        :matrix numpy matrix of genotypes
        """
        self.genotypematrix = matrix
        self.n = matrix.shape[0]
        self.p = matrix.shape[1]
        self.index = np.arange(0, self.n)
        self.geno_index = np.arange(0, self.p)


    def simple_phenotype(self, causal, hera, liability=None):
        """simulates a phenotypes (continues or binary)
        If a liability threshold is used, the method generates a binary phenotype

        :causal Number of causal SNPs
        :hera: Heritability
        :liability: Liability Threshold
        :returns: Vector of the phenotype

        """

        self.causal = self.__define_causal(causal)
        geffect = scale(np.dot(self.genotypematrix, self.causal))

        if liability is None:
            pheno = np.sqrt(hera)*geffect + np.sqrt(1-hera)*np.random.normal(0,
                    np.sqrt(1), self.n)
            return pheno

        elif isinstance(liability, tuple) and len(liability)==3:
            threshold = liability[0]
            ncases = liability[1]
            ncontrols = liability[2]

            cases, controls = self.__liability_model(ncases, ncontrols,
                    threshold, hera, geffect)
            pheno = np.append(np.repeat(1, len(cases)), np.repeat(0, len(controls)))
            self.index = np.append(cases, controls)
            return pheno

        else:
            sys.exit('values are missing')


    def __define_causal(self, causal):
        """Simulates a causal vector

        :causal: Number of causal SNPs
        :returns: Vector of Positions of causal SNPs

        """

        if type(causal) is float:
            causal = np.random.binomial(1, causal, self.p)
            return causal

        if type(causal) is np.ndarray:
            if ((any(x not in set(causal) for x in [0,1])) or
                    (len(set(causal))>2) or
                    len(causal)!=self.p):
                sys.exit('invalid causal vector')
            else:
                return causal


    def __liability_model(self, num_cases, num_controls,
            threshold, hera, geffect, max_iter=10000):
        """simulates cases and controls

        :num_cases: number of cases
        :num_controls: number of controls
        :threshold: liability threshold
        :returns: subject IDs of cases and controls

        """
        container_cases = []
        container_controls = []

        for item in range(max_iter):
            pheno = np.sqrt(hera)*geffect + np.sqrt(1-hera)*np.random.normal(0,
                    np.sqrt(1), self.n)
            if len(container_cases) < num_cases:
                container_cases.append(np.argwhere(pheno >=threshold))

            if len(container_controls) < num_controls:
                container_controls.append(np.argwhere(pheno <threshold))

            if (len(container_cases)>=num_cases and
                    len(container_controls) >=num_controls):
                break

        # unlist
        container_cases = [item for sublist in container_cases for item in sublist]
        container_controls = [item for sublist in container_controls for item in sublist]

        # remove overhanging samples
        container_cases = container_cases[0:(num_cases)]
        container_controls = container_controls[0:(num_controls)]

        return container_cases, container_controls

    def estimateF(self, x, G, lamb, B):
        """estimator function to estiamte scalar matrices

        :x: value to find
        :G: genetic effect matrix
        :lamb: adjacency matrix
        :B: genetic effect matrix
        :returns: variance matrix of phenotypes

        """
        t = lamb.shape[0]
        invert = np.linalg.inv(np.identity(t) - lamb)
        D = np.diag(x)
        C = invert.dot(B); F = invert.dot(D)

        outcov = C.dot(np.cov(G.T)).dot(C.T) + F.dot(np.identity(t)).dot(F.T)
        return outcov.diagonal() - np.ones(t)

    def __compute_multi_pheno(self, x, G, lamb, B):
        t = lamb.shape[0]
        error = np.zeros((t, G.shape[0]))
        for i in range(t):
            error[i, :] = np.random.normal(scale=np.sqrt(1), size=G.shape[0])
        invert = np.linalg.inv(np.identity(t) - lamb)
        D = np.diag(x)
        temp = (B.dot(G.T) + D.dot(error))
        return invert.dot(temp)

    def __multiple_effect_vectors(self, t, num_causal_snps):
        Beta = np.zeros([self.p, t])
        causal_index = [np.random.randint(low=0, high=self.p, size=k)
                for k in num_causal_snps]
        # simualte effect sizes for each variant
        for i in range(Beta.shape[1]):
            Beta[causal_index[i],i] = 1

        return Beta

    def multi_phenotype(self, lamb, B, num_causal):
        """Generates multiple inter-related phenotypes

        :lamb: adjacency matrix
        :B: genetic effect matrix
        :num_causal: number of causal variants (randomly choosen)
        :returns: matrix of size n*t

        """
        if np.all(lamb.shape!=B.shape):
            raise NameError("""dimensions of transition matrix and
                    genetic effect is not the same""")

        t = lamb.shape[0]
        Beta = self.__multiple_effect_vectors(t, num_causal)
        G = self.genotypematrix.dot(Beta)
        G = (G - G.mean(axis=0)) / G.std(axis=0)

        sol = root(lambda k:self.estimateF(k, G=G, lamb=lamb, B=B),
             np.array([0.3, 0.3,0.3]), method='krylov')
        phenotypes = self.__compute_multi_pheno(sol.x, G, lamb, B)
        return phenotypes
