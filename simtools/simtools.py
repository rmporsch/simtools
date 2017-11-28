import pandas as pd
import numpy as np
import sys
from sklearn.preprocessing import scale
import pymp
import statsmodels.api as sm



class Simtools(object):

    """Containing multiple functions to simulate phenotypes"""

    def __init__(self, matrix):
        """Initiate a simtools object
        
        :matrix pandas DataFrame object with the genotypes
        """
        self.genotypematrix = matrix
        self.n = matrix.shape[0]
        self.p = matrix.shape[1]
        self.index = np.arange(0, self.n)
        self.geno_index = np.arange(0, self.p)
            

    def simple_phenotype(self, causal, hera, liability=None):
        """simulates a phenotypes (continues or binary)

        :causal Number of causal SNPs
        :hera: Heritability
        :liability: Liability Threshold
        :returns: Vector of the phenotype

        """

        self.causal = self.define_causal(causal)
        geffect = scale(np.dot(self.genotypematrix, self.causal))

        if liability is None:
            pheno = np.sqrt(hera)*geffect + np.sqrt(1-hera)*np.random.normal(0, np.sqrt(1), self.n)
            return pheno

        elif isinstance(liability, tuple) and len(liability)==3:
            threshold = liability[0]
            ncases = liability[1]
            ncontrols = liability[2]

            cases, controls = self.liability_model(ncases, ncontrols, threshold, hera, geffect)
            pheno = np.append(np.repeat(1, len(cases)), np.repeat(0, len(controls)))
            self.index = np.append(cases, controls)
            return pheno

        else:
            sys.exit('values are missing')

        
    def define_causal(self, causal):
        """Simulates a causal vector

        :causal: Number of causal SNPs
        :returns: Vector of Positions of causal SNPs

        """

        if type(causal) is float:
            causal = np.random.binomial(1, causal, self.p)
            return causal

        if type(causal) is np.ndarray:
            if (any(x not in set(causal) for x in [0,1])) or (len(set(causal))>2) or len(causal)!=self.p:
                sys.exit('invalid causal vector')
            else:
                return causal
        

    def liability_model(self, num_cases, num_controls, threshold, hera, geffect, max_iter=10000):
        """simulates cases and controls

        :num_cases: number of cases
        :num_controls: number of controls
        :threshold: liability threshold
        :returns: subject IDs of cases and controls

        """
        container_cases = []
        container_controls = []

        for item in range(max_iter):
            pheno = np.sqrt(hera)*geffect + np.sqrt(1-hera)*np.random.normal(0, np.sqrt(1), self.n)
            if len(container_cases) < num_cases:
                container_cases.append(np.argwhere(pheno >=threshold))

            if len(container_controls) < num_controls:
                container_controls.append(np.argwhere(pheno <threshold))

            if len(container_cases)>=num_cases and len(container_controls) >=num_controls:
                break;

        # unlist
        container_cases = [item for sublist in container_cases for item in sublist]
        container_controls = [item for sublist in container_controls for item in sublist]

        # remove overhanging samples
        container_cases = container_cases[0:(num_cases)]
        container_controls = container_controls[0:(num_controls)]

        return container_cases, container_controls

    def gwas(self, phenotypes, genotypematrix=None, num_threads=1):
        """Computes summary statistics 

        :phenotypes: Vector of phenotypes
        :genotypematrix: Matrix of genotypes
        :returns: pandas DataFrame with the summary statistics

        """

        if genotypematrix is None:
            genotypematrix = self.genotypematrix[self.index,:]
            output = pymp.shared.array((self.p, 3))
            print('running GWAS on %i individuals and %i SNPs' % genotypematrix.shape)
        else:
            output = pymp.shared.array((genotypematrix.shape[1], 3))
            print('running GWAS on %i individuals and %i SNPs' % genotypematrix.shape)

        if len(np.unique(phenotypes))==2:
            # logistic regression
            with pymp.Parallel(num_threads) as th:
                for p in th.range(genotypematrix.shape[1]):
                    model = sm.GLM(phenotypes, sm.add_constant(genotypematrix[:,p]), family=sm.families.Binomial()).fit()
                    output[p,:] = model.params[1], model.bse[1], model.pvalues[1]
        else:
            # normal regression
            with pymp.Parallel(num_threads) as th:
                for p in th.range(genotypematrix.shape[1]):
                    model = sm.GLM(phenotypes, sm.add_constant(genotypematrix[:,p]), family=sm.families.Binomial()).fit()
                    output[p,:] = model.params[1], model.bse[1], model.pvalues[1]

        output = pd.DataFrame(output, columns=['beta', 'std_err', 'p_value'], index=self.geno_index)

        return output
