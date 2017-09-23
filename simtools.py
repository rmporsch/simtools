import pandas as pd
import numpy as np
import sys
from sklearn.preprocessing import scale

class Simtools(object):

    """Containing multiple functions to simulate phenotypes"""

    def __init__(self, matrix):
        """TODO: to be defined1.
        
        :matrix TODO
        """
        self.genotypematrix = matrix
        self.n = matrix.shape[0]
        self.p = matrix.shape[1]
        self.index = np.arange(0, self.n - 1)
            

    def simple_phenotype(self, causal, hera, liability=None):
        """simulates a phenotypes (continues or binary)

        :hera: TODO
        :liability: TODO
        :returns: TODO

        """

        self.causal = self.define_causal(causal)
        geffect = scale(np.dot(self.genotypematrix, self.causal))
        
        if liability is None:
            pheno = hera*geffect + np.sqrt(1-hera)*np.random.normal(0, np.sqrt(1), self.n)
            return pheno

        elif liability is tuple and len(liability)==3:
            threshold = liability[0]
            ncases = liability[1]
            ncontrols = liability[2]

            cases, controls = self.liability_model(ncases, ncontrols, threshold, hera, geffect)
            pheno = np.array(np.repeat(1, len(cases)), np.repeat(0, len(controls))).flatten()
            self.index = np.array(cases, controls).flatten()
            return pheno

        else:
            sys.exit('values are missing')

        
    def define_causal(self, causal):
        """Simulates a causal vector

        :causal: TODO
        :returns: TODO

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

        :num_cases: TODO
        :num_controls: TODO
        :threshold: TODO
        :returns: TODO

        """
        container_cases = []
        container_controls = []

        for item in range(max_iter):
            pheno = hera*geffect + np.sqrt(1-hera)*np.random.normal(0, np.sqrt(1), self.n)
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
        container_cases = container_cases[0:(num_cases -1)]
        container_controls = container_controls[0:(num_controls -1)]

        return container_cases, container_controls

