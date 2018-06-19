import numpy as np
from pyplink import PyPlink
import os


class ReadPlink(object):
    """docstring for ReadPlink."""

    def __init__(self, plinkstem):
        """Plink init functions."""
        # super(ReadPlink).__init__()
        self._plinkstem = plinkstem
        self._bim_path = os.path.basename(self._plinkstem)+'.bim'
        self._bed_path = os.path.basename(self._plinkstem)+'.bed'
        self._fam_path = os.path.basename(self._plinkstem)+'.fam'

        self.plinkfile = PyPlink(self._plinkstem)
        self.fam = self.plinkfile.get_fam()
        self.bim = self.plinkfile.get_bim()
        self.N = self.fam.shape[0]
        self.P = self.bim.shape[0]
        self.subject = self.fam['iid'].values
        self.variants = self.bim.index.values

    def get_gentoypematrix(self, marker=None, subjects=None):
        """Read bed file.

        :param marker: list of SNPs
        :param subjects: list of subjects
        :returns: genotype-matrix of size subjects*marker

        """
        if marker is None:
            p_size = self.P
            marker = self.variants
        else:
            p_size = len(marker)

        if subjects is None:
            n_size = self.N
            subjects = self.fam.index.values
        else:
            n_size = len(subjects)
            subjects = self.fam[self.fam.iid.isin(subjects)].index.values

        genotypematrix = np.zeros((n_size, p_size), dtype=np.int8)

        j = 0
        for m, g in self.plinkfile.iter_geno_marker(marker):
            genotypematrix[:, j] = g[subjects]
            j += 1

        genotypematrix[genotypematrix < 0] = 0

        return genotypematrix
