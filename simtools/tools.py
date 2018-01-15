"""
Functions to perform various tasks often needed when using simulated data
"""
import pandas as pd
import re
import numpy as np
import pymp
import scipy
import os
import glob
from matplotlib import pyplot as plt
import statsmodels.api as sm
from subprocess import Popen
from simtools.genotypes import ReadPlink


def qqplot(dat, grouping='pheno', pvalue='pvalue'):
    """Plots QQ-Plot of GWAS summary statistics

    :dat: pandas data input or numpy array of p-values
    :returns: qqplot

    """
    if isinstance(dat, np.ndarray):
        n = dat.shape[0]
        dat = -1*np.log(np.sort(dat))
        null = -1*np.log(np.arange(1, n+1)/n)
        plt.plot(null, dat, '.')
        plt.plot(null, null)

    elif isinstance(dat, pd.DataFrame):
        for name, group in dat.groupby(grouping):
            temp = group[['pvalue']]
            n = temp.shape[0]
            temp = temp.sort_values(pvalue)
            temp['null'] = np.arange(1, n+1)/ n
            temp = temp.apply(lambda x: -1*np.log(x))
            plt.plot(temp.null, temp.pvalue,'.', label=str(name))
            plt.plot(temp.null, temp.null)
        plt.legend(loc='upper left')

    plt.xlabel('Expected')
    plt.ylabel('Observed')
    plt.title('QQ plot')
    plt.show()

def gwas(phenotypes, genotypematrix, num_threads=1, verbose=False):
    """Computes summary statistics 

    :phenotypes: Vector of phenotypes
    :genotypematrix: numpy matrix of genotypes
    :num_threads: number of threads to use
    :verbose: print status of the GWAS
    :returns: pandas DataFrame with the summary statistics

    """

    output = pymp.shared.array((genotypematrix.shape[1], 3))
    if verbose:
        print('running GWAS on %i individuals and %i SNPs' % genotypematrix.shape)

    if len(np.unique(phenotypes))==2:
        # logistic regression
        with pymp.Parallel(num_threads) as th:
            for p in th.range(genotypematrix.shape[1]):
                model = sm.GLM(phenotypes, sm.add_constant(genotypematrix[:,p]),
                        family=sm.families.Binomial()).fit()
                output[p,:] = model.params[1], model.bse[1], model.pvalues[1]
    else:
        # normal regression
        with pymp.Parallel(num_threads) as th:
            for p in th.range(genotypematrix.shape[1]):
                model = sm.GLM(phenotypes, sm.add_constant(genotypematrix[:,p]),
                        family=sm.families.Gaussian()).fit()
                output[p,:] = model.params[1], model.bse[1], model.pvalues[1]

    output = pd.DataFrame(output, columns=['beta', 'std_err', 'p_value'],
            index=range(genotypematrix.shape[1]))

    return output


def inflation_factor(s, stats_type='pval', rounding=3):
    """Computes the inflation factor lambda for a given set of p-values

    :s: array of z-statistics, p-values or chisquare statistics
    :stats_type: type of statistic (values can be: pval, chisq, or z
    :rounding: number of decimal places lambda will be given
    :returns: genomic inflation factor lambda

    """
    if stats_type=='pval':
        z = scipy.stats.norm.ppf(s/2)
    if stats_type=='chisq':
        z = np.sqrt(s)
    if stats_type=='z':
        z = s

    lamb = np.round(np.median(z**2)/0.454, rounding)
    return lamb


def causal_models(causal_id, num_pheno):
    """Generates adjacency matrices for specific causal models

    :causal_id: model number
    :num_pheno: number of phenotypes
    :returns: adjacency matrix of size num_pheno*num_pheno

    """
    lamb = np.zeros((num_pheno, num_pheno))
    if causal_id == 1:
        num_edges = num_pheno -1
        arr = np.arange(num_pheno)
        np.random.shuffle(arr)
        i = 0
        while i <num_edges:
            lamb[arr[i+1], arr[i]] = 1
            i += 1
        return lamb
    if causal_id == 2:
        arr = np.arange(num_pheno)
        np.random.shuffle(arr)
        root = arr[0]
        to_edghe = np.delete(arr, 0)
        for i in to_edghe:
            lamb[i, root] = 1
        return lamb
    if causal_id == 3:
        num_pairs = 1
        arr = np.arange(num_pheno)
        np.random.shuffle(arr)
        if num_pheno % 2 !=0:
            arr = np.delete(arr,0)
        edges = np.split(arr, len(arr)/2) 
        for i in edges:
            lamb[i[0], i[1]] = 1
        return lamb
    if causal_id == 4:
        arr = np.arange(num_pheno)
        np.random.shuffle(arr)
        arr = np.append(arr, arr[0])
        num_edges = num_pheno
        i = 0
        while i <num_edges:
            lamb[arr[i+1], arr[i]] = 1
            i += 1
        return lamb

def causal_model_adjustment(lamb):
    """Evaluates an adjacency matrix and adjust it if necessary
    to ensure that evaluation remains fair.

    For example, if A -> B -> C then this would be the same as A -> C

    :lamb: adjacency matrix
    :returns: adjusted adjacency matrix

    """
    # identify if causal re-direction is required
    rowsum = np.sum(lamb, axis=1)
    colsum = np.sum(lamb, axis=0)
    if ((np.any(rowsum > 1) & np.any(colsum >1)) | np.sum(colsum == 0) > 1):
        return lamb
    else:
        start = np.where(rowsum==0)[0]
        end = np.where(colsum==0)[0]
        lamb[end, start] = 1.0
        return lamb 

class Plink(object):

    """Docstring for Plink. """

    def __init__(self, plink_stem, plink_path='auto'):
        """

        :plink_stem: plink stem file
        :plink_path: plink location

        """
        self._plink_stem = plink_stem

        if plink_path=='auto':
            self._bin_plink = '/usr/bin/plink'
        else:
            self._bin_plink=plink_path

        self._plink = ReadPlink(self._plink_stem)

    def gwas_plink(self, phenotype, subjects=None):
        """Calls plink to compute the summary statistics

        :phenotype: phenotype
        :returns: DataFrame with the results

        """
        model_paramters = 'hide-covar'

        # write input files for plink
        pheno = pd.DataFrame()

        if  phenotype.ndim > 1:
            pheno = pd.DataFrame(phenotypes,
                    columns=['p'+str(k+1) for k in range(phenotype.shape[1])])
        else:
            pheno = pd.DataFrame(phenotype, columns=['pheno'])
        if subjects is None:
            subjects = self._plink.fam.iid.values

        fam = pd.DataFrame({'fid':subjects, 'iid':subjects})
        pheno = pd.concat([fam, pheno], axis=1)
        pheno.to_csv('/tmp/temp_pheno', index=False, sep=' ')

        if os.path.isfile('/tmp/plink_temp.log'):
            for filename in glob.glob('/tmp/plink_temp*'):
                os.remove(filename)

        output_location = '/tmp/plink_temp'
        with open(os.devnull, 'w') as fp:
            plink_run = Popen([self._bin_plink,'--bfile', self._plink_stem,
                    "--allow-no-sex",
                    '--pheno', '/tmp/temp_pheno',
                    '--all-pheno', 
                    "--linear", model_paramters,
                    '--out', output_location], stdout=fp)
            plink_run.wait()

        assoc_files = glob.glob('/tmp/plink_temp.*.assoc.*')
        results = pd.DataFrame()
        if len(assoc_files) == 0:
            ValueError('Plink failed to generate any resutls')
        elif len(assoc_files) > 1:
            out_results = []
            for f in assoc_files:
                pheno_name = re.sub(r'(/tmp/plink_temp\.|\.assoc\.|\.linear|\.logistic)',
                        '',f)
                print(pheno_name)
                temp = pd.read_table(f, delim_whitespace=True)
                temp['Phenotype'] = pheno_name
                out_results.append(temp)
            results = pd.columns(out_results, axis=0)
        else:
            results = pd.read_table(assoc_files[0], delim_whitespace=True)

        if results.shape[0] < 1:
            ValueError('Plink failed to generate any resutls or result were not imported')

        return results

    def prs(self, weights, subjects=None):
        """Computes the prs score

        :weights: pandas dataframe with rsid, allele, weights
        :subjects: optional subjects
        :returns: pandas dataframe with results

        """
        if not isinstance(weights, pd.DataFrame):
            ValueError('weights must be a pandas dataframe with rsid, allele and weight')

        output_location = '/tmp/plink_prs'

        subjectfile = '/tmp/plink_filter_subjects.fam'
        if subjects is None:
            self._plink.fam.to_csv(subjectfile, index=False, sep=' ')
        else:
            pd.DataFrame({'iid': subjects,
                'fid': subjects}).to_csv(subjectfile, index=False, sep=' ')

        # write score file
        scorefile = '/tmp/plink_prs_weights.score'
        weights.to_csv(scorefile, index=False, sep=' ', header=False)

        with open(os.devnull, 'w') as fp:
            plink_run = Popen([self._bin_plink, '--bfile', self._plink_stem,
                '--allow-no-sex', '--score', scorefile,
                '--keep', subjectfile, '--out', output_location], stdout=fp)

            plink_run.wait()

        results = pd.read_table(output_location+'.profile', delim_whitespace=True)
        return results

    def clumping(self, summary_stats, p1=0.0001, p2=0.01, r2=0.5, kp=250):
        """Performs clumping with given summary stats

        :summary_stats: TODO
        :p1: TODO
        :p2: TODO
        :r2: TODO
        :kp: TODO
        :returns: TODO

        """
        output_location = '/tmp/plink_clump'
        var_names = summary_stats.columns
        summary_stats.to_csv(output_location+'.report', index=False, sep=' ')


        with open(os.devnull, 'w') as fp:
            plink_run = Popen([self._bin_plink, '--bfile', self._plink_stem,
                '--allow-no-sex', '--clump', output_location+'.report',
                '--clump-p1', p1, '--clump-p2', p2,
                '--clump-r2', r2, 
                '--clump-kp', kp, 
                '--out', output_location], stdout=fp)

            plink_run.wait()

        results = pd.read_table(output_location+'.clumped', delim_whitespace=True)
        return results

    def pruning(self, kp, step, r2):
        """Performs variant pruning and returns a set of variants in linkage 
        equilibrium

        :kp: TODO
        :step: TODO
        :r2: TODO
        :returns: TODO

        """

        output_location = '/tmp/plink_prune'

        with open(os.devnull, 'w') as fp:
            plink_run = Popen([self._bin_plink, '--bfile', self._plink_stem,
                '--allow-no-sex',
                '--indep-pairwise', kp+' '+step+' '+r2,
                '--out', output_location], stdout=fp)

            plink_run.wait()

        results = pd.read_table(output_location'.prune.in', delim_whitespace=True)
        return results
