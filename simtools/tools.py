"""
Functions to perform various tasks often needed when using simulated data
"""
import pandas as pd
import numpy as np
import scipy
from matplotlib import pyplot as plt


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

def gwas(phenotypes, genotypematrix, num_threads=1):
    """Computes summary statistics 

    :phenotypes: Vector of phenotypes
    :genotypematrix: numpy matrix of genotypes
    :num_threads: number of threads to use
    :returns: pandas DataFrame with the summary statistics

    """

    output = pymp.shared.array((genotypematrix.shape[1], 3))
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
                        family=sm.families.Binomial()).fit()
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
