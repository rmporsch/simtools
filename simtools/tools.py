import pandas as pd
import numpy as np
from matplotlib import pyplot as plt


def ggplot(self, dat, grouping='pheno', pvalue='pvalue'):
    """Plots QQ-Plot of GWAS summary statistics

    :dat: pandas data input or numpy array of p-values
    :returns: ggplot

    """
    if isinstance(dat, np.ndarray):
        n = dat.shape[0]
        dat = -1*np.log(np.sort(dat))
        null = -1*np.log(np.arange(1, n+1)/n)
        plt.plot(null, dat, '.')
        plt.plot(null, null)

    else if isinstance(dat, pd.DataFrame):
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
