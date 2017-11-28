#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pymp
import statsmodels.api as sm
import numpy as np
import pandas as pd
from itertools import product


# generate some data
p = 10000
n = 1000
print('Generationg Data')
y = np.random.normal(0, np.sqrt(1), n)
dat = pd.DataFrame(columns=['V'+str(x) for x in range(p)], index=range(n))

for item in range(p):
    dat[item] = np.random.normal(0, np.sqrt(1), n)

def regression(y, x):
    """computes a simple linear regression

    :y: TODO
    :x: TODO
    :returns: TODO

    """
    model = sm.regression.linear_model.OLS(y, sm.add_constant(x))
    results = model.fit()
    return results.params[1], results.bse[1], results.pvalues[1]


def info(title):
    print(title)
    print('module name:', __name__)
    print('parent process:', os.getppid())
    print('process id:', os.getpid())

#print(regression(y, dat[1]))

print('Computing Summary Stats')
output_array = pymp.shared.array((p, 3), dtype='float')

with pymp.Parallel(5) as thread:
    thread.print(thread.num_threads, thread.thread_num)
    for index in thread.range(p):
        output_array[index,:] = regression(y, dat[index])


output = pd.DataFrame(output_array, columns=['slope', 'se', 'pvalues'])
print(output)


