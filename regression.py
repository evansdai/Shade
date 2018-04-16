# use regression to determine

import pandas as pd
import numpy as np
from sklearn import datasets, linear_model
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
from scipy import stats
import hdf
#
# diabetes = datasets.load_diabetes()
# X = diabetes.data
# y = diabetes.target

df = hdf.read(r"C:\Users\evans\Dropbox\Shade\raw\jiaoshenmehaoen.h5")[0]
y = df[df.columns[-8:]].replace(np.nan, 0.0).as_matrix()[1, ]
X = np.array([0, 0, 5, 5, 15, 15, 30, 30])


def add_polynormial(x, dim=1):
    res = np.matrix(x).T
    for i in range(dim + 1):
        if i > 1:
            res = np.concatenate([res, np.matrix(x**(i)).T], axis=1)
    res = sm.add_constant(res)
    return res


def my_polyfit(X, y, dim):
    X2 = add_polynormial(X, dim=dim)
    est = sm.OLS(y, X2)
    est2 = est.fit()
    log_p = np.log(est2.pvalues)
    llf = est2.llf
    par = est2.params
    return np.append(llf, np.concatenate([log_p, par]))
    # return [est2.fvalue,


if __name__ == '__main__':
    my_polyfit(X, y, dim=2)

    np.array([1, 2, 3, np.nan])
    np.replace()

    test = {5: 1, 'b': 2, 'c': 3}
    sorted(test)
