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

X2 = sm.add_constant(X)

est = sm.OLS(y, X2)

est2 = est.fit()
print(est2.summary())

est2.fvalue
