import hdf
import pandas as pd

import numpy as np
sigs = pd.read_csv(r'C:\Users\evans\Dropbox\Shade\database\sigs.csv', index_col=0)
# expr=hdf.read(r"C:\Users\evans\Dropbox\Shade\raw\jiaoshenmehaoen.h5")[0][['L1','L2','S5_1','S5_2','S15_1','S15_2','S30_1','S30_2']]

expr = hdf.read(r"C:\Users\evans\Dropbox\Shade\raw\jiaoshenmehaoen.h5")['prot_RMP_nopoint_dupmean'][['L1', 'L2', 'S5_1', 'S5_2', 'S15_1', 'S15_2', 'S30_1', 'S30_2']]


def cal_regulate(df):
    def divide_filter(s1, s2):
        '''
        input two dataframe Series
        output one dataframe Series
        '''
        # o = np.log(s1 / s2)
        o = s1 / s2
        o.name = str(s1.name) + '/' + str(s2.name)
        return o
    ratio5 = divide_filter(df[['S5_1', 'S5_2']].mean(
        axis=1), (df[['L1', 'L2']].mean(axis=1)))
    ratio15 = divide_filter(df[['S15_1', 'S15_2']].mean(
        axis=1), (df[['L1', 'L2']].mean(axis=1)))
    ratio30 = divide_filter(df[['S30_1', 'S30_2']].mean(
        axis=1), (df[['L1', 'L2']].mean(axis=1)))
    return pd.DataFrame({'ratio5': ratio5, 'ratio15': ratio15, 'ratio30': ratio30}).dropna()


cal_regulate(expr.loc[sigs[sigs['label'].isin([1, 4])].index])

cal_regulate(expr.loc[set(sigs[sigs.columns[0]].tolist())]
             ).to_csv('ratios.csv')


cal_downregulate(expr.loc[sigs[sigs['label'].isin([2, 5])].index])

sm.anova.anova_lm()
