import pandas as pd
import funcs.hdf as hdf
from summary import change_index
from dup import dup_index_to_unique
import numpy as np
import scipy.stats as stats
from funcs.fig.scatter import plot
import matplotlib.pyplot as plt


def remove_multihit_peptides(df):
    return df.loc[~df.index.str.contains(';')]  # remove multihit peptides


def ati_merge_subtype(df):
    df.index = [i.split('.')[0] for i in df.index]
    return dup_index_to_unique(df, met='mean')

# TODO 编写函数计算比值


def cal_fold_change(df):
    ratio5 = (df[['S5_1', 'S5_2']].mean(axis=1)) / \
        (df[['L1', 'L2']].mean(axis=1))
    ratio15 = (df[['S15_1', 'S15_2']].mean(axis=1)) / \
        (df[['L1', 'L2']].mean(axis=1))
    ratio30 = (df[['S30_1', 'S30_2']].mean(axis=1)) / \
        (df[['L1', 'L2']].mean(axis=1))
    return pd.DataFrame({'ratio5': ratio5, 'ratio15': ratio15, 'ratio30': ratio30}).dropna()


def cal_flipped_change(df):
    def flip_divide(s1, s2):
        '''
        input two dataframe Series
        output one dataframe Series
        '''
        o = s1 / s2
        o[o < 1] = 1 / o
        o.name = str(s1.name) + '/' + str(s2.name)
        return o
    ratio5 = flip_divide(df[['S5_1', 'S5_2']].mean(
        axis=1), (df[['L1', 'L2']].mean(axis=1)))
    ratio15 = flip_divide(df[['S15_1', 'S15_2']].mean(
        axis=1), (df[['L1', 'L2']].mean(axis=1)))
    ratio30 = flip_divide(df[['S30_1', 'S30_2']].mean(
        axis=1), (df[['L1', 'L2']].mean(axis=1)))
    return pd.DataFrame({'ratio5': ratio5, 'ratio15': ratio15, 'ratio30': ratio30}).dropna()


def cal_upregulate(df):
    def divide_filter(s1, s2):
        '''
        input two dataframe Series
        output one dataframe Series
        '''
        o = s1 / s2 - 1
        o = o[o > 0]
        o.name = str(s1.name) + '/' + str(s2.name)
        return o
    ratio5 = divide_filter(df[['S5_1', 'S5_2']].mean(
        axis=1), (df[['L1', 'L2']].mean(axis=1)))
    ratio15 = divide_filter(df[['S15_1', 'S15_2']].mean(
        axis=1), (df[['L1', 'L2']].mean(axis=1)))
    ratio30 = divide_filter(df[['S30_1', 'S30_2']].mean(
        axis=1), (df[['L1', 'L2']].mean(axis=1)))
    return pd.DataFrame({'ratio5': ratio5, 'ratio15': ratio15, 'ratio30': ratio30}).dropna()


def cal_downregulate(df):
    def divide_filter(s1, s2):
        '''
        input two dataframe Series
        output one dataframe Series
        '''
        o = 1 - s1 / s2
        o = o[o > 0]
        o.name = str(s1.name) + '/' + str(s2.name)
        return o
    ratio5 = divide_filter(df[['S5_1', 'S5_2']].mean(
        axis=1), (df[['L1', 'L2']].mean(axis=1)))
    ratio15 = divide_filter(df[['S15_1', 'S15_2']].mean(
        axis=1), (df[['L1', 'L2']].mean(axis=1)))
    ratio30 = divide_filter(df[['S30_1', 'S30_2']].mean(
        axis=1), (df[['L1', 'L2']].mean(axis=1)))
    return pd.DataFrame({'ratio5': ratio5, 'ratio15': ratio15, 'ratio30': ratio30}).dropna()


# test2 = (df[['S5_1', 'S5_2']].mean(axis=1)) / (df[['L1', 'L2']].mean(axis=1))
# test2.to_csv('test_of_go_brick.csv')
h = hdf.read(r"C:\Users\evans\Dropbox\Shade\raw\jiaoshenmehaoen.h5")

h.keys

df1 = h['prot_RMP_nopoint_dupmean']  # count in multihit peptides
fc1 = cal_downregulate(df1)
fc1.to_csv('fcminus1_of_protein_downreg.csv')
