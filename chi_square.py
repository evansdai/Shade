from scipy.stats import chisquare
import pandas as pd
import hdf
import numpy as np

# def my_chisquare(obs1,obs2):
#
#     obs1=(np.array(obs1)+0.0)/np.sum(obs1)
#     obs2=(np.array(obs2)+0.0)/np.sum(obs2)
#     return chisquare(obs1,obs2)


cheat = hdf.read(r"C:\Users\evans\Dropbox\Shade\raw\jiaoshenmehaoen.h5")[0]
up_down = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\tsne\tsne_5.csv", index_col=0)
up_down['Name2'] = [cheat.loc[up_down.index[i]]['Name2']
                    for i in range(len(up_down))]
in_nuc = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\GO\in_nucleo.csv", index_col=0)

print up_down['Name2'].isin(in_nuc.index).sum()
print up_down.loc[up_down['label'].isin([1, 4]), 'Name2'].isin(in_nuc.index).value_counts()
print up_down.loc[up_down['label'].isin([2, 5]), 'Name2'].isin(in_nuc.index).value_counts()
print up_down.loc[up_down['label'].isin([3]), 'Name2'].isin(in_nuc.index).value_counts()
chisquare([49, 31], [42, 17])
chisquare([49, 31], [20, 8])
chisquare([42, 17], [20, 8])
chisquare([20, 8], [102, 56])

chisquare([42, 17], [20, 8])
102
11 / 7.0
chisquare([10000, 7], [102, 56])
chisquare()
chisquare([102, 56], [798, 1311 - 798])
102 / 56.0
798 / (1311 - 798.0)
1311 - 798
