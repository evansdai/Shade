from scipy.stats import chisquare
import pandas as pd
import hdf
cheat = hdf.read(r"C:\Users\evans\Dropbox\Shade\raw\jiaoshenmehaoen.h5")[0]
up_down = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\tsne\tsne_7.csv", index_col=0)
up_down['Name2'] = [cheat.loc[up_down.index[i]]['Name2']
                    for i in range(len(up_down))]
in_nuc = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\GO\in_nucleo.csv", index_col=0)

print up_down['Name2'].isin(in_nuc.index).sum()
print up_down.loc[up_down['label'].isin([7, 2, 4]), 'Name2'].isin(in_nuc.index).value_counts()
print up_down.loc[up_down['label'].isin([6, 3, 1]), 'Name2'].isin(in_nuc.index).value_counts()

chisquare([64, 29], [102, 56])

chisquare([47, 27], [102, 56])
102
11 / 7.0
chisquare([11, 7], [102, 56])
chisquare()
chisquare([102, 56], [798, 1311 - 798])
102 / 56.0
798 / (1311 - 798.0)
1311 - 798
