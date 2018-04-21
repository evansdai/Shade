from scipy.stats import chisquare
import pandas as pd
import hdf
import numpy as np
from funcs.fig.distribution import distribution
import matplotlib.pyplot as plt
from scipy import stats
import glob
import re


def my_chisquare(obs1, obs2):
    r = (np.sum(obs1) + 0.0) / np.sum(obs2)
    obs2 = np.array(obs2) * r
    return chisquare(obs1, obs2)


cheat = hdf.read(r"C:\Users\evans\Dropbox\Shade\raw\jiaoshenmehaoen.h5")[0]
up_down = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\tsne\tsne_2_0.1.csv", index_col=0)
up_down['Name2'] = [cheat.loc[up_down.index[i]]['Name2']
                    for i in range(len(up_down))]
up_down['Name'] = [cheat.loc[up_down.index[i]]['Name']
                   for i in range(len(up_down))]
# up_down.to_csv(r"C:\Users\evans\Dropbox\Shade\tsne\tsne2_0.1_name.csv")
expr = hdf.read(r"C:\Users\evans\Dropbox\Shade\raw\jiaoshenmehaoen.h5")[0]
expr = expr[-expr['Name'].str.contains(';')]
up_down.to_csv(r'C:\Users\evans\Dropbox\Shade\figures\FIGURE2\sigs.csv')
name_cheat = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\database\go.csv", sep='\t', index_col=0)
useful = [re.findall('\w+.bmp', i)[0][0:-4] for i in glob.glob(r'C:\Users\evans\Dropbox\Shade\figures\profiles\*\*')]
# cheat.loc[up_down2.loc[up_down2.index.isin(useful)].index]['Annotations']
# set([i for i in name_cheat.loc[up_down2['Name']]['Annotations'] if isinstance(i,str) and i[0:2]!='AT'])

in_nuc = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\GO\in_nucleo.csv", index_col=0)
in_loss = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\database\loss_of_function_vegetative.tsv", sep='\t', index_col=0)
in_chr = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\database\chr_weight_po.csv", index_col=0)
inter = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\database\interaction.csv", sep='\t', index_col=0)
print up_down['Name'].isin(inter.index).sum()

print up_down.loc[up_down['label'].isin([1]), 'Name2'].isin(in_nuc.index).value_counts()

[i in in_nuc.index for i in set(
    up_down.loc[up_down['label'].isin([1]), 'Name2'].tolist())].count(True)
[i in in_nuc.index for i in set(
    up_down.loc[up_down['label'].isin([1]), 'Name2'].tolist())].count(False)
[i in in_nuc.index for i in set(
    up_down.loc[up_down['label'].isin([2]), 'Name2'].tolist())].count(True)
[i in in_nuc.index for i in set(
    up_down.loc[up_down['label'].isin([2]), 'Name2'].tolist())].count(False)
[i in in_nuc.index for i in set(cheat['Name2'].tolist())].count(False)
[i in in_nuc.index for i in set(cheat['Name2'].tolist())].count(True)

print up_down.loc[up_down['label'].isin([2]), 'Name2'].isin(in_nuc.index).value_counts()

cheat['Name2'].isin(in_nuc.index).value_counts()
print up_down.loc[up_down['label'].isin([1]), 'Name2'].isin(in_loss.index).value_counts()
print up_down.loc[up_down['label'].isin([2]), 'Name2'].isin(in_loss.index).value_counts()

my_chisquare([104, 45], [798, 518])
my_chisquare([96, 53], [798, 518])
all = list(set(up_down['Name2']))
up_first = {}
up_second = {}
down_first = {}
down_second = {}
up = {}
down = {}
backgroup = {}
inter['second'].tolist().count(all[0])
for i in all:
    if i in up_down.loc[up_down['label'] == 1, 'Name2'].tolist():
        up_first.update({i: inter['first'].tolist().count(i)})
        up_second.update({i: inter['second'].tolist().count(i)})
        up.update({i: inter['first'].tolist().count(
            i) + inter['second'].tolist().count(i)})

    elif i in up_down.loc[up_down['label'] == 2, 'Name2'].tolist():
        down_first.update({i: inter['first'].tolist().count(i)})
        down_second.update({i: inter['second'].tolist().count(i)})
        down.update({i: inter['first'].tolist().count(
            i) + inter['second'].tolist().count(i)})
    else:
        raise RuntimeError(i)

for i in list(set(expr['Name2'])):
    backgroup.update({i: inter['first'].tolist().count(
        i) + inter['second'].tolist().count(i)})

backgroup = pd.Series(backgroup)
distribution(backgroup, namelist=None, title='',
             xlabel='', ylabel='KDE', kde=True)
plt.show()
up_first = pd.Series(up_first)
up_second = pd.Series(up_second)
down_first = pd.Series(down_first)
down_second = pd.Series(down_second)
down = pd.Series(down)
up = pd.Series(up)
p = distribution([down, up], namelist=['Down', 'Up'])
p.set(ylim=[0, 0.06])
p.figure

p.figure.savefig(r"C:\Users\evans\Dropbox\Shade\figures\FIGURE2\distribution_of_interaction_node_0.1.png", dpi=900)
plt.show()
stats.ks_2samp(down, up)

up_down.loc[up_down['label'].isin([1]), 'Name']
print up_down.loc[up_down['label'].isin([2]), 'Name'].isin(inter.index).value_counts()
my_chisquare([41, 27], [60, 39])


my_chisquare([60, 39], [798, 513])
