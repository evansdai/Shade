import hdf
import pandas as pd
from funcs.fig.scatter import plot
import dup
import numpy as np
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
sigs = pd.read_csv(r'C:\Users\evans\Dropbox\Shade\database\sigs.csv', index_col=0)
trans = pd.read_csv(r'C:\Users\evans\Dropbox\Shade\database\transcript_FPKM.csv', index_col=0)
phos = hdf.read(r"C:\Users\evans\Dropbox\Shade\raw\jiaoshenmehaoen.h5")[0]
phos = phos[-phos['Name'].str.contains(';')]

common_gene = [i for i in phos['Name2'] if i in trans.index]
sig_gene = [i for i in phos['Name2'] if i in sigs.index]
trans = trans.loc[common_gene]
trans.columns
phos = dup.dup_index_to_unique(phos.set_index('Name2'), met='mean')
trans = dup.dup_index_to_unique(trans, met='sum')
phos = phos[['L1', 'L2', 'S5_1', 'S5_2', 'S15_1', 'S15_2', 'S30_1', 'S30_2']]

trans.columns
test = pd.concat([phos['L1'], trans['col_0WL-1']], axis=1)

trans.columns
phos_l = np.mean(phos[['L1', 'L2']], axis=1)
phos_s = np.mean(phos[['S30_1', 'S30_2']], axis=1)
trans_s = np.mean(trans[['col_0SH-1', 'col_0SH-2', 'col_0SH-3']], axis=1)
trans_l = np.mean(trans[['col_0WL-1', 'col_0WL-2', 'col_0WL-3']], axis=1)

test1 = pd.concat([phos_l, phos_s], axis=1)
test2 = pd.concat([trans_l, phos_l], axis=1)
test3 = pd.concat([trans_s - trans_l, phos_s - phos_l], axis=1)
test4 = pd.concat([trans_s / trans_l, phos_s / phos_l], axis=1)
import delete
test4 = delete.del_row_with_zero_and_inf(test4)
test4.columns = ['phosphorylation_FC', 'transcription_FC']
p = plot(test4, drop_zero=True, log=True, spearman=True,
         savepath=None, summary=None, pointplot=False, postpros=None)
p.savefig('fkjdslf.png', dpi=900, width=6, height=4)
test = pd.concat([test4, trans_l, trans_s], axis=1)
test.columns = ['phos_fc', 'trans_fc', 'trans_l', 'trans_s']
formula = 'phos_fc~trans_l+trans_s+trans_fc'
anova_results = anova_lm(ols(formula, data=test).fit())
print ols(formula, data=test).fit().summary()
print anova_results
