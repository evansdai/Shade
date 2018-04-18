import hdf
import pandas as pd
from funcs.fig.distribution import distribution
import matplotlib.pyplot as plt
import re
from statsmodels.stats.anova import anova_lm
from statsmodels.formula.api import ols
from funcs.fig.color import colors10


db = hdf.read(r"C:\Users\evans\Dropbox\Shade\database\phosphat_20160120.h5")[0]
expr = hdf.read(r"C:\Users\evans\Dropbox\Shade\raw\jiaoshenmehaoen.h5")[0]
expr = expr[-expr['Name'].str.contains(';')]
full_length = pd.read_csv(r'C:\Users\evans\Dropbox\Shade\database\full_length.csv', index_col=0)
sty = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\figures\FIGURE1\sty_counts.csv", index_col=0)
exist_peptides = pd.read_csv(r'C:\Users\evans\Dropbox\Shade\database\peptide_existed.csv', index_col=0)
exist_transcripts = db['code'].tolist()
transcripts = expr['Name'].value_counts().index

exist_peptides.loc[sty.loc[sty['y'] > 0].index]['existed'].value_counts()
full_length
