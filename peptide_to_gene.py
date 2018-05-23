import pandas as pd
import hdf

allgene = pd.read_csv(r'C:\Users\evans\Dropbox\Shade\raw\allgene.csv', index_col='Accession')
loc = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\raw\Profile_all2076_B_withlocation.csv", index_col=0)
test = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\heatmap\ja.csv")
loc.columns
test
temp = allgene.loc[allgene['Name2'].isin(test['Name2'].tolist())]
temp
loc.loc[temp.index, 'Phospho'] + loc.loc[temp.index, 'Location']
