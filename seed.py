import pandas as pd
import hdf
sigs = hdf.read(r"C:\Users\evans\Dropbox\Shade\raw\sigs.h5")[0]
sigs['Name2'][0]
seed = pd.read_csv(r'C:\Users\evans\Dropbox\Shade\database\seed.txt', sep='\t', index_col=0)
seed.columns
seed_sig = pd.DataFrame()
for i in sigs.index:
    temp = seed.loc[seed['genotype'].str.contains(sigs.loc[i, 'Name2'])]
    if not temp.empty:
        temp['anno'] = sigs.loc[i, 'Annotations']
        temp['gene'] = sigs.loc[i, 'Name2']
        seed_sig = pd.concat([seed_sig, temp])

seed_sig.to_csv(r'C:\Users\evans\Dropbox\Shade\database\seed_sig.csv')
