import pandas as pd
import hdf

if __name__ == '__main__':
    sigs = hdf.read(r"C:\Users\evans\Dropbox\Shade\raw\sigs.h5")[0]
    relation = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\database\arabidopsis_kinase_substrate.txt", sep='\t')
    substrate = sigs['Name'][0]

    res = pd.DataFrame()

    substrate
    kinase = list(set(relation['Kinase'].tolist()))
    for i in kinase:
        onekinase = relation.loc[relation['Kinase'] == i]
        onepair = onekinase.loc[onekinase['Target'] == substrate]
        if len(onepair) > 0:
            print onepair

    relation_sigs = relation.loc[relation['Target'].isin(
        sigs['Name'].tolist())]

    relation_sigs['Family'].value_counts()
    test['Kinase'].value_counts()
    test['Family'].value_counts()
