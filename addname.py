
import pandas as pd
import hdf

if __name__ == '__main__':
    name_cheat
    test = pd.read_csv(r'C:\Users\evans\Dropbox\Shade\raw\allgene.csv', index_col='Accession')
    name_cheat = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\database\go.csv", index_col=0)
    temp = [name_cheat.loc[name_cheat['Gene'] == i, 'Annotations']
            for i in test['Name2']]

    temp2 = []
    for i in temp:
        if isinstance(i, str):
            temp2.append(i)

        elif isinstance(i, pd.Series):
            try:
                temp2.append(i[0])
            except:
                temp2.append('unknown')
        else:
            raise RuntimeError(i)
    test['Annotations'] = temp2
    test.to_csv(r"C:\Users\evans\Dropbox\Shade\network\background.csv")
