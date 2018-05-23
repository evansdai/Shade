import pandas as pd
import hdf
import re
import numpy as np
from funcs.fig.distribution import distribution
from scipy.stats import chisquare
import matplotlib.pyplot as plt
all_length = pd.read_csv(r'C:\Users\evans\Dropbox\Shade\database\full_length.csv', index_col=0)
sigs = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\tsne\tsne_2_0.1.csv", index_col=0)
expr = hdf.read(r"C:\Users\evans\Dropbox\Shade\raw\jiaoshenmehaoen.h5")[0]
expr = expr[-expr['Name'].str.contains(';')]
sigs['motif21'] = all_length.loc[all_length.index.isin(
    sigs.index), 'Motif21_3']
sigs['Name'] = expr.loc[sigs.index, 'Name']

expr['motif21'] = all_length.loc[all_length.index.isin(
    expr.index), 'Motif21_3']
motifs = hdf.read(r"C:\Users\evans\Dropbox\Shade\motifx\phoasphat_motifs.h5")[1]
kinase_substrate = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\database\arabidopsis_kinase_substrate.txt", sep='\t', index_col='Target')
kinase_substrate = kinase_substrate[-kinase_substrate['Kinase'].isin([np.nan])]
name_cheat = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\database\go.csv", index_col=0)
sigs['Annotations'] = [j if isinstance(j, str) else j[0] for j in [
    name_cheat.loc[i, 'Annotations'] for i in sigs['Name']]]


def phosphatase_motif_to_re(s):
    s = s.replace(' ', '')
    s = re.sub('-', '', s, count=0, flags=0)
    # s=s.replace('X','\w')
    s = re.sub('X', '\w', s, count=0, flags=0)
    return s


def my_chisquare(obs1, obs2):
    r = (np.sum(obs1) + 0.0) / np.sum(obs2)
    obs2 = np.array(obs2) * r
    return chisquare(obs1, obs2)


def list_upper(l):
    return [i.upper() for i in l]


list_upper(sigs.loc[sigs['label'] == 1].index)

motif = '[MILV]\\w'


def motif_chi(motif):
    found_1 = [re.findall(motif, i, flags=0) != []
               for i in sigs.loc[sigs['label'] == 1]['motif21']].count(True)

    notfound_1 = [re.findall(motif, i, flags=0) != []
                  for i in sigs.loc[sigs['label'] == 1]['motif21']].count(False)

    found_2 = [re.findall(motif, i, flags=0) != []
               for i in sigs.loc[sigs['label'] == 2]['motif21']].count(True)

    notfound_2 = [re.findall(motif, i, flags=0) != []
                  for i in sigs.loc[sigs['label'] == 2]['motif21']].count(False)

    [found_1, notfound_1, found_2, notfound_2]
    return my_chisquare([found_1, notfound_1], [found_2, notfound_2])[-1], found_1, notfound_1, found_2, notfound_2


def batch_motif_chi():

    res = []
    motif_list = map(phosphatase_motif_to_re, motifs.index)
    motifs.index = motif_list
    for i in motifs.index:
        kinase = motifs.loc[i]['kinase']
        p, found1, notfound1, found2, notfound2 = motif_chi(i)
        res.append({'motif': i, 'kinase': kinase, 'p': p, 'found1': found1,
                    'notfound1': notfound1, 'found2': found2, 'notfound2': notfound2})

    return res


def single_motif_peptides(motif):
    # print len([i for i in sigs.loc[sigs['label'] == 2].index if len(re.findall(motif, sigs.loc[i]['motif21'])) > 0])
    # print len([i for i in sigs.loc[sigs['label'] == 1].index if len(re.findall(motif, sigs.loc[i]['motif21'])) > 0])

    temp = pd.Series([i for i in sigs['motif21']
                      if len(re.findall(motif, i)) > 0])
    return sigs.loc[sigs['motif21'].isin(temp)]


def count_na(series):
    return series.tolist()


def kinase_family(sigs_or_subset):
    # test1=sigs.loc[sigs['label']==1]
    # test1=pd.DataFrame(kinase_family(test1))
    # test1['family'].value_counts()
    # test2=sigs.loc[sigs['label']==2]
    # test2=pd.DataFrame(kinase_family(test2))
    # test2['family'].value_counts()
    # test_all=pd.DataFrame(kinase_family(expr))
    # test_all['family'].value_counts()
    res = []
    namelist = list(set(sigs_or_subset['Name']))
    for transcript in namelist:
        # try:
            # print kinase_substrate.loc[transcript]
        # except:
            # pass
        if transcript in kinase_substrate.index and isinstance(kinase_substrate.loc[transcript], pd.DataFrame):
            first_kinase = kinase_substrate.loc[transcript]['Kinase'].value_counts(
            ).index[0]
            first_family = kinase_substrate.loc[kinase_substrate['Kinase']
                                                == first_kinase]['Family'].value_counts().index[0]
            res.append({'Name': transcript, 'kinase': first_kinase,
                        'family': first_family})
        elif transcript in kinase_substrate.index and isinstance(kinase_substrate.loc[transcript], pd.Series):
            first_kinase = kinase_substrate.loc[transcript].name
            first_family = kinase_substrate.loc[transcript]['Family']
            res.append({'Name': transcript, 'kinase': first_kinase,
                        'family': first_family})

        else:
            res.append(
                {'Name': transcript, 'kinase': 'None', 'family': 'None'})
    return res


def my_chisquare(obs1, obs2):
    r = (np.sum(obs1) + 0.0) / np.sum(obs2)
    obs2 = np.array(obs2) * r
    return chisquare(obs1, obs2)


if __name__ == '__main__':
    test = batch_motif_chi()
    test = pd.DataFrame(test).set_index('motif')
    test
    test[test['p'] < 0.05]

    test
    # test.to_csv('fsjfl.csv')
    # test2=test.dropna()
    # distribution(test2['p'])
    # plt.show()
    # single_motif_peptides('[LIV]\w[KRH]\w\w[st]\w\w\w[LIV]')
    # single_motif_peptides('Gs')
    # test=single_motif_peptides('[RK]\w\ws')
    # test=single_motif_peptides('sP[RK]')

    single_motif_peptides('s\w[DE]').to_csv('fdsfa.csv')

    # kinase_substrate=pd.read_csv(r"C:\Users\evans\Dropbox\Shade\database\arabidopsis_kinase_substrate.txt",sep='\t',index_col='Target')
    # test=kinase_substrate.loc[kinase_substrate.index.isin(sigs['Name'])]
    # test1=kinase_substrate.loc[kinase_substrate.index.isin(sigs.loc[sigs['label']==1,'Name'])]
    # test2=kinase_substrate.loc[kinase_substrate.index.isin(sigs.loc[sigs['label']==2,'Name'])]
    # test.columns
    # test1['Family'].value_counts()
    # test1['Pathway'].value_counts()
    # test2['Pathway'].value_counts()
    # test2['Family'].value_counts()
