import pandas as pd
import hdf
from drawlines import batch_draw_timecourse
sigs = hdf.read(r"C:\Users\evans\Dropbox\Shade\raw\sigs.h5")[0]


go_up = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\GO\tsne0.1_bp.txt", sep='\t', index_col='Term')
go_down = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\GO\tsne0.1_bp_down.txt", sep='\t', index_col='Term')


def string_to_list(s):
    return s.split(', ')


def list_of_string_to_list(l):
    res = []
    for i in l:
        res.extend(string_to_list(i))
    return res


def go_figure4(goterm):
    aba_up = list_of_string_to_list(
        go_up.loc[go_up.index.str.contains(goterm), 'Genes'].values)
    aba_down = list_of_string_to_list(
        go_down.loc[go_down.index.str.contains(goterm), 'Genes'].values)
    aba_down.extend(aba_up)
    aba_genes = list(set(aba_down))
    test = sigs.loc[sigs['Name2'].isin(aba_genes)]
    batch_draw_timecourse(test, outputroot=r'C:\Users\evans\Dropbox\Shade\figures\FIGURE4\{}\\'.format(goterm))
    return test


if __name__ == '__main__':
    go_figure4('ethylene')
