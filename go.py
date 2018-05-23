import pandas as pd
import hdf
import seaborn as sns
from drawlines import batch_draw_timecourse
from scipy import stats
import os
sns.set(font_scale=2.5, style="white", color_codes=True)
# sigs = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\network\sigs.csv", index_col=0)
sigs = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\network\sigs0.3.csv", index_col=0)
# TODO not significant
# sigs = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\network\background.csv", index_col=0)

# name_cheat = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\database\go.csv", sep='\t', index_col=0)
name_cheat = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\database\identifier_mappings.txt", sep='\t')
# tes.loc[tes['Source']=='Gene Name']
go_up = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\GO\tsne0.1_bp.txt", sep='\t', index_col='Term')
go_down = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\GO\tsne0.1_bp_down.txt", sep='\t', index_col='Term')

# TODO gobp , cc, mf
# go_all = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\GO\tsne0.1_bp_control.txt", sep='\t', index_col='Term')
# go_all = pd.read_csv(r'C:\Users\evans\Dropbox\Shade\GO\tsne0.1_mf_control.txt', sep='\t', index_col='Term')

# go_all
go_ALL = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\database\go.csv")
go_ALL['Term'].value_counts()
# go_all.loc[go_all.index.str.contains('brass')]


def string_to_list(s):
    return s.split(', ')


def list_of_string_to_list(l):
    res = []
    for i in l:
        res.extend(string_to_list(i))
    return res


def go_figure4(goterm):
    # TODO add go_all
    # aba_up = list_of_string_to_list(
    #     go_up.loc[go_up.index.str.contains(goterm), 'Genes'].values)
    # aba_down = list_of_string_to_list(
    #     go_down.loc[go_down.index.str.contains(goterm), 'Genes'].values)
    # aba_down.extend(aba_up)
    df_go_term = go_all.loc[go_all.index.str.contains(goterm), 'Genes']
    go_genes = list_of_string_to_list(
        df_go_term.values)
    go_genes = list(set(go_genes))
    if len(go_genes) == 0:
        pass
    else:
        try:
            os.mkdir(r'C:\Users\evans\Dropbox\Shade\figures\FIGURE4\{}'.format(goterm))
        except:
            pass
        # print go_genes
        # print sigs['Gene']
        test = sigs.loc[sigs['Gene'].isin(go_genes)]
        # print test
        batch_draw_timecourse(test, df_pos=True, outputroot=r'C:\Users\evans\Dropbox\Shade\figures\FIGURE4\{}\\'.format(goterm))
        s = ''
        for i in set(test['Gene'].values):
            print i
            s = s + i + '\n'
        print s

        with open(r'C:\Users\evans\Dropbox\Shade\figures\FIGURE4\{}\\{}.txt'.format(goterm, goterm), mode='w') as file:
            file.write(s)
        df_go_term.to_csv(r'C:\Users\evans\Dropbox\Shade\figures\FIGURE4\{}\\{}_go.csv'.format(goterm, goterm))
        test.to_csv(r'C:\Users\evans\Dropbox\Shade\figures\FIGURE4\{}\\{}.csv'.format(goterm, goterm))

    return 0

# TODO use go_ALL instead of go_all
# go_ALL.loc[go_ALL['Term'].str.contains('splicing')].columns


def new_go_figure4(goterm):
    df_go_term = go_ALL.loc[go_ALL['Term'].str.contains(goterm)]
    go_genes = list(set(df_go_term['Gene']))
    if len(go_genes) == 0:
        pass
    else:
        try:
            os.mkdir(r'C:\Users\evans\Dropbox\Shade\figures\FIGURE4\{}'.format(goterm))
        except:
            pass
        test = sigs.loc[sigs['Gene'].isin(go_genes)]

        batch_draw_timecourse(test, df_pos=True, outputroot=r'C:\Users\evans\Dropbox\Shade\figures\FIGURE4\{}\\'.format(goterm))
        s = ''
        for i in set(test['Gene'].values):
            print i
            s = s + i + '\n'
        print s

        with open(r'C:\Users\evans\Dropbox\Shade\figures\FIGURE4\{}\\{}.txt'.format(goterm, goterm), mode='w') as file:
            file.write(s)
        df_go_term.to_csv(r'C:\Users\evans\Dropbox\Shade\figures\FIGURE4\{}\\{}_go.csv'.format(goterm, goterm))
        test.to_csv(r'C:\Users\evans\Dropbox\Shade\figures\FIGURE4\{}\\{}.csv'.format(goterm, goterm))

    return 0


def go_label(goterm):
    aba_up = list_of_string_to_list(
        go_up.loc[go_up.index.str.contains(goterm), 'Genes'].values)
    aba_down = list_of_string_to_list(
        go_down.loc[go_down.index.str.contains(goterm), 'Genes'].values)
    aba_down.extend(aba_up)
    aba_genes = list(set(aba_down))
    test = sigs.loc[sigs['Name2'].isin(aba_genes)]
    sigs[goterm] = [1 if i in test.index else 0 for i in sigs.index]
    return sigs


class count_sig_in_goterm(object):
    """docstring for count_sig_in_goterm."""

    def __init__(self):
        super(count_sig_in_goterm, self).__init__()
        self.go = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\database\go.csv", index_col=0)
        self.sig = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\network\sigs.csv", index_col=0)
        self.phos = pd.read_csv(r'C:\Users\evans\Dropbox\Shade\raw\allgene.csv', index_col=0)
        self.sig_list = list(set(self.sig['Protein'].tolist()))
        # TODO one gene could be count in both up list and down list, which is fine
        self.sig_up_list = list(
            set(self.sig.loc[self.sig['label'] == 2, 'Protein'].tolist()))
        self.sig_down_list = list(
            set(self.sig.loc[self.sig['label'] == 1, 'Protein'].tolist()))
        self.phos_list = list(set(self.phos['Name'].tolist()))

    def count(self, keyword=None):
        if not keyword is None:
            go_key = self.go.loc[self.go['Term'].str.contains(keyword)]
            go_key_sig = go_key.loc[go_key.index.isin(self.sig_list)]
            go_key_all = go_key.loc[go_key.index.isin(self.phos_list)]
            return go_key_sig
        else:
            print len(self.go.loc[self.sig['Protein']])
            print len(self.go.loc[self.phos['Name']])

    def chi_one(self, goterm):
        go_one = self.go.loc[self.go['Term'] == goterm]
        in_sig = go_one.loc[go_one.index.isin(self.sig_list)]
        in_sig_up = go_one.loc[go_one.index.isin(self.sig_up_list)]
        in_sig_down = go_one.loc[go_one.index.isin(self.sig_down_list)]
        in_phos = go_one.loc[go_one.index.isin(self.phos_list)]
        not_in_phos = go_one.loc[~go_one.index.isin(self.phos_list)]
        all_in_phos = self.go.loc[self.go.index.isin(self.phos_list)]
        all_not_in_phos = self.go.loc[~self.go.index.isin(self.phos_list)]

        print len(in_phos)
        print len(not_in_phos)
        print len(all_in_phos)
        print len(all_not_in_phos)
        print stats.chisquare([len(in_sig_up), len(in_sig_down)], [(len(in_sig_up) + len(in_sig_down)) / 2.0, (len(in_sig_up) + len(in_sig_down)) / 2.0])
        print stats.chisquare([len(in_phos), len(not_in_phos)], [(len(in_phos) + len(not_in_phos)) / 2.0, (len(in_phos) + len(not_in_phos)) / 2.0])


if __name__ == '__main__':
    pass
    # TODO given a list of gene, draw profiles
    # csig = count_sig_in_goterm()
    # test=csig.count('GTP')
    # test
    # gtp=sigs.loc[sigs['Protein'].isin(test.index.tolist())]
    # gtp=gtp.replace('UK/UPRT1','UK_UPRT1')
    # batch_draw_timecourse(gtp, df_pos=True, outputroot=r"C:\Users\evans\Dropbox\Shade\Endomembrane system\sigs0.1//")

    # csig.chi_one('GTPase')
    # sns.set(font_scale=2.5, style="white", color_codes=True)
    # goterm = 'cadmium'
    # with open(r'C:\Users\evans\Dropbox\Shade\figures\FIGURE4\{}\\{}.txt'.format(goterm, goterm), mode='w') as file:
    #     file.write('dfsd')
    # for i in ['actin', 'abscisic', 'light', 'cadmium', 'chloroplast', 'chromophore', 'cold', 'cytoskeleton', 'jasmonic', 'phospho', 'photo', 'processing', 'splicing', 'sugar', 'translational', 'vesicle', 'viral']:
    #     # sigs = go_label(i)
    #     go_figure4(i)
    #
    # sigs.to_csv(r"C:\Users\evans\Dropbox\Shade\network\sigs.csv")

    # TODO new bp
    # for i in ['nucleus organization','RNA secondary structure unwinding','cytokinin','autophosphorylation','cytoskeleton','defense response to bacterium','defense response to fungus, incompatible interaction','jasmonic','brassinosteroid','chloroplast avoidance movement','phototropism','trichome morphogenesis','innate immune response','endocytosis','chloroplast accumulation movement','circadian rhythm','negative regulation of long-day photoperiodism, flowering']:
    #     new_go_figure4(i)

    # TODO new mf

    # for i in ['translation initiation factor activity', 'protein binding', 'mRNA binding', 'actin filament binding', 'GTPase', 'protein phosphatase inhibitor activity']:
    #     new_go_figure4(i)
    new_go_figure4('GTP')
    # for i in ['microtubule associated complex', 'microtubule', 'intracellular', 'plasma membrane', 'cytoskeleton']:
    #     new_go_figure4(i)

    # change annotations for sigs csv
    # name_cheat.head()
    # sigs = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\network\sigs0.3.csv")
    # sigs.head()
    # res = []
    # for i in sigs['Gene']:
    #     print i
    #     name = name_cheat.loc[(name_cheat['Preferred_Name'] == i) & (
    #         name_cheat['Source'] == 'Gene Name'), 'Name']
    #     if len(name) == 0:
    #         res.append(i)
    #     elif len(name) == 1:
    #         res.append(name.values[0])
    #     else:
    #         raise RuntimeError(i)
    # res = pd.Series(res)
    # res[res != sigs['Annotations']]
    # sigs['Annotations'] = res
    # sigs.to_csv(r"C:\Users\evans\Dropbox\Shade\network\sigs0.3_newannotations.csv")
    # sigs[res != sigs['Annotations']]
