import pandas as pd
import hdf
import subprocess
from delete import del_row_with_zero
from go import count_sig_in_goterm
# def copy2clip(txt):
#     if isinstance(txt, str):
#         s=txt
#     elif isinstance(txt, list):
#         s=''
#         for i in txt:
#             s=s+i
#
#     else:
#         raise ValueError
#     cmd='echo '+s.strip()+'|clip'
#     print cmd
#     return subprocess.check_call(cmd, shell=True)


class foldchange(object):
    """docstring for foldchange."""

    def __init__(self):
        super(foldchange, self).__init__()
        self.allgene = pd.read_csv(r'C:\Users\evans\Dropbox\Shade\raw\allgene.csv')
        self.allgene = del_row_with_zero(self.allgene)

    def fc(self, time, fc):
        temp = self.allgene.loc[abs(1 - (self.allgene['S{}_1'.format(time)] + self.allgene['S{}_2'.format(
            time)]) / (self.allgene['L1'.format(time)] + self.allgene['L2'.format(time)])) > fc - 1]
        # copy2clip(temp['Name2'].tolist())
        return temp


if __name__ == '__main__':
    test = foldchange()

    for i in ['5', '15', '30']:
        temp = test.fc(i, fc=1.5)
        # temp.to_csv('fdsf30.csv')

        test2 = count_sig_in_goterm(input=temp)
        # test2.count(keyword=None)
        # test2.chi_one('response to abscisic acid')

        test3 = test2.batch_chi()
        test3
        test3.columns
        test3.dropna().to_csv('search_allgo_{}.csv'.format(i))
