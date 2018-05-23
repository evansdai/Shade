import pandas as pd
import hdf
import subprocess
from delete import del_row_with_zero

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
    temp = test.fc('30', fc=1.5)
    len(temp)
    temp.to_csv('fdsf30.csv')
