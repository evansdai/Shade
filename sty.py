import hdf
import pandas as pd
import numpy as np
expr = hdf.read(r"C:\Users\evans\Dropbox\Shade\raw\jiaoshenmehaoen.h5")[0]
expr = expr[-expr['Name'].str.contains(';')]
expr.index.value_counts()
expr['Name'].value_counts()
expr['Name2'].value_counts()


def count_minor_letters(s):
    num_s = s.count('s')
    num_t = s.count('t')
    num_y = s.count('y')
    return [num_s, num_t, num_y]


def count_num_of_sites_on_peptide():
    sty_counts = pd.DataFrame(map(count_minor_letters, expr.index))
    sty_counts.index = expr.index
    sty_counts.columns = ['s', 't', 'y']
    np.sum(sty_counts['s'])
    np.sum(sty_counts['t'])
    np.sum(sty_counts['y'])
    sty_counts['s'].value_counts()


def pie():
    p = sty.sum().plot(kind='pie', colors=['#f37a5a', '#f3d45a', '#f3b35a'])
    plt.axis('equal')
    p.figure.savefig('dfs.png', dpi=900)
    plt.show()


if __name__ == '__main__':
    pass
