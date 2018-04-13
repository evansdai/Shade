# rhdf5 in R cannot save the rownames when saving dataframe as h5
# this function is aimed to find name for those dataframe that lose index

import hdf
import numpy as np


def find_name():
    origin = hdf.read(r"C:\Users\evans\Dropbox\Shade\raw\jiaoshenmehaoen.h5")[0]
    lost_index = hdf.read(r"C:\Users\evans\Dropbox\Shade\maSigPro\sig_by_perm_quadratic.h5")[0]
    found_name = []
    origin = origin.replace(np.nan, 0.0)
    lost_index = lost_index.replace(np.nan, 0.0)
    for j in range(len(lost_index)):
        matched = {}
        for i in range(len(origin)):
            matched.update({origin.index[i]: (
                abs(lost_index.iloc[j] - origin[origin.columns[-8::]].iloc[i]) < 0.0001).sum().sum()})
        res = [k for k, v in matched.iteritems() if v == 8]
        if len(res) == 1:
            found_name.append(res[0])
        else:
            print res
            print j
            print lost_index.iloc[j]
            raise RuntimeError
    lost_index.index = found_name
    return lost_index


if __name__ == '__main__':
    test = find_name()
    test.to_hdf(r"C:\Users\evans\Dropbox\Shade\maSigPro\sig_by_perm_quadratic1.h5", key='one')
    hdf.read(r"C:\Users\evans\Dropbox\Shade\maSigPro\sig_by_perm.h5")[0].to_csv(r"C:\Users\evans\Dropbox\Shade\maSigPro\sig_by_perm.csv")
