import funcs.hdf as hdf
import pandas as pd
from funcs.fig.scatter import plot
import matplotlib.pyplot as plt
import seaborn as sns


def timecource_feature_rearrange(df, ts=[0, 0, 5, 5, 15, 15, 30, 30]):
    o = pd.DataFrame()
    for i in range(len(df)):
        o = o.append(pd.DataFrame(
            {'time': ts, 'expr': df.iloc[i].tolist(), 'index': [df.iloc[i].name] * len(ts)}))
    o = o[['time', 'expr', 'index']]
    return o


class draw_timecourse(object):
    """docstring for draw_timecourse."""

    def __init__(self, which):
        super(draw_timecourse, self).__init__()
        if isinstance(which, str):
            self.which = [which]
        elif isinstance(which, list):
            self.which = which
        else:
            raise TypeError
        self.h5file = r"C:\Users\evans\Dropbox\Shade\raw\timecourse.h5"

        def query_timecourse():
            tc = hdf.read(self.h5file)[0]
            return tc.loc[tc['index'].isin(self.which)]
        self.timecourse = query_timecourse()

    def draw(self, color=None, legend_loc=None, rm_legend=False):
        g = plot(self.timecourse, color=color, pointplot=True)
        if legend_loc is not None:
            g.axes.legend(loc='center right', bbox_to_anchor=(1.6, 0.5))
        if rm_legend == True:
            g.axes.legend().remove()
        self.fig = g

# h = hdf.read(r"C:\Users\evans\Dropbox\Shade\raw\jiaoshenmehaoen.h5")[0]
# test = h[h.columns[-8:]]


# def query_timecourse(which=['lSSLsLNLSNQPAAIAAR']):
#     tc=hdf.read(r"C:\Users\evans\Dropbox\Shade\raw\timecourse.h5")[0]
#     return tc.loc[tc['index'].isin(which)]


if __name__ == '__main__':

    dtc = draw_timecourse(
        which=['lSSLsLNLSNQPAAIAAR', 'sFADIFQADMGHPVVQQPsRPASR'])

    dtc.draw(color='#5833c2', legend_loc='1', rm_legend=False)
