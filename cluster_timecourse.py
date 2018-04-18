import pandas as pd
import funcs.hdf as hdf
import scipy.cluster.hierarchy as hy
import sklearn.mixture.gaussian_mixture as gm
import sklearn.manifold.t_sne as tsne
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import normalize
from sklearn.cluster import KMeans
from funcs.fig.scatter import plot
import matplotlib.pyplot as plt
from drawlines import draw_timecourse
import numpy as np
import seaborn as sns
import statsmodels.api as sm


class cluster_timecourse(object):
    """docstring for cluster_timecourse."""

    def __init__(self, num_c):
        super(cluster_timecourse, self).__init__()
        self.num_c = num_c
        # all peptides
        # self.df = hdf.read(r"C:\Users\evans\Dropbox\Shade\raw\jiaoshenmehaoen.h5")[0]
        # self.df = self.df[self.df.columns[-8:]].replace(np.nan, 0.0)
        # self.X = self.df[self.df.columns[-8:]].as_matrix()

        # significant peptides
        self.df = hdf.read(r"C:\Users\evans\Dropbox\Shade\maSigPro\sig_by_perm.h5")[0].replace(np.nan, 0.0)
        # TODO a abnormal value : 'iTESEQAAGDsDEGVDSITTGGR'
        # L1       1.272
        # L2       0.736
        # S5_1     1.016
        # S5_2     8.807
        # S15_1    1.259
        # S15_2    1.296
        # S30_1    1.443
        # S30_2    1.108
        # Name: iTESEQAAGDsDEGVDSITTGGR, dtype: float64
        self.df = self.df.drop('iTESEQAAGDsDEGVDSITTGGR')
        self.X = self.df.as_matrix()

        def add_polynormial(x, dim=1):
            res = np.matrix(x).T
            for i in range(dim + 1):
                if i > 1:
                    res = np.concatenate([res, np.matrix(x**(i)).T], axis=1)
            res = sm.add_constant(res)
            return res

        def my_polyfit(x, y, dim):
            X = add_polynormial(x, dim=dim)
            est = sm.OLS(y, X)
            est2 = est.fit()
            log_p = np.log(est2.pvalues)
            llf = est2.llf
            par = est2.params
            return np.append(llf, np.concatenate([log_p, par]))

        # quadratic features
        self.X2 = np.concatenate((self.X, self.X**2), axis=1)
        self.fit = np.matrix([my_polyfit(x=np.array([0, 0, 5, 5, 15, 15, 30, 30]),
                                         y=self.X[i, ], dim=1) for i in range(self.X.shape[0])])

        self.fit2 = np.matrix([my_polyfit(x=np.array([0, 0, 5, 5, 15, 15, 30, 30]),
                                          y=self.X[i, ], dim=2) for i in range(self.X.shape[0])])

        self.fit3 = np.matrix([my_polyfit(x=np.array([0, 0, 5, 5, 15, 15, 30, 30]),
                                          y=self.X[i, ], dim=3) for i in range(self.X.shape[0])])

        self.ratio = self.cal_ratio_of_change_of_mean(self.df).as_matrix()

        self.rank = self.sort_mean(self.df)

        self.all = np.concatenate((self.rank, self.fit), axis=1)

        def normalize(array):
            res = (array - np.min(array)) / (np.max(array) - np.min(array))
            res[res == np.inf] = 0
            res[res == -np.inf] = 0
            return np.nan_to_num(res, copy=True)

        self.all = np.concatenate([normalize(self.all[:, i])
                                   for i in range(self.all.shape[1])], axis=1)

    def cal_ratio_of_change_of_mean(self, df):
        def divide_filter(s1, s2):
            '''
            input two dataframe Series
            output one dataframe Series
            '''
            o = np.log2(s1 / s2)
            o.name = str(s1.name) + '/' + str(s2.name)
            return o
        ratio5 = divide_filter(df[['S5_1', 'S5_2']].mean(
            axis=1), (df[['L1', 'L2']].mean(axis=1)))
        ratio15 = divide_filter(df[['S15_1', 'S15_2']].mean(
            axis=1), (df[['L1', 'L2']].mean(axis=1)))
        ratio30 = divide_filter(df[['S30_1', 'S30_2']].mean(
            axis=1), (df[['L1', 'L2']].mean(axis=1)))
        res = pd.DataFrame(
            {'ratio5': ratio5, 'ratio15': ratio15, 'ratio30': ratio30})
        res = res.replace(np.inf, 0.0)
        res = res.replace(-np.inf, 0.0)
        res = res.replace(np.nan, 0.0)
        return res

    def sort_mean(self, df):
        l = df[['L1', 'L2']].mean(axis=1)
        s5 = df[['S5_1', 'S5_2']].mean(axis=1)
        s15 = df[['S15_1', 'S15_2']].mean(axis=1)
        s30 = df[['S30_1', 'S30_2']].mean(axis=1)
        res = []
        for i in range(len(l)):
            temp = {0: l.iloc[i], 5: s5.iloc[i],
                    15: s15.iloc[i], 30: s30.iloc[i]}
            res.append([i[0] for i in sorted(temp.items(), cmp=None,
                                             key=lambda (k, v): v, reverse=False)])
        return np.matrix(res)

    def cluster(self):
        # g = gm.GaussianMixture(n_components=self.num_c, covariance_type='full', tol=1e-3, reg_covar=1e-6, max_iter=100, n_init=1, init_params='kmeans',
                               # weights_init=None, means_init=None, precisions_init=None, random_state=None, warm_start=False, verbose=0, verbose_interval=10)

        # g = DBSCAN(eps=0.05, min_samples=2, metric='euclidean',
                   # metric_params=None, algorithm='auto', leaf_size=30, p=None, n_jobs=1)
        g = KMeans(n_clusters=self.num_c, init='k-means++', n_init=10, max_iter=300, tol=1e-4,
                   precompute_distances='auto', verbose=0, random_state=None, copy_x=True, n_jobs=1, algorithm='auto')
        # self.model = g.fit(self.all)
        # self.label = self.model.predict(self.all).tolist()
        self.label = g.fit(self.all).labels_

    def tsne(self, legend=True):
        X_embedded = tsne.TSNE(n_components=2).fit_transform(self.all, y=None)
        dX_embedded = pd.DataFrame(X_embedded)
        dX_embedded.columns = ['Component1', 'Component2']
        dX_embedded.index = self.df.index
        dX_embedded['label'] = self.label
        dX_embedded['label'] = dX_embedded['label'] + 1
        self.dX_embedded = dX_embedded

        self.fig = plot(self.dX_embedded, drop_zero=False, log=False,
                        spearman=False, savepath=None, summary=None, legend=legend, fit_reg=False)

    def draw_all_lines(self, savefile=None):
        cluster = self.dX_embedded
        # TODO color palette
        # colors = sns.color_palette('husl', self.num_c).as_hex()
        colors = sns.color_palette().as_hex()
        for i in range(self.num_c):
            try:
                dtc = draw_timecourse(
                    which=cluster.loc[cluster['label'] == i + 1].index.tolist())
                dtc.draw(color=colors[i], legend_loc=None, rm_legend=True)
                if savefile is not None:
                    dtc.fig.get_figure().savefig(savefile.format(i))
                del dtc
            except:
                print i


if __name__ == '__main__':

    n_cluster = 2
    ctc = cluster_timecourse(n_cluster)
    ctc.cluster()
    ctc.tsne(legend=True)
    ctc.fig.savefig(r"C:\Users\evans\Dropbox\Shade\tsne\tsne.png", dpi=900)
    ctc.draw_all_lines(savefile=None)
    ctc.dX_embedded.to_csv(r"C:\Users\evans\Dropbox\Shade\tsne\tsne_5.csv")
