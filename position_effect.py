import hdf
import pandas as pd
from funcs.fig.distribution import distribution
import matplotlib.pyplot as plt
import re
from statsmodels.stats.anova import anova_lm
from statsmodels.formula.api import ols
sty = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\figures\FIGURE1\sty_counts.csv", index_col=0)
test3 = pd.read_csv(r"C:\Users\evans\Dropbox\Shade\database\relative_pos.csv", index_col=0)
test4 = pd.concat([test3, sty], axis=1, join='inner')

p = distribution([test4.loc[test4['s'] > 0, 'relative_pos'], test4.loc[test4['t'] > 0, 'relative_pos'],
                  test4.loc[test4['y'] > 0, 'relative_pos']], namelist=['s', 't', 'y'], kde=False)
p.set(xlim=[0, 1])
p.figure
plt.show()
p.figure.savefig('fdjslfd.png', dpi=900)
# res=ols('relative_pos~y',data=test4).fit()
# anova_lm(res,type=1)
