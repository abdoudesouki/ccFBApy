from scipy.stats import pearsonr
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import pearsonr

def scatterPlot(x,y,xlbl='x',ylbl='y',fs=(6, 6)):
    corr, _ = pearsonr(x, y)
    print("Pearson corr R=",round(corr,3),"R^2=",round(corr**2,3))

    Mx=max(max(y),max(y))*1.1
    f, ax = plt.subplots(figsize=fs)
    ax.set(xlim=(0, Mx), ylim=(0, Mx))
    diag_line, = ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c="red")
    ax.scatter(x, y, c = 'blue')
    #Plot line
    ax.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)))
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    plt.show()