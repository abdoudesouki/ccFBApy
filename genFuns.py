from scipy.stats import pearsonr
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import pearsonr

def scatterPlot(x,y,xlbl='x',ylbl='y',fs=(6, 6),labels=None):
    # labels: optional list of point labels (e.g. carbon source names), same length as x and y
    corr, _ = pearsonr(x, y)
    print("Pearson corr R=",round(corr,3),"R^2=",round(corr**2,3))

    Mx=max(max(x),max(y))*1.1  # common limit for both axes so no point is cut off
    f, ax = plt.subplots(figsize=fs)
    ax.set(xlim=(0, Mx), ylim=(0, Mx))
    diag_line, = ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c="red")
    ax.scatter(x, y, c = 'blue')
    if labels is not None:
        for xi, yi, lbl in zip(x, y, labels):
            ax.annotate(lbl, (xi, yi), textcoords="offset points", xytext=(3, 3), fontsize=7)
    #Plot line
    ax.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)))
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    plt.show()