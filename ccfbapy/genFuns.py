# General helper functions of ccFBApy.
import numpy as np
from scipy.stats import pearsonr


def scatterPlot(x, y, xlbl='x', ylbl='y', fs=(6, 6), labels=None):
    """Scatter plot of measured (x) against predicted (y) values (e.g. growth rates): the line y = x (red,
    dashed), a linear fit, and the Pearson correlation R and R^2 (printed).
    labels: optional list with one label per point (e.g. carbon source names). fs: figure size.
    Needs matplotlib (pip install ccFBApy[plot])."""
    import matplotlib.pyplot as plt
    corr, _ = pearsonr(x, y)
    print("Pearson corr R=", round(corr, 3), "R^2=", round(corr ** 2, 3))

    Mx = max(max(x), max(y)) * 1.1  # common limit for both axes so no point is cut off
    f, ax = plt.subplots(figsize=fs)
    ax.set(xlim=(0, Mx), ylim=(0, Mx))
    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c="red")
    ax.scatter(x, y, c='blue')
    if labels is not None:
        for xi, yi, lbl in zip(x, y, labels):
            ax.annotate(lbl, (xi, yi), textcoords="offset points", xytext=(3, 3), fontsize=7)
    # linear fit
    ax.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)))
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    plt.show()
