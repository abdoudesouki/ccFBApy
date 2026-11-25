# Utilities

Python helper functions that have no R help page.

## The R data sets (no rdata needed)

The data sets of sybilccFBA (`R/sybilccFBA/data/*.RData`) are shipped with the package as SBML and CSV
files, read with [`ccfbapy.datasets`](../data.md#data-sets-in-the-package) (`load_model`, `load_dataset`). The
package does not read `.RData` files and does not need the package `rdata`.

The files were exported from the `.RData` files with the development tool `tools/make_package_data.py` of the
repository. It uses `tools/sybil2cobra.py` (`read_RData`, `sybil2cobra`: sybil `modelorg` → `cobra.Model`)
and needs `pip install rdata`. The tool reads each export again and compares it with the R data set (ids, bounds,
stoichiometry, GPR rules, objective, FBA optimum; tables: exact values). See
[Data sets](../data.md#how-the-data-sets-were-made).

## scatterPlot

*`ccfbapy/genFuns.py`*

```python
from ccfbapy import scatterPlot
scatterPlot(x, y, xlbl='x', ylbl='y', fs=(6, 6), labels=None)
```

Scatter plot of measured (`x`) against predicted (`y`) values, e.g. growth rates. It draws the
line `y = x` (red, dashed) and a linear fit, and prints the Pearson correlation `R` and `R^2`.
`labels` is an optional list with one label per point (e.g. carbon source names). `fs` is the
figure size. Needs matplotlib.

```python
from ccfbapy import scatterPlot
scatterPlot([0.66, 0.47, 0.29], [0.47, 0.52, 0.32], xlbl='Measured growth rate (1/h)',
            ylbl='Predicted growth rate (1/h)', labels=['glc', 'glyc', 'ac'])
```

## remove_enzyme_layer

*`ccfbapy/cfba_moment.py`*

```python
from ccfbapy import remove_enzyme_layer
remove_enzyme_layer(model)
```

Removes the variables and constraints that an earlier call of [`cfba_moment`](cfba_moment.md) or
[`ccFBA_mr`](cfba_moment_mr.md) added to the solver problem of `model`, so they can be added
again. `cfba_moment` and `ccFBA_mr` call it on their copy of the model. You need it only when you
reuse the `model` they return (`res['model']`).

- variables: `g_<gene>`, `gr_<gene>_<rxn>`, `cpx_*`;
- constraints: `<rxn>`, `<rxn>_b`, `<rxn>_kcat`, `<rxn>_b_kcat`, `<rxn>_<gene>`, `Aux_*`, `mfe_<gene>`, `MW`.

Reaction variables and metabolite constraints are never removed.
