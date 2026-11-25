# addGlcTrns: glucose transport constraint

*R help page: `man/addGlcTrns.Rd`. Python: `ccfbapy/addGlcTrns.py`, functions `addGlcTrns()` and
`glc_trns_capacity()`.*

## Description

Adds a glucose transport constraint to the problem: an upper bound on glucose consumption, based on the
kinetics of the *S. cerevisiae* hexose transporters.

The R function is not complete. The code that adds the transporter columns is commented out, so
`trnsCol` is `NULL` and the new row cannot be built. The Python version adds the columns, as that
commented-out code intended.

## Usage

```python
from ccfbapy import addGlcTrns, glc_trns_capacity

addGlcTrns(prob, mod2=None, glc_rxn='R_GLCt1', Si=0.2, transporters=None, trns_ub=1000, trnsCol=None)
glc_trns_capacity(Si=0.2, transporters=None)
```

R:

```r
addGlcTrns(prob, mod2)
```

## Arguments

| Argument | Description |
|---|---|
| `prob` | optlang problem, changed in place: e.g. `res['prob']` of [`cfba_moment`](cfba_moment.md) / [`ccFBA_mr`](cfba_moment_mr.md), `solve_ccFBA_mat(...)['prob']`, or `model.solver` of a cobra model (FBA) |
| `mod2` | accepted for compatibility with R, not needed: the flux variable is found by name |
| `glc_rxn` | id of the glucose transport reaction, default `'R_GLCt1'` (iMM904 from the R data set; `'GLCt1'` in `models/iMM904.xml`). The variable `glc_rxn`, `'R_' + glc_rxn` (columns of `getccFBA_mat`) or `glc_rxn + '_f'` is used: the forward (uptake) flux. |
| `Si` | glucose concentration, in the same unit as the Km values. Default 0.2, as in R. |
| `transporters` | transporters to include, default Hxt1–Hxt7 (as in R). Add `'Gal2'` when galactose is present. |
| `trns_ub` | upper bound of the new transporter variables. Default 1000, as in the commented-out R code. |
| `trnsCol` | optional `dict` transporter → name of an existing variable of `prob`, used instead of a new variable |

## Details

For each transporter `t`, the capacity is `kcat_t * Si / (Si + Km_t)`:

| Transporter | Gene | kcat | Km | Capacity at Si = 0.2 |
|---|---|---|---|---|
| Hxt1 | YHR094C | 41 | 107 | 0.0765 |
| Hxt2 | YMR011W | 16.1 | 2.9 | 1.0387 |
| Hxt3 | YDR345C | 18.5 | 29 | 0.1267 |
| Hxt4 | YHR092C | 12 | 6.2 | 0.3750 |
| Hxt5 | YHR096C | 14 | 10 | 0.2745 |
| Hxt6 | YDR343C | 11.4 | 1.5 | 1.3412 |
| Hxt7 | YDR342C | 11.7 | 1.3 | 1.5600 |
| Gal2 | YLR081W | 17.5 | 1.5 | 2.0588 (only with galactose) |

The kcat and Km values are those of the R code. The R code does not state their units. The genes
are the standard names of these transporters. They all appear in the GPR rule of `R_GLCt1`.

The function adds, for each transporter without an entry in `trnsCol`, a variable
`trns_<t>` with `0 <= x_t <= trns_ub`, and the constraint

```
glcTrns:   v(glc_rxn) - sum_t capacity_t * x_t <= 0
```

With the default `trns_ub = 1000`, the limit is 1000 × 4.79 = 4793, so the constraint does not limit
the uptake. With `trns_ub = 1` the glucose uptake is at most 4.79. The R code asks whether the
transporters should also take part in the enzyme budget (`# Add to crowding constraint (same budget),
Molecular weights required?`). That is not done, in R or in Python. To link the transporters to
existing variables, use `trnsCol`.

## Value

`prob`, with the new variables and the constraint `glcTrns` (as in R). `glc_trns_capacity()` returns a
`dict` transporter → capacity.

## See also

[`cfba_moment`](cfba_moment.md), [`getccFBA_mat`](getccFBA_mat.md), data set [`iMM904`](../data.md#imm904).

## Examples

iMM904 (the R data set), FBA with glucose uptake up to 1000, then with the transport constraint
(`trns_ub = 1`):

```python
from ccfbapy import datasets
from ccfbapy import addGlcTrns, glc_trns_capacity

model = datasets.load_model('iMM904')           # R: data(iMM904)
model.reactions.get_by_id('R_EX_glc_e_').lower_bound = -1000
print('capacities:', glc_trns_capacity())

print('FBA:', model.slim_optimize(), 'glucose transport:', model.optimize().fluxes['R_GLCt1'])
m = model.copy()                                 # the constraint is added to the solver of the copy
prob = addGlcTrns(m.solver, trns_ub=1)
sol = m.optimize()
print('with glcTrns:', sol.objective_value, 'glucose transport:', sol.fluxes['R_GLCt1'])
print(prob.constraints['glcTrns'])
```

On a ccFBA problem (yeast molecular weights; all reactions get the same kcat, 25.5 1/s):

```python
import pandas as pd
from ccfbapy import getccFBA_mat, solve_ccFBA_mat

model = datasets.load_model('iMM904')
mw = datasets.load_dataset('yst_mw')[['gene', 'mw']]   # R: data(yst_mw)
kcat = pd.DataFrame({'rxn_id': [], 'val': [], 'dirxn': []})
mat = getccFBA_mat(model, Kcat=kcat, MW=mw, RHS=0.27, medval=3600*25.5, verboseMode=0)
Res = solve_ccFBA_mat(mat, verboseMode=0)
prob = addGlcTrns(Res['prob'], trns_ub=1)       # column 'R_R_GLCt1'
prob.optimize()
print(Res['obj'], '->', prob.objective.value)
```
