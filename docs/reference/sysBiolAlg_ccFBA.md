# sysBiolAlg_ccFBA: problem object of ccFBA

*R help pages: `man/sysBiolAlg_ccFBA-class.Rd` (class `sysBiolAlg_ccFBA`, aliases
`sysBiolAlg_ccFBA`, `ccFBA`) and `man/initialize-methods.Rd` (its `initialize` method). R source:
`R/sybilccFBA/R/sysBiolAlg_ccFBAClass.R`.*

!!! note "Replaced in Python"
    There is no class. [`solve_ccFBA_mat`](getccFBA_mat.md#solve_ccfba_mat) builds an optlang
    problem from the result of [`getccFBA_mat`](getccFBA_mat.md) and solves it. This replaces
    `sysBiolAlg(..., algorithm="ccFBA")` followed by `optimizeProb()`.

## Description (R)

The class `sysBiolAlg_ccFBA` holds an object of class `optObj`, made to meet the requirements of
the ccFBA algorithm. It extends class `sysBiolAlg` directly, and no methods are defined with it in
the signature.

Objects are created by calls of the form

```r
sysBiolAlg(model, algorithm = "ccFBA", ...)
```

### Slots

| Slot | Content |
|---|---|
| `problem` | object of class `optObj` containing the problem object |
| `algorithm` | name of the algorithm |
| `nr`, `nc` | number of rows and columns of the problem |
| `fldind` | columns (variables) of the fluxes: variable `fldind[i]` is reaction `i` of the network |
| `alg_par` | named list of algorithm specific parameters |

## initialize method

```r
## S4 method for signature 'sysBiolAlg_ccFBA'
initialize(.Object, model, LHS, rlb, rub, rtype, lb, ub, obj, lpdir,
           cnames = NULL, rnames = NULL, pname = NULL, scaling = NULL,
           writeProbToFileName = NULL, ...)
```

| Argument | Description | Python (`getccFBA_mat` result `mat`) |
|---|---|---|
| `.Object` | the problem object | |
| `model` | `modelorg` (the irreversible model). The bounds of the first `react_num(model)` columns are replaced by its bounds. | `mod2` argument of `solve_ccFBA_mat` (default `mat['model']`) |
| `LHS` | constraint matrix: S matrix and ccFBA constraints | `mat['LHS']` |
| `rlb`, `rub` | bounds of the rows (constraints) | `mat['rlb']`, `mat['rub']` |
| `rtype` | row type (`E`, `U`, `L`, ...) | `mat['rtype']` |
| `lb`, `ub` | bounds of the columns (variables) | `mat['clb']`, `mat['cub']` |
| `obj` | objective coefficients | `mat['obj_cf']` |
| `lpdir` | direction, `"max"` or `"min"`. Default `"max"`. | `mat['obj_dir']` |
| `cnames` | variable names. Default: the reaction ids of `model`. | `mat['cnames']` |
| `rnames` | constraint names. Default: the metabolite ids of `model`. | `mat['rnames']` |
| `pname` | name of the problem object | |
| `scaling` | scaling of the constraint matrix (see `scaleProb`). Default: no scaling. | |
| `writeProbToFileName` | file name to write the problem in LP format | `open(f, 'w').write(Res['prob'].to_lp())` |
| `...` | passed to the `initialize` method of `sysBiolAlg`: `solver`, `method`, `solverParm` | `interface` argument of `solve_ccFBA_mat` |

## Details

The problem object is built to perform cost constrained flux balance analysis with a given model.
This is a linear programming problem:

```
max  c^T v
s.t. S v = 0
     a_i <= v_i <= b_i      for i = 1, ..., n
     + the ccFBA constraints (kcat, complexes, multi-functional enzymes, budget)
```

`S` is the stoichiometric matrix, `a_i` and `b_i` are the lower and upper bounds of flux
(variable) `i`, and `n` is the number of variables. The solution is a flux distribution that
maximizes the objective `c^T v` in the given environment, at steady state.

## Example

R:

```r
prob <- sysBiolAlg(mod2, LHS=mrMat$LHS, rlb=mrMat$rlb, rub=mrMat$rub, rtype=mrMat$rtype,
                   lb=mrMat$clb, ub=mrMat$cub, obj=mrMat$obj_cf, lpdir='max',
                   cnames=mrMat$cnames, algorithm="ccFBA", solver="glpkAPI")
sol <- optimizeProb(prob)
```

Python:

```python
import cobra
import pandas as pd
from ccfbapy import getccFBA_mat, solve_ccFBA_mat

model = cobra.io.read_sbml_model('models/iAF1260.xml')
mw = pd.read_csv('data/Ecoli_MW.csv')
kcat = pd.read_csv('data/Ecoli_kcat.csv')
mrMat = getccFBA_mat(model, Kcat=kcat, MW=mw, RHS=0.27, medval=3600*22.6, verboseMode=0)
sol = solve_ccFBA_mat(mrMat)          # uses mrMat['model'] for the flux bounds
print(sol['status'], sol['obj'])
prob = sol['prob']                    # optlang problem: add constraints, change bounds, ...
```

## References

Gelius-Dietrich G, Desouki AA, Fritzemeier CJ, Lercher MJ (2013). sybil – Efficient
constraint-based modelling in R. *BMC Systems Biology* 7(1): 125.

Adadi R, et al. (2012), see [`cfba_moment`](cfba_moment.md#references).
