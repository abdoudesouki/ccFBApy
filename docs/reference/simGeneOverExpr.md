# simGeneOverExpr: over expression of one gene

*No R help page. R source: `R/OverExpr/simGeneOverExpr.R` (test: `R/OverExpr/test_simGeneOverExpr.R`,
kcat file `R/OverExpr/kcat2`). Python: `ccfbapy/simGeneOverExpr.py`, function `simGeneOverExpr()`.*

## Description

Simulates the over expression of one gene with ccFBA (Tin Y Pang). The gene `gx` must use at
least the fraction `y_pct` of the budget C. Then the fluxes of the enzymatic reactions are fixed
and the total flux is minimized.

## Usage

```python
from ccfbapy import simGeneOverExpr

simGeneOverExpr(gx, y_pct, model, mod2=None, Kcat=None, MW=None, selected_rxns=None, verboseMode=2,
                objVal=None, RHS=None, solver=None, medval=None, FUN=None)
```

R:

```r
simGeneOverExpr(gx, y_pct, model, mod2=NULL, Kcat, MW=NULL, selected_rxns=NULL, verboseMode=2,
                objVal=NULL, RHS=NULL, solver=SYBIL_SETTINGS("SOLVER"), medval=NULL,
                runFVA=FALSE, fvaRxn=NULL, FUN=cfba_moment)
```

## Arguments

| Argument | Description |
|---|---|
| `gx` | the gene to over express, e.g. `'b2323'` |
| `y_pct` | fraction of C that gene `gx` must use (at most 1), e.g. `0.05` for 5% |
| `FUN` | ccFBA function: [`cfba_moment`](cfba_moment.md) (default), [`ccFBA_mr`](cfba_moment_mr.md) or [`getccFBA_sol`](getccFBA_sol.md). It must return a `dict` with `obj`, `status`, `rxnFlx` and `geneConc`. |
| `model` | `cobra.Model` |
| `mod2` | irreversible model, only passed to `FUN` when `FUN` has this argument (`getccFBA_sol`). Default there: `mod2irrev(model)`. |
| `Kcat`, `MW`, `selected_rxns`, `verboseMode`, `medval` | passed to `FUN` |
| `RHS` | the budget C, **required** |
| `solver` | passed to `FUN` when not `None` |
| `objVal` | accepted for compatibility with R, not used |

## Details

Steps (as in R):

1. `WTSol`: call `FUN` without over expression. The protein cost (MC) of `gx` is its concentration × MW.
2. If `gx` already uses at least `y_pct` of C, `OExSol = WTSol` and `OverExprEff` is `False`.
   Otherwise `FUN` is called again with the lower bound `RHS * y_pct / MW` on the
   concentration of `gx` (`OverExprEff` is `True`).
3. Minimum total flux (sybil's algorithm `"mtf"`): the fluxes of the reactions that have at least
   one gene with a non-zero molecular weight, and of the objective reaction, are fixed to the
   fluxes of `OExSol`. Then the sum of the absolute fluxes of all reactions is minimized.

The fluxes are computed with the bounds of `model` (`cfba_moment`, `ccFBA_mr`) or of `mod2`
(`getccFBA_sol`). The minimum total flux step uses the bounds of `model`; the model is restored
afterwards.

If `gx` has no concentration variable in the problem of `FUN`, the over expressed protein cannot
catalyze any reaction. This happens with `getccFBA_sol` when the gene is not in the chosen
isoenzyme of any reaction. In that case the protein only takes its part of the budget:
`FUN` is called with `RHS * (1 - y_pct)`. R stops with an error in this case.

With `ccFBA_mr`, the sum of the parts of a multi-functional enzyme is `<=` its total
concentration in Python (`=` in R). Part of the over expressed protein may therefore stay
unused.

## Value

A `dict`:

| Key | Content |
|---|---|
| `Sol` | minimum total flux solution: `status`, `obj` (sum of absolute fluxes), `rxnFlx` |
| `OExSol` | result of `FUN` with over expression of `gx` |
| `WTSol` | wild type result, before over expression |
| `OverExprEff` | `True` when the over expression changed the problem (`gx` used less than `y_pct` of C) |
| `objVal` | flux of the objective reaction in `OExSol` (a `dict` when there are several) |
| `gxMC_WT`, `gxMC_OEx` | MC (g/gDW) of gene `gx` in `WTSol` and `OExSol` (not returned in R) |

## See also

[`cfba_moment`](cfba_moment.md), [`ccFBA_mr`](cfba_moment_mr.md), [`getccFBA_sol`](getccFBA_sol.md).
`test_simGeneOverExpr.ipynb` is the Python version of `R/OverExpr/test_simGeneOverExpr.R`.

## Examples

Over expression of b2323 (fabB) to 5–50% of C with MOMENT, iAF1260 on glucose:

```python
from ccfbapy import datasets
from ccfbapy import simGeneOverExpr
from ccfbapy import cfba_moment

kcat = datasets.load_dataset('kcat')         # R: data(kcat)
mw = datasets.load_dataset('mw')             # R: data(mw)
model = datasets.load_model('iAF1260')       # R: data(iAF1260)
model.reactions.get_by_id('R_EX_glc_e_').lower_bound = -1000
model.reactions.get_by_id('R_EX_o2_e_').lower_bound = -1000
model.reactions.get_by_id('R_ATPM').lower_bound = 0

for pct in [0.05, 0.25, 0.5]:
    RR = simGeneOverExpr('b2323', pct, model, Kcat=kcat, MW=mw, verboseMode=0, RHS=0.27,
                         medval=3600*22.6, FUN=cfba_moment)
    print(pct, RR['OverExprEff'], 'growth: %.4f' % RR['objVal'],
          'MC of b2323: %.4f -> %.4f' % (RR['gxMC_WT'], RR['gxMC_OEx']))
```

With `getccFBA_sol` (multi-functional enzymes, cheapest isoenzyme). The porin b0929 (ompF) is not
in the chosen isoenzyme of any reaction, so it has no concentration variable, and the budget of the
other proteins is reduced by 5%:

```python
from ccfbapy import mod2irrev
from ccfbapy import getccFBA_sol

mod2 = mod2irrev(model)
RR = simGeneOverExpr('b0929', 0.05, model, mod2, Kcat=kcat, MW=mw, verboseMode=1, RHS=0.27,
                     medval=3600*22.6, FUN=getccFBA_sol)
print(RR['OverExprEff'], RR['objVal'])
```
