# simulate_EColi: growth of *E. coli* on several carbon sources

*R help page: `man/simulate_EColi.Rd`. Python: `ccfbapy/simulate_EColi.py`, function `simulate_EColi()`.*

## Description

Shows how to use the ccFBA functions to simulate the metabolism of an *E. coli* model on a list of
carbon sources. The matrix is built once with [`getccFBA_mat`](getccFBA_mat.md). For each carbon
source, the uptake of all the other carbon sources is closed and the problem is solved.

## Usage

```python
from ccfbapy import simulate_EColi

simulate_EColi(model, mod2, kcat, mw, budget_C, CSList, atpmz=True, trns_rxns=None,
               cpx_stoich=None, solver=None, verboseMode=1)
```

R:

```r
simulate_EColi(model, mod2, kcat, mw, budget_C, CSList, atpmz=TRUE, trns_rxns=NULL,
               cpx_stoich=NULL, solver="glpkAPI")
```

## Arguments

| Argument | Description |
|---|---|
| `model` | `cobra.Model` |
| `mod2` | irreversible model, `mod2irrev(model, rev=...)` (see [mod2irrev](mod2irrev.md)). The uptake reactions in `CSList` must exist in it, so the exchange reactions of the carbon sources must be split. It is not modified (a copy is used). |
| `kcat` | kcat values in 1/s: DataFrame `rxn_id`, `val`, `dirxn`, optionally `src`. Reactions without kcat get `3600 * median(kcat['val'])`. |
| `mw` | molecular weights in g/mmol: first column gene, second column MW |
| `budget_C` | the budget C, 0.27 for *E. coli* |
| `CSList` | carbon sources to simulate: ids of reactions of `mod2` (the uptake direction, e.g. `'EX_glc__D(e)_b'`) |
| `atpmz` | set the ATPM reaction to zero. Default `True`. Only a reaction with the id `ATPM` (or `ATPM_f`/`ATPM_b`) is found. |
| `trns_rxns` | transport reactions to leave out of the budget (may also contain other reactions): booleans in the order of `model.reactions`, or a list of reaction ids |
| `cpx_stoich` | stoichiometry of the complexes: DataFrame with the columns `genes` and `stoich`. Default `None`: all coefficients are one. |
| `solver` | accepted for compatibility with R, not used (cobra's default solver is used) |
| `verboseMode` | `0` prints less |

## Details

The budget constraint covers all reactions with a GPR except `trns_rxns`. For each carbon source
`cs` (numbered from 1, as in R), the upper bounds of all reactions in `CSList` are set to 0, the
bound of `cs` is set to 1000, and the problem is solved. The problem is built once and reused,
so only the bounds change.

## Value

A `dict`:

| Key | Content |
|---|---|
| `gr_data` | rows of `all_flx` for the objective reaction: the growth rate (column `flx`) on each carbon source |
| `uptake` | rows of `all_flx` for the uptake reaction of each carbon source |
| `all_flx` | `cs`, `csname`, `rxn_id`, flux (`flx`), `kcat`, `kcat_src`, enzyme concentration (`enzConc`), `obj`, `ub`, `ubR` for all reactions of `mod2` and all carbon sources |
| `rg_MC` | reaction-gene molecular crowding (MC): `all_flx_MC` sorted by decreasing `rxnMC` for each carbon source (columns `rxnMCorder`, `rgMCorder`) |
| `all_flx_MC` | reaction-gene pairs with `flx`, `grConc`, `enzConc`, `geneConc`, `rgMC` (MC of the gene part) and `rxnMC` (MC of the reaction), not sorted. See [`flx_MC_tables`](getccFBA_sol.md#flx_mc_tables). |
| `mrMat` | the matrix from `getccFBA_mat` (not returned in R) |

In R the enzyme concentration is divided by the stoichiometry of each gene. In Python it is
`x[exprCol] / r1gstoich`, as in the kcat constraint (see [`flx_MC_tables`](getccFBA_sol.md#flx_mc_tables)).

## References

Desouki AA (2016). *Algorithms for improving the predictive power of flux balance analysis.* PhD
dissertation.

## See also

[`getccFBA_mat`](getccFBA_mat.md), [`cfba_moment_mr`](cfba_moment_mr.md),
[`getccFBA_sol`](getccFBA_sol.md).

## Examples

Python version of the example of the R help page: iML1515, kcat_max of `allKcats_upd34_dd_h.csv`,
complex stoichiometry, C = 0.27, transport reactions of the outer membrane (`tex`) left out of the
budget. The same code, with plots and tables, is in `test_simulate_EColi.ipynb` (Spearman
correlation with the measured growth rates 0.37, Pearson 0.47).

```python
import os
import pandas as pd
from scipy.stats import spearmanr
from ccfbapy import datasets
from ccfbapy import mod2irrev
from ccfbapy import simulate_EColi

pkg = os.path.join('R', 'sybilccFBA')
# 1- get metabolic model
model = datasets.load_model('iML1515')            # R: data(iML1515)
# 2- get molecular weights
mw = datasets.load_dataset('mw_iML1515')          # R: data(mw_iML1515)
mw = pd.concat([mw, pd.DataFrame({'Synonym': ['s0001'], 'mw': [0.001]})], ignore_index=True)
mw.columns = ['gene', 'mw']
# 3- get kcat list
kl = pd.read_csv(os.path.join(pkg, 'inst', 'extdata', 'allKcats_upd34_dd_h.csv'))
kl = kl[kl['ijo_id'].notna()]
kcat = pd.DataFrame({'rxn_id': kl['ijo_id'], 'val': kl['kcat_max'], 'dirxn': kl['dirxn'], 'src': kl['src']})
kcat = kcat[kcat['rxn_id'].isin([r.id for r in model.reactions])]
kcat.loc[kcat['src'].isna(), 'src'] = 'Max'
# 4- get complex stoichiometry if used
cpx_stoich = pd.read_csv(os.path.join(pkg, 'inst', 'extdata', 'cpx_stoich_me.csv'))
# 5- identify carbon sources to be tested
csl = ["EX_glc__D(e)_b", "EX_glyc(e)_b", "EX_ac(e)_b", "EX_fru(e)_b", "EX_pyr(e)_b", "EX_gal(e)_b",
       "EX_lac__L(e)_b", "EX_malt(e)_b", "EX_mal__L(e)_b", "EX_fum(e)_b", "EX_xyl__D(e)_b", "EX_man(e)_b",
       "EX_tre(e)_b", "EX_mnl(e)_b", "EX_g6p(e)_b", "EX_succ(e)_b", "EX_gam(e)_b", "EX_sbt__D(e)_b",
       "EX_glcn(e)_b", "EX_rib__D(e)_b", "EX_gsn(e)_b", "EX_ala__L(e)_b", "EX_akg(e)_b", "EX_acgam(e)_b"]
msrd = [0.66, 0.47, 0.29, 0.54, 0.41, 0.24, 0.41, 0.52, 0.55, 0.47, 0.51, 0.35, 0.48, 0.61,
        0.78, 0.50, 0.40, 0.48, 0.68, 0.41, 0.37, 0.24, 0.24, 0.61]
# 6- get irreversible model: split the exchange reactions of the carbon sources
mod2 = mod2irrev(model, rev=[c[:-2] for c in csl])
for c in csl:                                   # close uptake and secretion of all carbon sources
    mod2.uppbnd[mod2.react_id.index(c)] = 0
    mod2.uppbnd[mod2.react_id.index(c[:-2] + '_f')] = 0
mod2.uppbnd[mod2.react_id.index('EX_o2(e)_b')] = 1000
trns_rxns = [r.id.endswith('tex') for r in model.reactions]

## Call function
tmp_res = simulate_EColi(model, mod2, mw=mw, budget_C=0.27, kcat=kcat, cpx_stoich=cpx_stoich,
                         atpmz=False, trns_rxns=trns_rxns, CSList=csl)
bm = tmp_res['gr_data']
print(spearmanr(bm['flx'], msrd))
```
