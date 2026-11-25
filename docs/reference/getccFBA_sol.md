# getccFBA_sol: solve the ccFBA problem

*R help page: `man/getccFBA_sol.Rd`. Python: `ccfbapy/getccFBA_sol.py`, functions `getccFBA_sol()` and
`flx_MC_tables()`.*

(The R help page has the title and description of `getccFBA_mat`. This function builds the ccFBA
matrix with [`getccFBA_mat`](getccFBA_mat.md) and solves it.)

## Description

Formulates the ccFBA problem (maximize biomass) from the kcat values (1/s), molecular weights and,
optionally, the stoichiometry of protein complexes, and solves it. It uses the same improvements
over MOMENT as [`getccFBA_mat`](getccFBA_mat.md): multi-functional enzymes, the cheapest
isoenzyme and complex stoichiometry. It can also take a matrix built earlier (`mrMat`), so the
matrix is not built again in each call. It returns the fluxes, gene concentrations and protein cost
of each reaction.

## Usage

```python
from ccfbapy import getccFBA_sol

getccFBA_sol(model, mod2=None, Kcat=None, MW=None, selected_rxns=None, verboseMode=2, objVal=None,
             RHS=None, solver=None, medval=None, geneConcLB=None, geneConcUB=None,
             cpx_stoich=None, C_mu_coef=0, mrMat=None, CSList=None, cs=1)
```

R:

```r
getccFBA_sol(model, mod2=NULL, Kcat, MW=NULL, selected_rxns=NULL, verboseMode=2, objVal=NULL,
             RHS=NULL, solver=SYBIL_SETTINGS("SOLVER"), medval=NULL,
             geneConcLB=rep(0,length(allGenes(model))), geneConcUB=rep(1000,length(allGenes(model))),
             cpx_stoich=NULL, C_mu_coef=0, mrMat=NULL, CSList=NULL, cs=1)
```

## Arguments

`model`, `mod2`, `Kcat`, `MW`, `selected_rxns`, `verboseMode`, `RHS`, `medval`, `geneConcLB`,
`geneConcUB` and `cpx_stoich` are passed to [`getccFBA_mat`](getccFBA_mat.md#arguments).

| Argument | Description |
|---|---|
| `mrMat` | a matrix from `getccFBA_mat`. When given, the matrix is not built again; only `geneConcLB`/`geneConcUB` are applied to it, and the flux bounds are taken from `mod2` (default `mrMat['model']`). |
| `objVal`, `C_mu_coef` | accepted, but **not passed** to `getccFBA_mat` (as in R). Call `getccFBA_mat` yourself and give the result as `mrMat` to use them. |
| `solver` | accepted for compatibility with R, not used. |
| `CSList`, `cs` | list of carbon sources and index of the current one. Accepted for compatibility with R. R only uses them in a table that it does not return. |

## Value

A `dict`:

| Key | Content |
|---|---|
| `sol` | result of [`solve_ccFBA_mat`](getccFBA_mat.md#solve_ccfba_mat) (`obj`, `status`, `x`, ...) |
| `obj`, `status` | objective value and solver status |
| `metflx` | Series: net flux of the reactions of `model` |
| `metflx_ir` | Series: fluxes of the irreversible model `mod2` |
| `rxnFlx` | `metflx` as a `dict` (same key as in `cfba_moment` and `ccFBA_mr`) |
| `geneConc` | `dict` gene → concentration (mmol/gDW), genes that have a column. R returns a data frame `gene`, `geneconc`. |
| `flx_MC` | protein cost (molecular crowding, MC) of the reaction-gene pairs, see `flx_MC_tables` below |
| `rg_MC` | `flx_MC` sorted by decreasing `rxnMC` |
| `mrMat` | the matrix |
| `prob` | the optlang problem (R returns it as `prop`) |

## flx_MC_tables

```python
from ccfbapy import flx_MC_tables
flx_MC, rg_MC = flx_MC_tables(mrMat, x)
```

Protein cost of the reaction-gene pairs for a solution `x` (values of all columns) of the matrix
`mrMat`. `flx_MC` has the columns of `mrMat['rxnGeneCol']` and:

| Column | Content |
|---|---|
| `flx` | flux of the reaction (of `mod2`) |
| `grConc` | concentration of the gene part used by the reaction (mmol/gDW) |
| `enzConc` | enzyme concentration of the reaction: `x[exprCol] / r1gstoich`. NaN for the genes of isoenzymes that were not chosen. |
| `geneConc` | total concentration of the gene |
| `rgMC` | `gMW * grConc`: cost of the gene part (g/gDW) |
| `rxnMC` | MW of the enzyme × `enzConc`: cost of the reaction (g/gDW) |

`rg_MC` is `flx_MC` sorted by decreasing `rxnMC`, with the columns `rxnMCorder` and `rgMCorder`
(order by `rgMC`).

In R the enzyme concentration is divided by the stoichiometry of each gene. That gives `Inf`
for the genes of isoenzymes that were not chosen, and a wrong cost for complexes with
stoichiometry. In Python the costs of the reactions add up to the budget when each reaction of the
original model is counted once. The `_f` and `_b` rows of a reaction, and the rows of the genes
of a complex, share the same enzyme:
`flx_MC.groupby('revRxn')['rxnMC'].max().sum()` equals `RHS` when the budget is fully used.

## Details

Variable names (columns):

- `g1_x_1`: gene `x` catalyzes one reaction (one-to-one);
- `gr_x_r`: gene `x` catalyzes more than one reaction; this variable is the part used by reaction `r`;
- `gmr_x_cnt`: gene `x` catalyzes `cnt` reactions; this variable is the sum of the parts;
- `cpx_r`: the complex that catalyzes reaction `r`.

## References

Desouki AA (2016). *Algorithms for improving the predictive power of flux balance analysis.* PhD
dissertation.

## See also

[`getccFBA_mat`](getccFBA_mat.md), [`getGpr1iso`](getGpr1iso.md),
[`simulate_EColi`](simulate_EColi.md), [`simGeneOverExpr`](simGeneOverExpr.md).

## Examples

```python
import cobra
import pandas as pd
from ccfbapy import getccFBA_mat, mod2irrev
from ccfbapy import getccFBA_sol

model = cobra.io.read_sbml_model('models/iAF1260.xml')
mw = pd.read_csv('data/Ecoli_MW.csv')
kcat = pd.read_csv('data/Ecoli_kcat.csv')
mod2 = mod2irrev(model)

res = getccFBA_sol(model, mod2=mod2, Kcat=kcat, MW=mw, RHS=0.27, medval=3600*22.6, verboseMode=0)
print(res['status'], res['obj'])
print(res['rg_MC'][['rxn_id', 'gene', 'flx', 'enzConc', 'rxnMC', 'rgMC']].head())

# reuse the matrix, with a lower bound on the concentration of b2323
res2 = getccFBA_sol(model, mod2=mod2, RHS=0.27, mrMat=res['mrMat'], verboseMode=0,
                    geneConcLB={'b2323': 0.001})
print(res2['obj'], res2['geneConc']['b2323'])
```
