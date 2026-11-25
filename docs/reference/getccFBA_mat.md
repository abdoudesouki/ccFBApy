# getccFBA_mat: ccFBA problem as a constraint matrix

*R help page: `man/getccFBA_mat.Rd`. Python: `ccfbapy/getccFBA_mat.py`, functions `getccFBA_mat()` and
`solve_ccFBA_mat()`.*

## Description

Given kcat values (1/s), molecular weights and, optionally, the stoichiometry of protein
complexes, builds a ccFBA problem (maximize the objective, e.g. biomass) as a constraint matrix.
Improvements over MOMENT:

1. multi-functional enzymes: a gene that catalyzes `n > 1` reactions gets one variable per
   reaction, and their sum equals the gene concentration;
2. one enzyme (the cheapest isoenzyme, see [`getGpr1iso`](getGpr1iso.md)) per reaction. This
   removes the under estimation of the cost of rules like `(A and B) or (A and C)`;
3. optional protein complex stoichiometry.

[`ccFBA_iso`](ccFBA_iso.md) removes the under estimation of point 2 without removing isoenzymes, by
allocating each gene to the complexes that use it.

The matrix is built once and can be solved many times with other flux bounds, e.g. for several
carbon sources (see [`simulate_EColi`](simulate_EColi.md)). Use
[`solve_ccFBA_mat`](#solve_ccfba_mat) to solve it.

## Usage

```python
from ccfbapy import getccFBA_mat, solve_ccFBA_mat

getccFBA_mat(model, mod2=None, Kcat=None, MW=None, selected_rxns=None, verboseMode=2, objVal=None,
             RHS=None, solver=None, medval=None, geneConcLB=None, geneConcUB=None,
             cpx_stoich=None, C_mu_coef=0)
```

R:

```r
getccFBA_mat(model, mod2=NULL, Kcat, MW=NULL, selected_rxns=NULL, verboseMode=2, objVal=NULL,
             RHS=NULL, solver=SYBIL_SETTINGS("SOLVER"), medval=NULL,
             geneConcLB=rep(0,length(allGenes(model))), geneConcUB=rep(1000,length(allGenes(model))),
             cpx_stoich=NULL, C_mu_coef=0)
```

## Arguments

| Argument | Description |
|---|---|
| `model` | `cobra.Model` |
| `mod2` | irreversible model, `mod2irrev(model)` (see [mod2irrev](mod2irrev.md)). Give it to avoid building it again, or to split more reactions (`rev=`). It must be made from `model`. Default: `mod2irrev(model)`. |
| `Kcat` | **required**. kcat values in 1/s: DataFrame with columns `rxn_id`, `val`, `dirxn`, optionally `src`. Reaction ids with or without the `R_` prefix. |
| `MW` | **required** (it is used to choose the isoenzymes). Molecular weights in g/mmol: first column gene, second column MW. |
| `selected_rxns` | reaction ids of `model` to constrain. They must have a GPR. Default: all reactions with a GPR (including `s0001`, which has no cost). |
| `verboseMode` | `0` silent, `1` status, `2` more details. Default `2`. |
| `objVal` | used when `RHS` is `None`: the objective value to reach. Default: FBA optimum of `mod2`. |
| `RHS` | the budget C in g/gDW, e.g. 0.27. When `None`, the problem is to find the minimum budget that gives `objVal`. |
| `solver` | accepted for compatibility with R, not used. |
| `medval` | kcat in **1/h** for reactions without kcat. Default: median of `Kcat['val']`. |
| `geneConcLB`, `geneConcUB` | gene concentration bounds (mmol/gDW): `dict` gene → value, or a list in the order of `model.genes`. Defaults 0 and 1000. |
| `cpx_stoich` | stoichiometry of the complexes: DataFrame with the columns `genes` and `stoich` (see [conventions](../conventions.md#protein-complex-stoichiometry-cpx_stoich)). Default `None`: all coefficients are one. |
| `C_mu_coef` | the budget as a linear function of the growth rate `mu`: `C = RHS + C_mu_coef * mu`. Default 0 (not used). |

## Details

The matrix is `[S(mod2) 0; ccFBA rows]`. The columns are the fluxes of `mod2` followed by the enzyme
variables.

**Columns (variables)**

| Name | Meaning |
|---|---|
| `R_<rxn>` | flux of reaction `<rxn>` of `mod2` (bounds from `mod2`) |
| `g1_<gene>_1` | gene that catalyzes one reaction (one-to-one) |
| `gmr_<gene>_<n>` | gene that catalyzes `n` reactions: total concentration |
| `gr_<gene>_<rxn>` | part of the multi-functional gene `<gene>` that catalyzes `<rxn>` |
| `cpx_<rxn>` | the complex that catalyzes `<rxn>` |

**Rows (constraints)**, row types as in sybil: `E` equal, `U` `<= rub`, `L` `>= rlb`.

| Name | Type | Constraint |
|---|---|---|
| `S_<met>` | E | `S v = 0` |
| `rkcat_<rxn>` | U | `v - kcat/st * enz <= 0`, with `enz` the `g1_`, `gr_` or `cpx_` column of the reaction |
| `rbkcat_<rxn>_b` | U | the same for the backward direction, with the reverse kcat |
| `min_<gene>_<rxn>_<stoich>` | L | `gene - stoich * cpx >= 0` for each gene of a complex |
| `mfe_<gene>` | E | `sum_rxn gr_<gene>_<rxn> - gmr_<gene>_<n> = 0` |
| `cc` | U | `sum(MW_g * g_g) - C_mu_coef * obj <= RHS` (when `RHS` is given) |
| `ObjC_` | L | `obj >= objVal` (when `RHS` is `None`) |

`st` is the number of subunits when the enzyme is one gene (`v <= kcat * g / n` for an enzyme with `n`
identical subunits) and 1 for complexes. The stoichiometry of a complex is in the `min_` rows.

Only the genes of the chosen isoenzymes of the selected reactions get a column. Genes without a
molecular weight are not in the `cc` row.

The objective is the objective of `mod2` (direction `mod2.obj_dir`). When `RHS` is `None`, the
objective is to minimize `sum(MW_g * g_g)`.

## Value

A `dict`. All row and column indices are **0-based**.

| Key | Content |
|---|---|
| `LHS` | constraint matrix (`scipy.sparse`, rows × columns): S matrix and ccFBA constraints |
| `rlb`, `rub`, `rtype`, `rnames` | bounds, types and names of the rows |
| `clb`, `cub`, `cnames` | bounds and names of the columns |
| `obj_cf`, `obj_dir` | objective coefficients and direction (`'max'`/`'min'`) |
| `fluxind` | columns of the fluxes |
| `geneind` | Series gene → column of the gene (`<NA>` for genes without a column) |
| `rxnkcatRow` | kcat constraints: reaction, `vind`/`vind_b` (flux columns), `exprCol` (enzyme column), `kval`/`rkval` (kcat, 1/h), `r1gstoich`, `rowind`/`rowind_b`, `consNames`/`consNames_b` |
| `minGiRow` | complex constraints (`min_` rows) |
| `mfeRow` | multi-functional enzyme constraints (`mfe_` rows) |
| `rgst`, `mod_cpx_mw`, `gpr1iso` | output of [`getGpr1iso`](getGpr1iso.md) (`gpr1iso` with the column `exprCol`) |
| `rxnExprCol` | `flxCol`, `exprCol`, `rowind` of each kcat constraint |
| `rxnGeneCol` | reaction-gene pairs with their columns, kcat, MW and stoichiometry: used to compute the protein cost (molecular crowding) |
| `geneCol`, `mfeCol` | gene columns and multi-functional gene parts (low level, for debugging) |
| `kcat_ir`, `kcat_src_ir` | kcat (1/h) and its source for each reaction of `mod2` (`medianVal` when missing) |
| `geneMW` | molecular weight of each gene of the model (NaN when missing) |
| `model` | the irreversible model `mod2` |
| `cc_ind` | row of the capacity constraint (`None` when `RHS` is `None`) |
| `mod_objc_ind` | row of the constraint on the objective (`None` when `RHS` is given) |

## solve_ccFBA_mat

```python
solve_ccFBA_mat(mat, mod2=None, interface=None, verboseMode=1, prob=None)
```

Solves the problem returned by `getccFBA_mat` with optlang. In R this is
`sysBiolAlg(mod2, algorithm="ccFBA", LHS=..., ...)` followed by `optimizeProb()` (see
[sysBiolAlg_ccFBA](sysBiolAlg_ccFBA.md)).

| Argument | Description |
|---|---|
| `mat` | result of `getccFBA_mat` |
| `mod2` | irreversible model whose `lowbnd`/`uppbnd` are used for the flux columns. Default `mat['model']`. Change the bounds there before solving. |
| `interface` | optlang interface, e.g. `optlang.glpk_interface`. Default: cobra's default solver. |
| `verboseMode` | `0` silent |
| `prob` | the problem of an earlier call (`Res['prob']`) with the same `mat`. Only the flux bounds are updated, which is faster when many conditions are solved. |

Returns a `dict`: `obj`, `status`, `x` (values of all columns, NaN when not optimal), `rxnFlx`
(net flux of each reaction of the original model), `flx_ir` (fluxes of `mod2`), `geneConc` (gene →
value of its `g1_`/`gmr_` column), `prob` (the optlang problem).

## Differences from R

- The bounds of the fluxes must be set in `mod2` (or in `model` before `mod2irrev`). A reaction
  that was not split (lower bound ≥ 0 when `mod2` was made) cannot carry a backward flux later.
  Open the uptake of every carbon source you want to use before building `mod2`, or give these
  reactions in `mod2irrev(model, rev=[...])`.
- When `RHS` is `None`, the constraint on the objective is written to a new row of type `L`.
  R writes it to a wrong row with type `U`, so it has no effect there.
- Reactions that are not in `selected_rxns` get no enzyme column (`exprCol` is NA).
- `Kcat` does not need a `src` column.

## References

Desouki AA (2016). *Algorithms for improving the predictive power of flux balance analysis.* PhD
dissertation.

## See also

[`getGpr1iso`](getGpr1iso.md), [`ccFBA_iso`](ccFBA_iso.md), [`cfba_moment_mr`](cfba_moment_mr.md),
[`getccFBA_sol`](getccFBA_sol.md), [`simulate_EColi`](simulate_EColi.md),
[`mod2irrev`](mod2irrev.md).

## Examples

Build the matrix once, then solve it for two carbon sources (see also `test_ccFBA_mat.ipynb`):

```python
import cobra
import pandas as pd
from ccfbapy import getccFBA_mat, solve_ccFBA_mat

model = cobra.io.read_sbml_model('models/iAF1260.xml')
mw = pd.read_csv('data/Ecoli_MW.csv')
kcat = pd.read_csv('data/Ecoli_kcat.csv')

CSList = ['EX_glc__D_e', 'EX_glyc_e']
for rid in CSList:                      # open the uptake before the model is split
    model.reactions.get_by_id(rid).bounds = (-1000, 0)
model.reactions.get_by_id('EX_o2_e').lower_bound = -1000
model.reactions.get_by_id('ATPM').bounds = (0, 0)

mat = getccFBA_mat(model, Kcat=kcat, MW=mw, RHS=0.27, medval=3600*22.6, verboseMode=1)
mod2 = mat['model']
ub0 = mod2.uppbnd.copy()
prob = None
for cs in CSList:
    mod2.uppbnd[:] = ub0
    for other in CSList:                # close the uptake of the other carbon sources
        if other != cs:
            mod2.uppbnd[mod2.react_id.index(other + '_b')] = 0
    Res = solve_ccFBA_mat(mat, verboseMode=0, prob=prob)
    prob = Res['prob']
    print(cs, Res['status'], round(Res['obj'], 4), 'uptake:', round(-Res['rxnFlx'][cs], 3))
mod2.uppbnd[:] = ub0
```

With protein complex stoichiometry:

```python
cpx_stoich = pd.read_csv('R/sybilccFBA/inst/extdata/cpx_stoich_me.csv')
mat_st = getccFBA_mat(model, mod2=mod2, Kcat=kcat, MW=mw, RHS=0.27, medval=3600*22.6,
                      cpx_stoich=cpx_stoich, verboseMode=0)
print(mat_st['minGiRow'].loc[mat_st['minGiRow']['stoich'] > 1, ['rxn', 'gene', 'stoich']].head())
```

Minimum budget needed for a growth rate of 0.2 1/h (`RHS=None`). Without `objVal` the FBA optimum
of `mod2` is used. Here the uptake is not limited, so that optimum is not a realistic growth rate.

```python
mat_min = getccFBA_mat(model, mod2=mod2, Kcat=kcat, MW=mw, RHS=None, objVal=0.2, medval=3600*22.6,
                       verboseMode=0)
print('minimum budget (g/gDW):', solve_ccFBA_mat(mat_min, verboseMode=0)['obj'])   # about 0.2
```
