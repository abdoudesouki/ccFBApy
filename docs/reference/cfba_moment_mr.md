# cfba_moment_mr: MOMENT with multi-functional enzymes (`ccFBA_mr`)

*R help page: `man/cfba_moment_mr.Rd`. Python: `ccfbapy/ccFBA_mr.py`, function `ccFBA_mr()` (the R name is
`cfba_moment_mr`).*

## Description

Computes fluxes like [`cfba_moment`](cfba_moment.md) (GPR, kcat, molecular weights), but takes
multi-functional enzymes into account. In MOMENT, a protein `i` that catalyzes several reactions can
use its whole concentration `g_i` for each of them. Here, when protein `i` takes part in more than one
reaction, there is an auxiliary concentration variable `x_ij` for each of these reactions `j`. The
`x_ij` replace `g_i` in the kcat constraint of reaction `j`. Their sum is bounded by the total
concentration `g_i`, which is the variable in the budget constraint.

## Usage

```python
from ccfbapy import ccFBA_mr

ccFBA_mr(model, Kcat=None, MW=None, selected_rxns=None, verboseMode=2, objVal=None,
         RHS=None, solver="GLPK", medval=None, geneConcLB=None, geneConcUB=None)
```

R:

```r
cfba_moment_mr(model, mod2=NULL, Kcat, MW=NULL, selected_rxns=NULL, verboseMode=2, objVal=NULL,
               RHS=NULL, solver=SYBIL_SETTINGS("SOLVER"), C_mu_coef=0, medval=NULL,
               runFVA=FALSE, fvaRxn=NULL,
               geneConcLB=rep(0,length(allGenes(model))), geneConcUB=rep(1000,length(allGenes(model))))
```

## Arguments

| Argument | Description |
|---|---|
| `model` | `cobra.Model`. It is copied, so the input model is not changed. The bounds of its reactions are used. |
| `Kcat` | kcat values in 1/s: DataFrame `rxn_id`, `val`, `dirxn`. `val` must be the second column and `rxn_id` must have the `R_` prefix, see [Kcat](../conventions.md#kinetic-data-kcat). |
| `MW` | molecular weights in g/mmol: first column gene, second column MW. |
| `selected_rxns` | list of reaction ids to constrain. Default: all reactions with a GPR, except those whose GPR is only `s0001`. See [Known issues](#known-issues). |
| `verboseMode` | amount of output, default `2`. |
| `objVal` | accepted for compatibility with R, not used. |
| `RHS` | the budget C in g/gDW, e.g. 0.27. When `None` no budget constraint is added. |
| `solver` | accepted for compatibility with R, not used. |
| `medval` | kcat in **1/h** for reactions without a kcat. Default: median of `Kcat['val']`. |
| `geneConcLB`, `geneConcUB` | `dict` gene id → bound of the total concentration `g_<gene>` (mmol/gDW). Defaults 0 and 1000. |

The R arguments `mod2`, `C_mu_coef`, `runFVA` and `fvaRxn` are not available. For `C_mu_coef`,
use [`getccFBA_mat`](getccFBA_mat.md).

## Details

Main steps:

1. Add a variable `g_<gene>` for every gene that is in at least one selected reaction. If the gene
   is in `n > 1` selected reactions, also add one variable `gr_<gene>_<rxn>` for each of them.
2. For each selected reaction, parse the GPR rule.
3. Add the kcat constraints as in [`cfba_moment`](cfba_moment.md#details). For a
   multi-functional gene, the constraint uses `gr_<gene>_<rxn>` instead of `g_<gene>`.
4. Add the multi-functional enzyme constraints `mfe_<gene>`: `sum_rxn gr_<gene>_<rxn> - g_<gene> <= 0`.
5. Add the budget constraint `MW`: `sum(MW_g * g_g) <= RHS`.

In R, the `mfe_` constraints are equalities (`sum = g`). In Python they are `<=`. The growth
rate is the same, but with a lower bound on `g_<gene>` (over expression) the Python version
can leave part of the protein unused.

With the same `RHS`, `ccFBA_mr` predicts lower growth rates than `cfba_moment`, because the capacity
of multi-functional enzymes is shared between their reactions (iAF1260 on glucose in the test
notebooks: 0.26 vs 0.47 1/h).

## Value

A `dict`:

| Key | Content |
|---|---|
| `obj` | objective value |
| `status` | solver status |
| `rxnFlx` | `dict` reaction id → net flux |
| `geneConc` | `dict` gene id → total concentration `g_<gene>` (mmol/gDW) |
| `rxnkcatRow` | DataFrame of the kcat terms, as in [`cfba_moment`](cfba_moment.md#value) |
| `prob` | the optlang problem |
| `model` | the copy of `model` that has the enzyme layer |

The R function returns more: `sol`, `geneCol`, `metflx`, `metflx_ir`, `actual_C` (budget used),
`geneConc` (a data frame, with a row for the total and one per reaction for multi-functional
genes), `rxnMC` (protein cost of each reaction; they sum to the budget), `rxnGeneMC` (contribution
of each gene to each reaction), `lp_Mat`, `cNames`, `rNames`, and `fva` when `runFVA=TRUE`. To get
the protein cost of each reaction in Python, use [`getccFBA_sol`](getccFBA_sol.md) (`flx_MC`,
`rg_MC`).

## Known issues

- With `selected_rxns`, the problem is solved, but collecting `geneConc` fails with a `KeyError`
  for the genes that are not in any selected reaction (the result cannot be returned).

## References

Adadi R, et al. (2012), see [`cfba_moment`](cfba_moment.md#references).

Desouki AA (2016). *Algorithms for improving the predictive power of flux balance analysis.* PhD
dissertation.

## See also

[`cfba_moment`](cfba_moment.md), [`getccFBA_mat`](getccFBA_mat.md) (multi-functional enzymes
with the cheapest isoenzyme and complex stoichiometry), [`simGeneOverExpr`](simGeneOverExpr.md).

## Examples

Python version of the example of the R help page (iAF1260 from the R data set, glucose):

```python
from ccfbapy import datasets
from ccfbapy import ccFBA_mr

model = datasets.load_model('iAF1260')       # R: data(iAF1260)
mw = datasets.load_dataset('mw')             # R: data(mw)
kcat = datasets.load_dataset('kcat')         # R: data(kcat)

model.reactions.get_by_id('R_EX_glc_e_').lower_bound = -1000
model.reactions.get_by_id('R_EX_glyc_e_').lower_bound = 0
model.reactions.get_by_id('R_EX_ac_e_').lower_bound = 0
model.reactions.get_by_id('R_EX_o2_e_').lower_bound = -1000
model.reactions.get_by_id('R_ATPM').lower_bound = 0

sol_mr = ccFBA_mr(model, kcat, MW=mw, verboseMode=1, RHS=0.27, medval=3600*22.6)
bm_rxn = [r.id for r in model.reactions if r.objective_coefficient != 0][0]
print('biomass=%f' % sol_mr['rxnFlx'][bm_rxn])

# enzyme concentrations and the parts of a multi-functional enzyme
gconc = sol_mr['geneConc']
prob = sol_mr['prob']
parts = {v.name: v.primal for v in prob.variables if v.name.startswith('gr_b2323_') and v.primal > 1e-9}
print(gconc['b2323'], parts)
```

Growth on 24 carbon sources and the comparison with measured growth rates: `test_ccFBA_mr.ipynb`.
