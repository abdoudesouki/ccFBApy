# cfba_moment: MOMENT method

*R help page: `man/cfba_moment.Rd`. Python: `ccfbapy/cfba_moment.py`, function `cfba_moment()`.*

## Description

Uses the GPR rules, kcat values and molecular weights to compute fluxes with the MOMENT method
(Adadi et al. 2012). Each reaction with a GPR gets a kcat constraint on its flux. The total
mass of the enzymes is limited by the budget `RHS`. The objective of the model (e.g. growth) is
maximized.

## Usage

```python
from ccfbapy import cfba_moment

cfba_moment(model, Kcat=None, MW=None, selected_rxns=None, verboseMode=2, objVal=None,
            RHS=None, solver="GLPK", medval=None, geneConcLB=None, geneConcUB=None)
```

R:

```r
cfba_moment(model, mod2=NULL, Kcat, MW=NULL, selected_rxns=NULL, verboseMode=2, objVal=NULL,
            RHS=NULL, solver=SYBIL_SETTINGS("SOLVER"), medval=NULL, runFVA=FALSE, fvaRxn=NULL,
            geneConcLB=rep(0,length(allGenes(model))), geneConcUB=rep(1000,length(allGenes(model))))
```

## Arguments

| Argument | Description |
|---|---|
| `model` | `cobra.Model`. It is copied, so the input model is not changed. The bounds of its reactions are used. |
| `Kcat` | kcat values in 1/s: DataFrame `rxn_id`, `val`, `dirxn`. `val` must be the second column and `rxn_id` must have the `R_` prefix, see [Kcat](../conventions.md#kinetic-data-kcat). |
| `MW` | molecular weights in g/mmol: DataFrame, first column gene, second column MW. See [MW](../conventions.md#molecular-weights-mw). |
| `selected_rxns` | list of reaction ids to constrain. Default: all reactions with a GPR, except those whose GPR is only `s0001`. |
| `verboseMode` | amount of output: `0` almost nothing, `1` status messages, `>2` each reaction. Default `2`. |
| `objVal` | accepted for compatibility with R, **not used** (see Details). |
| `RHS` | the budget C in g/gDW, 0.27 for *E. coli*. When `None` no budget constraint is added. |
| `solver` | accepted for compatibility with R, not used. The solver of `model` is used. |
| `medval` | kcat in **1/h** for reactions without a kcat, e.g. `3600*22.6`. Default: median of `Kcat['val']` (in 1/h). |
| `geneConcLB` | `dict` gene id → lower bound of the gene concentration (mmol/gDW). Default 0. |
| `geneConcUB` | `dict` gene id → upper bound of the gene concentration. Default 1000. |

The R arguments `mod2`, `runFVA` and `fvaRxn` are not available.

## Details

Main steps:

1. Add a concentration variable `g_<gene>` for every gene of the model.
2. For each selected reaction, parse the GPR rule.
3. Add variables and constraints according to the rule (table below).
4. Add the budget ("solvent capacity") constraint `MW`: `sum(MW_g * g_g) <= RHS`, for the genes
   that have a molecular weight.

Variables and constraints for each type of GPR rule:

| GPR rule | Variables | Constraints |
|---|---|---|
| one gene `A` | | `<rxn>_kcat`: `v <= kcat * g_A` |
| complex `A and B` | `cpx_<rxn>` | `<rxn>_kcat`: `v <= kcat * cpx`; `<rxn>_A`: `cpx <= g_A`; `<rxn>_B`: `cpx <= g_B` |
| isoenzymes `A or B` | | `<rxn>`: `v <= kcat * (g_A + g_B)` |
| `(A and B) or C` | `cpx_<rxn>_<p>` for AND term `p` (0-based) | `<rxn>`: `v <= kcat * (cpx_<rxn>_0 + g_C)`; `Aux_<n>`: `cpx <= g` for each gene of the term |

A reversible reaction (lower bound < 0) gets the same constraint for its reverse variable,
with the reverse kcat (`dirxn = -1`): `<rxn>_b_kcat` or `<rxn>_b`.

An OR with more than two terms is written as one constraint, not pairwise as in the original
MATLAB implementation of MOMENT (see [`cfba_moment_pw`](cfba_moment_pw.md), which also works on the
irreversible model and removes `s0001` from the rules).

The constraints are added to the solver problem of a copy of `model`. They are not added as
cobra reactions. If the copy already has an enzyme layer from an earlier call, it is removed first
(see [`remove_enzyme_layer`](utilities.md#remove_enzyme_layer)).

In R, when `RHS` is `NULL` the problem becomes: find the minimum budget that gives `objVal`
(default: FBA optimum). The Python version does not implement this. With `RHS=None` it
solves FBA with the kcat constraints only. Use [`getccFBA_mat`](getccFBA_mat.md) with
`RHS=None` for the minimum budget.

## Value

A `dict`:

| Key | Content |
|---|---|
| `obj` | objective value (growth rate, 1/h) |
| `status` | solver status (`'optimal'`) |
| `rxnFlx` | `dict` reaction id → net flux (forward − reverse) |
| `geneConc` | `dict` gene id → concentration (mmol/gDW), all genes of the model |
| `rxnkcatRow` | DataFrame, one row per kcat term: `rsn` (reaction index), `rxn_id`, `dirxn`, `rowind` (position of the constraint in `prob.constraints`; `None` for OR rules), `colind` (position of the enzyme variable in `prob.variables`), `val` (kcat, 1/h) |
| `prob` | the optlang problem (`model.solver` of the copy) |
| `model` | the copy of `model` that has the enzyme layer |

R returns `sol` (sybil `optsol`), `prob` (`sysBiolAlg`), `geneCol` (gene → column) and, since
version 3.0.2, `geneConc`.

## References

Adadi R, Volkmer B, Milo R, Heinemann M, Shlomi T (2012). Prediction of microbial growth rate versus
biomass yield by a metabolic network with kinetic parameters. *PLoS Comput Biol* 8(7).
doi:10.1371/journal.pcbi.1002575

Gelius-Dietrich G, Desouki AA, Fritzemeier CJ, Lercher MJ (2013). sybil – Efficient
constraint-based modelling in R. *BMC Systems Biology* 7(1): 125.

## See also

[`ccFBA_mr`](cfba_moment_mr.md) (multi-functional enzymes), [`getccFBA_mat`](getccFBA_mat.md),
[`simGeneOverExpr`](simGeneOverExpr.md), [conventions](../conventions.md).

## Examples

Python version of the first example of the R help page: iAF1260 from the R data set, growth on
glucose with oxygen, ATPM lower bound 0. In R the bounds are changed in `mod2`. In Python they are
changed in the cobra model.

```python
from ccfbapy import datasets
from ccfbapy import cfba_moment

model = datasets.load_model('iAF1260')       # R: data(iAF1260)
mw = datasets.load_dataset('mw')             # R: data(mw)
kcat = datasets.load_dataset('kcat')         # R: data(kcat)

model.reactions.get_by_id('R_EX_glc_e_').lower_bound = -1000   # uppbnd(mod2)[R_EX_glc_e__b] = 1000
model.reactions.get_by_id('R_EX_glyc_e_').lower_bound = 0      # uppbnd(mod2)[R_EX_glyc_e__b] = 0
model.reactions.get_by_id('R_EX_ac_e_').lower_bound = 0
model.reactions.get_by_id('R_EX_o2_e_').lower_bound = -1000
model.reactions.get_by_id('R_ATPM').lower_bound = 0

sol = cfba_moment(model, kcat, MW=mw, verboseMode=1, RHS=0.27, medval=3600*22.6)
bm_rxn = [r.id for r in model.reactions if r.objective_coefficient != 0][0]
print('biomass=%f' % sol['rxnFlx'][bm_rxn])

# fluxes and enzyme concentrations
import pandas as pd
rxnflx = pd.Series(sol['rxnFlx'], name='flx')
geneConc = pd.Series(sol['geneConc'], name='conc')
print(geneConc.sort_values(ascending=False).head())
```

Growth on several carbon sources (SBML model, BiGG ids). The full example with 24 carbon sources
and the comparison with the measured growth rates is in `test_moment.ipynb`.

```python
import cobra
import pandas as pd
from ccfbapy import cfba_moment

model = cobra.io.read_sbml_model('models/iAF1260.xml')
mw = pd.read_csv('data/Ecoli_MW.csv')
kcat = pd.read_csv('data/Ecoli_kcat.csv')

CSList = ['EX_glc__D_e', 'EX_glyc_e', 'EX_ac_e']
for rid in CSList:
    model.reactions.get_by_id(rid).bounds = (0, 0)
model.reactions.get_by_id('EX_o2_e').lower_bound = -1000
model.reactions.get_by_id('ATPM').bounds = (0, 0)

for cs in CSList:
    with model:                                   # bounds are restored after the block
        model.reactions.get_by_id(cs).lower_bound = -1000
        res = cfba_moment(model, kcat, MW=mw, verboseMode=0, RHS=0.27, medval=3600*22.6)
    print(cs, res['status'], round(res['obj'], 4), 'uptake:', round(-res['rxnFlx'][cs], 3))
```
