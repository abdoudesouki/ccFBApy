# ccFBApy

Python implementation of the R package **sybilccFBA** (version 3.0.2), built on
[COBRApy](https://github.com/opencobra/cobrapy) / [optlang](https://github.com/opencobra/optlang)
instead of sybil. The original R package is in [`R/sybilccFBA`](R/sybilccFBA).

## Cost constrained flux balance analysis (ccFBA)

Flux balance analysis (FBA) finds fluxes in metabolic models at steady state
(Orth, Thiele & Palsson 2010). Cost constrained FBA also limits the fluxes by the
amount of enzyme the cell can afford. It needs kinetic data (kcat values), enzyme
molecular weights and the gene-protein-reaction (GPR) rules of the model:

- A flux can be at most kcat × enzyme concentration: `v_i <= kcat_i * g_i`
- The total enzyme mass is limited: `sum(MW_g * g_g) <= C`

The package implements two methods:

1. **MOMENT** (MetabOlic Modeling with ENzyme kineTics). Adadi R, Volkmer B, Milo R, Heinemann M,
   Shlomi T (2012). *Prediction of Microbial Growth Rate versus Biomass Yield by a Metabolic Network
   with Kinetic Parameters.* PLoS Comput Biol 8(7). [doi:10.1371/journal.pcbi.1002575](https://doi.org/10.1371/journal.pcbi.1002575)
2. **ccFBA with multi-functional enzymes**: an improvement of MOMENT in which an enzyme that catalyzes
   several reactions must share its concentration between them. Desouki, Abdelmoneim.
   *Algorithms for improving the predictive power of flux balance analysis.* PhD dissertation, 2016.

### How GPR rules become constraints

| GPR rule | Constraint |
|---|---|
| one gene `A` | `v <= kcat * g_A` |
| complex `A and B` | `cpx <= g_A`, `cpx <= g_B`, `v <= kcat * cpx` |
| isoenzymes `A or B` | `v <= kcat * (g_A + g_B)` |
| `(A and B) or C` | one auxiliary variable per AND term, summed as isoenzymes |

Reversible reactions get separate constraints for the forward and reverse flux. The reverse
direction uses the reverse kcat (`dirxn = -1`). Reactions without a measured kcat use `medval`.
Multiple ORs are evaluated in one constraint instead of pairwise as in the original MOMENT
implementation.

In **`ccFBA_mr`**, a gene `g` that catalyzes `n > 1` of the selected reactions gets one variable per
reaction, `g_r`, and `sum_r g_r <= g`. The molecular weight constraint is applied to `g`.

## Status of the Python implementation

| R (sybilccFBA) | Python | Status |
|---|---|---|
| `cfba_moment` | `cfba_moment.py`: `cfba_moment()` | Working, tested on iAF1260 ([`test_moment.ipynb`](test_moment.ipynb)) |
| `cfba_moment_mr` | `ccFBA_mr.py`: `ccFBA_mr()` | Working, tested on iAF1260 ([`test_ccFBA_mr.ipynb`](test_ccFBA_mr.ipynb)). See known issues for `selected_rxns` |
| `getGpr1iso` | `gpr1iso.py`: `getGpr1iso()` | Working, including complex stoichiometry (`cpx_stoich`) |
| `getccFBA_mat` | `getccFBA_mat.py`: `getccFBA_mat()` | Working, tested on iAF1260 ([`test_ccFBA_mat.ipynb`](test_ccFBA_mat.ipynb)). The same file has `mod2irrev()` (irreversible model) and `solve_ccFBA_mat()` (solves the matrix, replaces `sysBiolAlg`/`optimizeProb`) |
| `getRevFlux` | not needed | `rxnFlx` in the results already holds net fluxes (forward − reverse) |
| `simulate_EColi` | `simulate_EColi.py`: `simulate_EColi()` | Working, example of the R help page in [`test_simulate_EColi.ipynb`](test_simulate_EColi.ipynb) (iML1515) |
| `getccFBA_sol` | `getccFBA_sol.py`: `getccFBA_sol()` | Working: builds and solves the matrix; also returns `obj`, `status`, `rxnFlx`, `geneConc` (dict) like the other python functions |
| `simGeneOverExpr` (`R/Tin/simGeneOverExpr.R`) | `simGeneOverExpr.py`: `simGeneOverExpr()` | Working, test in [`test_simGeneOverExpr.ipynb`](test_simGeneOverExpr.ipynb) (python version of `R/Tin/test_simGeneOverExpr.R`) |
| `cfba_moment_pw`, `calc_MW`, `addGlcTrns`, `readmodel`, class `sysBiolAlg_ccFBA` | | Not ported yet |
| data sets (`data/*.RData`) | `sybil2cobra.py`: `read_RData()` | Reads the `.RData` files of the R package; sybil models become cobra models (ids kept as in sybil, e.g. `EX_glc__D(e)`) |
| | `genFuns.py`: `scatterPlot()` | Helper: measured vs predicted scatter plot with y = x line, linear fit and Pearson correlation |

Options of the R functions that are **not** available in Python yet: `runFVA`, `fvaRxn`, and in
`cfba_moment`/`ccFBA_mr` also `C_mu_coef` and `objVal` (accepted but unused; both work in `getccFBA_mat`). `solver` is also accepted but unused: the problem is solved
with the model's solver. Only GLPK, cobra's default, has been tested.

## Requirements

```
pip install cobra pandas numpy scipy matplotlib
pip install rdata    # only to read the .RData files of the R package (sybil2cobra.py)
```

Tested with Python 3.12, cobra 0.32.1, optlang 1.9.1 (GLPK solver). The notebooks were created with Python 3.11.

## Usage

```python
import cobra
import pandas as pd
import cfba_moment as cfba
import ccFBA_mr as ccfba

model = cobra.io.read_sbml_model('models/iAF1260.xml')
mw = pd.read_csv('data/ecoli_mw.csv')      # gene, molecular weight (g/mmol)
kcat = pd.read_csv('data/ecoli_kcat.csv')  # rxn_id, val (1/s), dirxn (1: forward, -1: reverse)

# MOMENT
Res = cfba.cfba_moment(model, kcat, MW=mw, verboseMode=1, RHS=0.27, medval=3600*22.6)
print(Res['obj'], Res['status'])

# ccFBA with multi-functional enzymes
Res_mr = ccfba.ccFBA_mr(model, kcat, MW=mw, verboseMode=1, RHS=0.27, medval=3600*22.6)
```

Run the scripts from the repository folder: `ccFBA_mr.py` imports `cfba_moment.py`, and `getccFBA_mat.py`
imports `gpr1iso.py`. The functions
work on a copy of the model, so the model you pass in is not modified. They can be called repeatedly
on the same model or on copies of it. If a model already contains the variables and constraints of
an earlier run, these are removed first.

### Parameters (`cfba_moment` and `ccFBA_mr`)

| Parameter | Description |
|---|---|
| `model` | cobra `Model` |
| `Kcat` | DataFrame with columns `rxn_id`, `val`, `dirxn`. `val` is in 1/s and converted to 1/h internally. `rxn_id` is the SBML id with the `R_` prefix (e.g. `R_ASPTA`); cobra removes this prefix from reaction ids, so the functions add it for the lookup (unless the id already starts with `R_`) |
| `MW` | DataFrame: first column gene id, second column molecular weight in g/mmol (kDa). Only genes with a molecular weight are included in the enzyme mass constraint |
| `RHS` | Upper limit `C` of enzyme mass in g protein/gDW (e.g. 0.27 for *E. coli*, Adadi et al. 2012). The mass constraint is only added when `RHS` is given |
| `medval` | kcat (**in 1/h**) used for reactions without a measured kcat. Default: median of `Kcat` |
| `selected_rxns` | list of reaction ids to constrain. Default: all reactions with a GPR |
| `geneConcLB`, `geneConcUB` | dictionaries gene id → lower/upper bound of the gene concentration (mmol/gDW), e.g. to simulate over-expression or knock-down |
| `verboseMode` | 0: silent, 1: status messages, 2: more details, 3: every reaction and its GPR |

### Returned dictionary

| Key | Content |
|---|---|
| `obj` | objective value (growth rate, 1/h) |
| `status` | solver status (e.g. `'optimal'`) |
| `rxnFlx` | dictionary reaction id → net flux (mmol/gDW/h) |
| `geneConc` | dictionary gene id → enzyme concentration (mmol/gDW) |
| `rxnkcatRow` | DataFrame with one row per reaction, direction (`dirxn`) and enzyme variable, with the kcat value used (`val`, 1/h) |
| `model` | the copy of the model that contains the enzyme constraints |
| `prob` | the optlang problem (`model.solver`) |

Names in the problem: variables `g_<gene>` (gene concentration), `cpx_<rxn>` / `cpx_<rxn>_<k>`
(complexes), `gr_<gene>_<rxn>` (share of a multi-functional enzyme, `ccFBA_mr`); constraints
`<rxn>_kcat`, `<rxn>_b_kcat`, `<rxn>`, `<rxn>_b`, `<rxn>_<gene>`, `Aux_<n>`, `mfe_<gene>` and `MW`.

### ccFBA as a constraint matrix (`getccFBA_mat`)

`getccFBA_mat` builds the ccFBA problem as a matrix, as in R: multi-functional enzymes as in `ccFBA_mr`, but each
reaction uses only its cheapest isoenzyme (`getGpr1iso`). The number of subunits of the genes of a complex can be
given with `cpx_stoich`: `v <= kcat * cpx`, `g_gene >= stoich * cpx`, and a one-gene enzyme with `n` subunits gives
`v <= kcat * g / n`.

```python
from getccFBA_mat import getccFBA_mat, solve_ccFBA_mat

mat = getccFBA_mat(model, Kcat=kcat, MW=mw, RHS=0.27, medval=3600*22.6)   # optional: cpx_stoich=DataFrame('genes','stoich')
Res = solve_ccFBA_mat(mat)
print(Res['obj'], Res['status'], Res['rxnFlx']['EX_glc__D_e'])

# reuse the matrix with other bounds (like simulate_EColi in R): change the irreversible model and solve again
mod2 = mat['model']
mod2.uppbnd[mod2.react_id.index('EX_glc__D_e_b')] = 5   # glucose uptake <= 5
Res = solve_ccFBA_mat(mat)
```

- The matrix is built on the irreversible model `mat['model']` (`mod2irrev(model)`, can also be given as `mod2`): a
  reaction with `lower_bound < 0` is split in `<id>_f` and `<id>_b`. Set the bounds of the cobra model **before**
  building the matrix, e.g. open the uptake (`lower_bound < 0`) of every carbon source you want to use later, or
  give these reactions in `mod2irrev(model, rev=[...])` (like setting `react_rev` to `TRUE` in R).
- `solve_ccFBA_mat(mat, mod2=..., prob=Res['prob'])` reuses the problem of an earlier call and only updates the flux
  bounds (faster when solving many conditions).
- Parameters as in `ccFBA_mr`, and: `cpx_stoich` (DataFrame with columns `genes`, sorted comma separated gene ids,
  and `stoich`, comma separated subunit numbers, e.g. `R/sybilccFBA/inst/extdata/cpx_stoich_me.csv`); `objVal`: when
  `RHS` is `None`, the problem is to find the minimum budget that gives this objective value (default: FBA optimum);
  `C_mu_coef`: budget as a function of growth, `C = RHS + C_mu_coef * mu`. `geneConcLB`/`geneConcUB` can also be
  lists in the order of `model.genes`. `selected_rxns` defaults to all reactions with a GPR (including `s0001`, which
  has no molecular weight and so no cost).
- Returned dictionary (same names as in R, but all row/column indices are **0-based**): `LHS` (scipy sparse),
  `rlb`, `rub`, `rtype`, `rnames` (rows; `rtype` as in sybil: `E` =, `U` <= `rub`, `L` >= `rlb`), `clb`, `cub`, `cnames`
  (columns), `obj_cf`, `obj_dir`, `rxnkcatRow`, `minGiRow`, `mfeRow`, `geneCol`, `mfeCol`, `gpr1iso`, `rgst`,
  `rxnExprCol`, `rxnGeneCol`, `kcat_ir`, `kcat_src_ir`, `geneMW`, `model`, `cc_ind`, `mod_objc_ind`.
- Names in the problem: columns `R_<rxn>`, `g1_<gene>_1` (gene of one reaction), `gmr_<gene>_<n>` (gene of `n`
  reactions), `gr_<gene>_<rxn>`, `cpx_<rxn>`; rows `S_<met>`, `rkcat_<rxn>`, `rbkcat_<rxn>_b`,
  `min_<gene>_<rxn>_<stoich>`, `mfe_<gene>`, `cc` (budget) or `ObjC_` (objective, when `RHS` is `None`).
- `solve_ccFBA_mat(mat)` returns `obj`, `status`, `x` (all columns), `rxnFlx` (net flux of each reaction of the cobra
  model), `flx_ir`, `geneConc`, `prob`.

Differences with the R version: the constraint on the objective when `RHS` is `None` is written to a new row of
type `L` (R writes it to a wrong row with type `U`, so it has no effect); reactions that are not in `selected_rxns`
get no enzyme column (`exprCol` is NA), and Kcat does not need a `src` column.

### Growth on several carbon sources (`simulate_EColi`)

`simulate_EColi(model, mod2, kcat, mw, budget_C, CSList, atpmz=True, trns_rxns=None, cpx_stoich=None)` builds the
matrix once (`getccFBA_mat`) and, for each carbon source in `CSList` (uptake reactions of `mod2`, e.g.
`'EX_glc__D(e)_b'`), closes the uptake of the other carbon sources and solves. `trns_rxns` (booleans in the order of
`model.reactions`, or reaction ids) are excluded from the budget. It returns `gr_data` (growth rates, column `flx`),
`uptake`, `all_flx`, `all_flx_MC`, `rg_MC` (protein cost of reactions, `rxnMC`, and genes, `rgMC`, sorted) as in R, and
`mrMat` (the matrix). `mod2` is not modified. Difference with R: the enzyme concentration of a reaction (`enzConc`) is
`x[exprCol] / r1gstoich` as in the kcat constraint; R divides by the stoichiometry of each gene, which gives `Inf`
for the genes of isoenzymes that were not chosen and a wrong cost for complexes with stoichiometry. With this
definition the costs of the reactions (`rxnMC`) add up to the budget.

### Over expression of a gene (`simGeneOverExpr`)

`simGeneOverExpr(gx, y_pct, model, mod2, Kcat, MW=mw, RHS=0.27, medval=..., FUN=cfba_moment)` makes gene `gx` use at
least the fraction `y_pct` of the budget: the wild type is solved with `FUN` (`cfba_moment`, `ccFBA_mr` or
`getccFBA_sol`); if the gene uses less than `y_pct`, the problem is solved again with the lower bound
`RHS*y_pct/MW` of its concentration. Then the fluxes of the reactions with a measured gene and of the objective are
fixed and the total flux is minimized (MTF). Returns `Sol` (MTF), `OExSol`, `WTSol`, `OverExprEff`, `objVal` and
`gxMC_WT`, `gxMC_OEx` (MC of the gene). When the gene has no concentration variable (with `getccFBA_sol`, a gene that is
not in the chosen isoenzyme of any reaction; R stops with an error), the protein cannot catalyze any reaction, so the
budget is reduced by `y_pct` instead. Bounds are taken from the cobra model (`cfba_moment`, `ccFBA_mr`) or from `mod2`
(`getccFBA_sol`).

## Data and models

| File | Content |
|---|---|
| `models/iAF1260.xml` | *E. coli* iAF1260 (SBML) |
| `models/iMM904.xml` | *S. cerevisiae* iMM904 (SBML) |
| `data/ecoli_kcat.csv` | *E. coli* kcat values used in MOMENT (same as R dataset `kcat`) |
| `data/ecoli_mw.csv` | *E. coli* molecular weights (similar to R dataset `mw`: 288 more genes, 7 values differ by up to 2%) |
| `R/Tin/` | `simGeneOverExpr.R`, its test `test_simGeneOverExpr.R` and the kcat file `kcat2` used by the test |
| `data/yeast_imm_kcat.csv`, `data/yst_mw.csv` | *S. cerevisiae* kcat values and molecular weights (from R datasets `yeast_kcat`, `yst_mw`) |

The R package contains more data: models iJO1366 and iML1515, molecular weights for iML1515, and
protein complex stoichiometry (`inst/extdata/cpx_stoich_me.csv`). The `.RData` files can be read with
`sybil2cobra.read_RData`, e.g. `read_RData('R/sybilccFBA/data/iML1515.RData')['iML1515']`.

## Test notebooks

Both notebooks simulate growth of iAF1260 on 24 carbon sources (oxygen uptake open, ATPM = 0,
`RHS=0.27`, `medval=3600*22.6`). They compare the predicted growth rates with the measured ones
(Adadi et al. 2012), report the uptake of each carbon source and draw a scatter plot.

| Notebook | Method | Pearson R (predicted vs measured growth) |
|---|---|---|
| [`test_moment.ipynb`](test_moment.ipynb) | `cfba_moment` | 0.412 |
| [`test_ccFBA_mr.ipynb`](test_ccFBA_mr.ipynb) | `ccFBA_mr` | 0.398 |
| [`test_ccFBA_mat.ipynb`](test_ccFBA_mat.ipynb) | `getccFBA_mat` | 0.394 |

[`test_simulate_EColi.ipynb`](test_simulate_EColi.ipynb) runs the example of the R help page of `simulate_EColi`:
iML1515 (from `iML1515.RData`), kcat_max of `allKcats_upd34_dd_h.csv`, complex stoichiometry, C = 0.27, `tex`
transport reactions excluded, ATPM not set to zero. Spearman correlation 0.37, Pearson R 0.47.

[`test_simGeneOverExpr.ipynb`](test_simGeneOverExpr.ipynb) runs `R/Tin/test_simGeneOverExpr.R`: over expression of
b2323 with `cfba_moment` (5–50% of C) and of b0929 with the three ccFBA functions, on iAF1260 (from `iAF1260.RData`).

With the same `RHS`, `ccFBA_mr` predicts lower growth rates (0.20–0.29 1/h vs 0.32–0.57 1/h for
MOMENT), because the capacity of multi-functional enzymes is shared between their reactions. The
notebooks use Windows paths (`models\\iAF1260.xml`); on Linux/macOS, use `/`.

## Known issues

- `ccFBA_mr` with `selected_rxns`: the problem is solved, but collecting `geneConc` fails with a
  `KeyError` for genes that catalyze none of the selected reactions (`ccFBA_mr.py`, line 346).
- `mod2irrev.py` (indentation error) .

## License

GPL-3 (same as sybilccFBA), see [LICENSE](LICENSE).
