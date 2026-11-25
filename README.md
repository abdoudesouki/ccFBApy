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

## Documentation

The folder [`docs`](docs) has one page for each R help page in
[`R/sybilccFBA/man`](R/sybilccFBA/man). Each page gives the Python usage, arguments, details,
value, differences from R and examples (the examples have been run). Start with
[`docs/index.md`](docs/index.md). The same manual, with the examples as code cells you can run (and their outputs),
is in the notebook [`reference_manual.ipynb`](reference_manual.ipynb), and in one HTML file,
[`reference_manual.html`](reference_manual.html) (open it in a browser). Both have a table of contents, and each
section links back to it. They are generated from `docs/`: after you change a page, run
`python tools/make_reference_notebook.py --execute --html`. To build the multi-page HTML site (output in `site/`):

```
pip install mkdocs mkdocs-material
mkdocs build        # or: mkdocs serve
```

## Status of the Python implementation

| R (sybilccFBA) | Python | Status |
|---|---|---|
| `cfba_moment` | `ccfbapy/cfba_moment.py`: `cfba_moment()` | Working, tested on iAF1260 ([`test_moment.ipynb`](test_moment.ipynb)) |
| `cfba_moment_mr` | `ccfbapy/ccFBA_mr.py`: `ccFBA_mr()` | Working, tested on iAF1260 ([`test_ccFBA_mr.ipynb`](test_ccFBA_mr.ipynb)). See known issues for `selected_rxns` |
| `getGpr1iso` | `ccfbapy/gpr1iso.py`: `getGpr1iso()` | Working, including complex stoichiometry (`cpx_stoich`) |
| `getccFBA_mat` | `ccfbapy/getccFBA_mat.py`: `getccFBA_mat()` | Working, tested on iAF1260 ([`test_ccFBA_mat.ipynb`](test_ccFBA_mat.ipynb)). The same file has `mod2irrev()` (irreversible model) and `solve_ccFBA_mat()` (solves the matrix, replaces `sysBiolAlg`/`optimizeProb`) |
| `getRevFlux` | not needed | `rxnFlx` in the results already holds net fluxes (forward − reverse) |
| `simulate_EColi` | `ccfbapy/simulate_EColi.py`: `simulate_EColi()` | Working, example of the R help page in [`test_simulate_EColi.ipynb`](test_simulate_EColi.ipynb) (iML1515) |
| `getccFBA_sol` | `ccfbapy/getccFBA_sol.py`: `getccFBA_sol()` | Working: builds and solves the matrix; also returns `obj`, `status`, `rxnFlx`, `geneConc` (dict) like the other python functions |
| `simGeneOverExpr` (`R/OverExpr/simGeneOverExpr.R`) | `ccfbapy/simGeneOverExpr.py`: `simGeneOverExpr()` | Working, test in [`test_simGeneOverExpr.ipynb`](test_simGeneOverExpr.ipynb) (python version of `R/OverExpr/test_simGeneOverExpr.R`) |
| `calc_MW` | `ccfbapy/calc_MW.py`: `calc_MW()` | Working: reproduces the data set `mw_iML1515` from `NC_000913.faa` (mw in g/mol, divide by 1000 for g/mmol) |
| `readmodel` | `ccfbapy/readmodel.py`: `readmodel()`, `read_moment_files()` | Working: solves the MOMENT model exported from MATLAB (`inst/extdata`), growth 0.508 with C = 0.27 |
| `addGlcTrns` | `ccfbapy/addGlcTrns.py`: `addGlcTrns()` | Completed (the R function cannot run: the transporter columns are not added); works on cobra, `cfba_moment` and `getccFBA_mat` problems |
| `cfba_moment_pw` | `ccfbapy/cfba_moment_pw.py`: `cfba_moment_pw()` | Working: MOMENT on the irreversible model with pairwise OR; iAF1260 growth 0.5083, as the MATLAB model (`readmodel`) |
| class `sysBiolAlg_ccFBA` | `ccfbapy/getccFBA_mat.py`: `solve_ccFBA_mat()` | Replaced |
| (new) | `ccfbapy/ccFBA_iso.py`: `ccFBA_iso()`, `getccFBA_iso_mat()` | New: isoenzymes without `getGpr1iso`. Each gene is allocated to the complexes that use it (`sum_t st*e_rt = g`), so shared subunits are not counted twice, and the LP chooses the isoenzymes (knock outs and gene limits are handled). See [`docs/reference/ccFBA_iso.md`](docs/reference/ccFBA_iso.md) |
| data sets (`R/sybilccFBA/data/*.RData`) | `ccfbapy/datasets.py`: `load_model()`, `load_dataset()` | The data sets of the R package, exported once to SBML (models, ids kept as in sybil, e.g. `EX_glc__D(e)`) and CSV (tables) with `tools/make_package_data.py`; the package does not need `rdata` |
| | `ccfbapy/genFuns.py`: `scatterPlot()` | Helper: measured vs predicted scatter plot with y = x line, linear fit and Pearson correlation |

Options of the R functions that are **not** available in Python yet: `runFVA`, `fvaRxn`, and in
`cfba_moment`/`ccFBA_mr` also `C_mu_coef` and `objVal` (accepted but unused; both work in `getccFBA_mat`). `solver` is also accepted but unused: the problem is solved
with the model's solver. Only GLPK, cobra's default, has been tested.

## Installation

ccFBApy is a Python package (`pyproject.toml`), published on PyPI as `ccFBApy` and imported as `ccfbapy`:

```
pip install ccFBApy              # the methods and the data sets (cobra, optlang, numpy, pandas, scipy)
pip install "ccFBApy[all]"       # + matplotlib for scatterPlot
```

From this repository: `pip install .`, or for development `pip install -e ".[all,test]"` and run the tests with
`pytest` (folder `tests`). Without installing, `import ccfbapy` works from the root folder of the repository.
Building and publishing the package is described in [RELEASING.md](RELEASING.md).

Tested with Python 3.12, cobra 0.32.1, optlang 1.9.1 (GLPK solver). The notebooks were created with Python 3.11.

## Usage

```python
import cobra
import pandas as pd
import ccfbapy

model = cobra.io.read_sbml_model('models/iAF1260.xml')
mw = pd.read_csv('data/Ecoli_MW.csv')      # gene, molecular weight (g/mmol)
kcat = pd.read_csv('data/Ecoli_kcat.csv')  # rxn_id, val (1/s), dirxn (1: forward, -1: reverse)

# MOMENT
Res = ccfbapy.cfba_moment(model, kcat, MW=mw, verboseMode=1, RHS=0.27, medval=3600*22.6)
print(Res['obj'], Res['status'])

# ccFBA with multi-functional enzymes
Res_mr = ccfbapy.ccFBA_mr(model, kcat, MW=mw, verboseMode=1, RHS=0.27, medval=3600*22.6)
```

The same data are in the package: `ccfbapy.datasets.load_kcat()`, `load_mw()`, `load_model('iAF1260')`
(sybil ids), `load_cpx_stoich()`, ... (see `docs/data.md`). The functions
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
from ccfbapy import getccFBA_mat, solve_ccFBA_mat

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
| `data/Ecoli_kcat.csv` | *E. coli* kcat values used in MOMENT (same as R dataset `kcat`) |
| `data/Ecoli_MW.csv` | *E. coli* molecular weights (similar to R dataset `mw`: 288 more genes, 7 values differ by up to 2%) |
| `ccfbapy/data/` | the data shipped with the package (`ccfbapy.datasets`): the R data sets as SBML and CSV (`sybilccFBA/`), the models of `models/` (`bigg/*.xml.gz`), the CSV tables of `data/` and the files of `R/sybilccFBA/inst/extdata`. Made with `tools/make_package_data.py` |
| `R/OverExpr/` | `simGeneOverExpr.R`, its test `test_simGeneOverExpr.R` and the kcat file `kcat2` used by the test |
| `data/yeast_imm_kcat.csv`, `data/yeast_mw.csv` | *S. cerevisiae* kcat values and molecular weights (from R datasets `yst_kcat`, `yst_mw`) |

The R package contains more data: models iJO1366 and iML1515, molecular weights for iML1515, and
protein complex stoichiometry (`inst/extdata/cpx_stoich_me.csv`). Its data sets are in the package as SBML and
CSV files: `ccfbapy.datasets.load_model('iML1515')`, `load_dataset('mw_iML1515')`, ... (the `.RData` files are only read
by the development tool `tools/make_package_data.py`, which needs `rdata`).

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
iML1515 (the R data set, `datasets.load_model('iML1515')`), kcat_max of `allKcats_upd34_dd_h.csv`, complex stoichiometry, C = 0.27, `tex`
transport reactions excluded, ATPM not set to zero. Spearman correlation 0.37, Pearson R 0.47.

[`test_simGeneOverExpr.ipynb`](test_simGeneOverExpr.ipynb) runs `R/OverExpr/test_simGeneOverExpr.R`: over expression of
b2323 with `cfba_moment` (5–50% of C) and of b0929 with the three ccFBA functions, on iAF1260 (the R data set, `datasets.load_model('iAF1260')`).

With the same `RHS`, `ccFBA_mr` predicts lower growth rates (0.20–0.29 1/h vs 0.32–0.57 1/h for
MOMENT), because the capacity of multi-functional enzymes is shared between their reactions. The
notebooks use Windows paths (`models\\iAF1260.xml`); on Linux/macOS, use `/`.

## Known issues

- `ccFBA_mr` with `selected_rxns`: the problem is solved, but collecting `geneConc` fails with a
  `KeyError` for genes that catalyze none of the selected reactions (`ccfbapy/ccFBA_mr.py`, line 346).
- `test_moment.ipynb` and `test_ccFBA_mr.ipynb` use Windows paths (`models\\iAF1260.xml`, `data\\ecoli_mw.csv`); on
  Linux/macOS use `models/iAF1260.xml` and `data/Ecoli_MW.csv`.

## License

GPL-3 (same as sybilccFBA), see [LICENSE](LICENSE).
