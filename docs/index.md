# ccFBApy documentation

Cost constrained flux balance analysis (ccFBA) in Python. This is a port of the R package
**sybilccFBA** (version 3.0.2, by Abdelmoneim Amer Desouki). It is built on
[COBRApy](https://github.com/opencobra/cobrapy) and [optlang](https://github.com/opencobra/optlang)
instead of [sybil](https://CRAN.R-project.org/package=sybil).

These pages are based on the R help pages in `R/sybilccFBA/man/*.Rd`. Each R help page has a page here that
describes the Python function, or says what replaces it. Each page
also gives the R usage for comparison.

## What the package does

Flux balance analysis (FBA) finds the fluxes of a metabolic model at steady state (Orth, Thiele &
Palsson 2010). Cost constrained FBA also limits the fluxes by the amount of enzyme the cell can afford.
It needs three inputs besides the model:

- kinetic data: the turnover number `kcat` of the enzyme of each reaction,
- the molecular weight (MW) of each protein,
- the gene-protein-reaction (GPR) rules of the model.

Two constraints are added to FBA:

```
v_i <= kcat_i * g_i        # flux <= turnover number x enzyme concentration
sum_g (MW_g * g_g) <= C    # total enzyme mass <= budget C (0.27 g/gDW for E. coli)
```

The package implements:

1. **MOMENT** (MetabOlic Modeling with ENzyme kineTics), Adadi et al. (2012):
   [`cfba_moment`](reference/cfba_moment.md), and [`cfba_moment_pw`](reference/cfba_moment_pw.md) written like the
   original MATLAB implementation.
2. **ccFBA**, an improvement of MOMENT (Desouki 2016):
    - multi-functional enzymes: an enzyme that catalyzes several reactions must share its
      concentration between them ([`ccFBA_mr`](reference/cfba_moment_mr.md));
    - only the cheapest isoenzyme of each reaction is used
      ([`getGpr1iso`](reference/getGpr1iso.md));
    - optional protein complex stoichiometry
      ([`getccFBA_mat`](reference/getccFBA_mat.md), [`getccFBA_sol`](reference/getccFBA_sol.md));
    - new in Python: [`ccFBA_iso`](reference/ccFBA_iso.md) keeps all the isoenzymes and allocates each gene
      to the complexes that use it, instead of choosing one isoenzyme with `getGpr1iso`.

## Installation

```
pip install ccFBApy              # the package (imported as ccfbapy), with the data sets of sybilccFBA
pip install "ccFBApy[all]"       # + matplotlib for scatterPlot
```

From a copy of the repository: `pip install .` (or `pip install -e ".[all,test]"` for development; run the
tests with `pytest`). Without installing, run Python from the root folder of the repository, where the
folder `ccfbapy` is.

The examples of this manual read the SBML models of `models/` and the tables of `data/`, or the data sets
shipped with the package ([`ccfbapy.datasets`](data.md#data-sets-in-the-package): SBML and CSV files). Run
them from the root folder of the repository.

## Quick start

```python
import cobra
import pandas as pd
from ccfbapy import cfba_moment

model = cobra.io.read_sbml_model('models/iAF1260.xml')   # run from the root folder of the repository
mw = pd.read_csv('data/Ecoli_MW.csv')      # gene, molecular weight (g/mmol)
kcat = pd.read_csv('data/Ecoli_kcat.csv')  # rxn_id, val (1/s), dirxn

res = cfba_moment(model, kcat, MW=mw, RHS=0.27, medval=3600*22.6, verboseMode=1)
print(res['status'], res['obj'])           # optimal 0.4552...
print(res['rxnFlx']['EX_glc__D_e'])        # net flux of glucose exchange
```

See [Inputs, units and conventions](conventions.md) for the formats of `Kcat`, `MW` and `cpx_stoich`,
and [Differences from R and known issues](differences.md) before you port R code.

## R help pages and their Python equivalents

| R help page (`man/*.Rd`) | Python | Status |
|---|---|---|
| `sybilccFBA-package.Rd` | this page | |
| [`cfba_moment.Rd`](reference/cfba_moment.md) | `ccfbapy.cfba_moment()` | ported |
| [`cfba_moment_mr.Rd`](reference/cfba_moment_mr.md) | `ccfbapy.ccFBA_mr()` | ported (no FVA) |
| [`cfba_moment_pw.Rd`](reference/cfba_moment_pw.md) | `ccfbapy.cfba_moment_pw()` | ported |
| [`getGpr1iso.Rd`](reference/getGpr1iso.md) | `ccfbapy.getGpr1iso()` | ported |
| [`getccFBA_mat.Rd`](reference/getccFBA_mat.md) | `ccfbapy.getccFBA_mat()`, `solve_ccFBA_mat()` | ported |
| [`getccFBA_sol.Rd`](reference/getccFBA_sol.md) | `ccfbapy.getccFBA_sol()` | ported |
| [`getRevFlux.Rd`](reference/getRevFlux.md) | `IrrevModel.net_flux()`, `rxnFlx` in the results | not needed |
| [`simulate_EColi.Rd`](reference/simulate_EColi.md) | `ccfbapy.simulate_EColi()` | ported |
| [`sysBiolAlg_ccFBA-class.Rd`](reference/sysBiolAlg_ccFBA.md) | `solve_ccFBA_mat()` | replaced |
| [`initialize-methods.Rd`](reference/sysBiolAlg_ccFBA.md#initialize-method) | `solve_ccFBA_mat()` | replaced |
| [`calc_MW.Rd`](reference/calc_MW.md) | `ccfbapy.calc_MW()` | ported |
| [`readmodel.Rd`](reference/readmodel.md) | `ccfbapy.readmodel()`, `read_moment_files()` | ported |
| [`addGlcTrns.Rd`](reference/addGlcTrns.md) | `ccfbapy.addGlcTrns()` | ported (completed) |
| [`iAF1260.Rd`](data.md#iaf1260), [`iJO1366.Rd`](data.md#ijo1366), [`iML1515.Rd`](data.md#iml1515), [`iMM904.Rd`](data.md#imm904) | `datasets.load_model()` (SBML), `models/*.xml` | available |
| [`kcat.Rd`](data.md#kcat), [`mw.Rd`](data.md#mw), [`mw_iML1515.Rd`](data.md#mw_iml1515) | `datasets.load_dataset()` (CSV), `data/*.csv` | available |
| [`yst_kcat.Rd`](data.md#yst_kcat), [`yst_mw.Rd`](data.md#yst_mw) | `datasets.load_dataset()` (CSV), `data/yeast_*.csv` | available |

Python functions that have no R help page:

| Python | Page |
|---|---|
| `ccfbapy.mod2irrev()`, class `IrrevModel` (sybil's `mod2irrev`) | [mod2irrev](reference/mod2irrev.md) |
| `ccfbapy.solve_ccFBA_mat()` | [getccFBA_mat](reference/getccFBA_mat.md#solve_ccfba_mat) |
| `ccfbapy.flx_MC_tables()` | [getccFBA_sol](reference/getccFBA_sol.md#flx_mc_tables) |
| `ccfbapy.ccFBA_iso()`, `getccFBA_iso_mat()`: isoenzymes without `getGpr1iso` (gene allocation) | [ccFBA_iso](reference/ccFBA_iso.md) |
| `ccfbapy.simGeneOverExpr()` (R source `R/OverExpr/simGeneOverExpr.R`) | [simGeneOverExpr](reference/simGeneOverExpr.md) |
| `ccfbapy.datasets`: the data sets of sybilccFBA as SBML and CSV | [Data sets](data.md#data-sets-in-the-package) |
| `ccfbapy.scatterPlot()`, `ccfbapy.remove_enzyme_layer()` | [Utilities](reference/utilities.md) |

## Package information

| | |
|---|---|
| R package | sybilccFBA 3.0.2 (2022-08-13) |
| Title | Cost Constrained Flux Balance Analysis (ccFBA): MetabOlic Modeling with ENzyme kineTics (MOMENT) |
| Author | Abdelmoneim Amer Desouki |
| License | GPL-3 |
| Python package | `ccFBApy` on PyPI, imported as `ccfbapy` |
| Python dependencies | cobra, optlang, pandas, numpy, scipy; optional matplotlib |

The R help page `sybilccFBA-package.Rd` still gives version 0.0.1 (2013-06-03); `DESCRIPTION` and
`NEWS` give 3.0.2.

## References

- Adadi R, Volkmer B, Milo R, Heinemann M, Shlomi T (2012). Prediction of microbial growth rate
  versus biomass yield by a metabolic network with kinetic parameters. *PLoS Comput Biol* 8(7):
  e1002575. [doi:10.1371/journal.pcbi.1002575](https://doi.org/10.1371/journal.pcbi.1002575)
- Desouki AA (2016). *Algorithms for improving the predictive power of flux balance analysis.*
  PhD dissertation.
- Gelius-Dietrich G, Desouki AA, Fritzemeier CJ, Lercher MJ (2013). sybil – Efficient
  constraint-based modelling in R. *BMC Systems Biology* 7: 125.
- Orth JD, Thiele I, Palsson BØ (2010). What is flux balance analysis? *Nat Biotechnol* 28:
  245–248.
