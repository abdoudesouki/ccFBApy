# ccFBApy

Cost constrained flux balance analysis (ccFBA) in Python, built on [COBRApy](https://github.com/opencobra/cobrapy)
and [optlang](https://github.com/opencobra/optlang). ccFBApy is the Python port of the R package
[sybilccFBA](https://CRAN.R-project.org/package=sybilccFBA) (version 3.0.2) by Abdelmoneim Amer Desouki.

Flux balance analysis (FBA) finds the fluxes of a metabolic model at steady state. Cost constrained FBA also
limits the fluxes by the amount of enzyme the cell can afford, using the turnover numbers (kcat) of the
enzymes, the molecular weights of the proteins and the gene-protein-reaction (GPR) rules of the model:

```
v_i <= kcat_i * g_i          flux <= kcat x enzyme concentration
sum_g MW_g * g_g <= C        total enzyme mass <= budget C (0.27 g/gDW for E. coli)
```

## Installation

```
pip install ccFBApy              # the methods and the data sets of sybilccFBA
pip install "ccFBApy[all]"       # + matplotlib for scatterPlot
```

The package is imported as `ccfbapy`. The LP solver is the one of cobra/optlang: GLPK (`swiglpk`) by default;
CPLEX, Gurobi, ... can be used when they are installed.

## Methods

| Function | Method |
|---|---|
| `cfba_moment` | MOMENT: MetabOlic Modeling with ENzyme kineTics (Adadi et al. 2012) |
| `cfba_moment_pw` | MOMENT on the irreversible model with pairwise OR, like the original MATLAB implementation |
| `ccFBA_mr` (R: `cfba_moment_mr`) | MOMENT with multi-functional enzymes: an enzyme shares its concentration between its reactions |
| `getccFBA_mat`, `solve_ccFBA_mat`, `getccFBA_sol` | ccFBA as a constraint matrix: multi-functional enzymes, the cheapest isoenzyme of each reaction (`getGpr1iso`), protein complex stoichiometry; build once, solve for many conditions |
| `ccFBA_iso`, `getccFBA_iso_mat` | ccFBA with gene allocation: all isoenzymes are kept and each gene is shared between the complexes that use it, so shared subunits are not counted twice and knocked out or limited genes are handled (new in ccFBApy) |
| `simulate_EColi` | growth on a list of carbon sources, with the protein cost of each reaction |
| `simGeneOverExpr` | over expression of one gene |
| `getGpr1iso`, `mod2irrev`, `calc_MW`, `readmodel`, `addGlcTrns` | cheapest isoenzyme, irreversible model, molecular weights from a genome (.faa), MOMENT model exported from MATLAB, glucose transport constraint |
| `datasets` | data sets of sybilccFBA as SBML and CSV files: models iAF1260, iJO1366, iML1515 (*E. coli*) and iMM904 (yeast), kcat values, molecular weights, complex stoichiometry, example files; also the BiGG versions of the models (`load_model(name, source='bigg')`) |

## Example

MOMENT and ccFBA for *E. coli* iAF1260 growing on glucose, with the data sets shipped with the package:

```python
import ccfbapy
from ccfbapy import datasets

model = datasets.load_model('iAF1260')          # cobra.Model (SBML), ids as in sybil (R_...)
kcat = datasets.load_dataset('kcat')            # kcat values (1/s): rxn_id, val, dirxn
mw = datasets.load_dataset('mw')                # molecular weights (g/mmol): gene, mw

model.reactions.get_by_id('R_EX_glc_e_').lower_bound = -1000   # glucose
model.reactions.get_by_id('R_EX_o2_e_').lower_bound = -1000    # oxygen
model.reactions.get_by_id('R_ATPM').lower_bound = 0

res = ccfbapy.cfba_moment(model, kcat, MW=mw, RHS=0.27, medval=3600*22.6, verboseMode=0)
print('MOMENT growth rate:', res['obj'])

res = ccfbapy.ccFBA_iso(model, Kcat=kcat, MW=mw, RHS=0.27, medval=3600*22.6, verboseMode=0)
print('ccFBA growth rate:', res['obj'])
print(res['rxnMC'].head())                      # protein cost of the reactions (g/gDW)

# knock out a gene: the other isoenzymes take over
res_ko = ccfbapy.ccFBA_iso(model, mat=res['mat'], geneConcUB={'b0118': 0}, verboseMode=0)
print('acnB knock out:', res_ko['obj'])
```

Inputs: `Kcat` is a DataFrame with the columns `rxn_id`, `val` (kcat in 1/s) and `dirxn` (1 forward, -1
reverse); `MW` has the gene in the first column and the molecular weight in g/mmol in the second; `RHS` is the
budget C in g protein/gDW; `medval` is the kcat (in 1/h) of the reactions without a value.

## Documentation

The reference manual, with one section per help page of sybilccFBA and the Python examples, is in the
repository of the project (`docs/`, `reference_manual.html`, `reference_manual.ipynb`).

## References

- Adadi R, Volkmer B, Milo R, Heinemann M, Shlomi T (2012). Prediction of microbial growth rate versus biomass
  yield by a metabolic network with kinetic parameters. PLoS Comput Biol 8(7): e1002575.
  [doi:10.1371/journal.pcbi.1002575](https://doi.org/10.1371/journal.pcbi.1002575)
- Desouki AA (2016). Algorithms for improving the predictive power of flux balance analysis. PhD dissertation.
- Gelius-Dietrich G, Desouki AA, Fritzemeier CJ, Lercher MJ (2013). sybil – Efficient constraint-based modelling
  in R. BMC Systems Biology 7: 125.

## License

GPL-3, as sybilccFBA.
