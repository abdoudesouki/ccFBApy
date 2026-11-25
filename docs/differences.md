# Differences from R and known issues

A summary of the notes on the function pages.

## General

| R (sybilccFBA) | Python |
|---|---|
| sybil `modelorg` | `cobra.Model` |
| `mod2 = mod2irrev(model)` | [`mod2irrev(model, rev=None)`](reference/mod2irrev.md) (`IrrevModel`). `cfba_moment` and `ccFBA_mr` have no `mod2` and use the bounds of the cobra model. |
| `uppbnd(mod2)[...] <- x` | `mod2.uppbnd[mod2.react_id.index(...)] = x` |
| `sysBiolAlg(..., algorithm="ccFBA")` + `optimizeProb()` | [`solve_ccFBA_mat(mat)`](reference/getccFBA_mat.md#solve_ccfba_mat) |
| `solver = "glpkAPI"` / `"cplexAPI"` | argument accepted but not used; set the solver in cobra/optlang |
| `verboseMode = 3` writes the LP file | not available; use `prob.to_lp()` |
| results: sybil `optsol` (`sol$fluxes`), `geneCol` | `dict` with `obj`, `status`, `rxnFlx`, `geneConc`, `prob` |
| 1-based indices | 0-based indices |
| `getRevFlux()` | `rxnFlx` in the results, `IrrevModel.net_flux()` |
| `data(iAF1260)`, `data(kcat)`, ... | `datasets.load_model('iAF1260')`, `datasets.load_dataset('kcat')`, ...: the data sets as SBML and CSV files (no `.RData`) |

## Per function

- [`cfba_moment`](reference/cfba_moment.md): no `mod2`, `runFVA`, `fvaRxn`. `objVal` is
  ignored and `RHS=None` gives no budget constraint, not the minimum budget.
- [`cfba_moment_mr` → `ccFBA_mr`](reference/cfba_moment_mr.md): also no `C_mu_coef`. The
  multi-functional enzyme constraints are `sum <= g` (R: `=`). Returns fewer results (no
  `rxnMC`, `rxnGeneMC`, `actual_C`, ...).
- [`getGpr1iso`](reference/getGpr1iso.md): arguments `(model, MW, cpx_stoich=None)`; R:
  `(model, cpx_stoich, MW)`.
- [`getccFBA_mat`](reference/getccFBA_mat.md): with `RHS=None`, the constraint on the objective
  has an effect (in R it is written to a wrong row). Reactions that are not selected get no enzyme
  column.
- [`getccFBA_sol`](reference/getccFBA_sol.md), [`simulate_EColi`](reference/simulate_EColi.md):
  the enzyme concentration of a reaction is `x[exprCol] / r1gstoich`. R divides by the
  stoichiometry of each gene.
- [`simGeneOverExpr`](reference/simGeneOverExpr.md): when the gene has no concentration
  variable, the budget is reduced by `y_pct` (R stops with an error).
- [`calc_MW`](reference/calc_MW.md): returns numbers (R: text). A `.faa` file with a single protein
  is handled. `nchrm` defaults to the number of files.
- [`readmodel`](reference/readmodel.md): the size of the problem is taken from the inputs (fixed in R).
  The budget `C` can be changed. The LP file is only written when `lp_fname` is given. Returns a
  `dict` (R: the vector of column values, `res['x']`). `read_moment_files` reads `mat.csv` without
  a header, so the first entry is not lost as with `read.table(header=T)`.
- [`addGlcTrns`](reference/addGlcTrns.md): completed. The R function cannot run because the
  transporter columns are never added. Python adds them (`trns_<t>`, bound `trns_ub`) and has options
  for the reaction, `Si`, the transporters and existing variables (`trnsCol`).
- [`cfba_moment_pw`](reference/cfba_moment_pw.md): built as a matrix and solved with `solve_ccFBA_mat`.
  `s0001` is removed from all OR rules. Nested rules are expanded. It works without `RHS`/`objVal`
  (no budget). The option `pairwise=False` writes one row per rule.
- All the functions of the R package are ported, except the class `sysBiolAlg_ccFBA` (replaced by
  `solve_ccFBA_mat`) and `getRevFlux` (not needed).

## Known issues

- `ccFBA_mr` with `selected_rxns`: the problem is solved, but collecting `geneConc` fails with a
  `KeyError` for the genes that are not in any selected reaction.
- `cfba_moment` and `ccFBA_mr` read the kcat value from the **second column** of `Kcat` and only
  match `rxn_id` with the `R_` prefix. A table without the prefix (e.g. from
  `allKcats_upd34_dd_h.csv`) gives `medval` for every reaction, with no warning. See
  [Kcat](conventions.md#kinetic-data-kcat).
- `cfba_moment` and `ccFBA_mr`: the pseudo gene `s0001` gets a concentration variable with no
  cost, so reactions with rules like `b0875 or s0001` are not limited in practice. R removes
  `s0001` from these rules (see [conventions](conventions.md#spontaneous-reactions-s0001)).
- `simulate_EColi(atpmz=True)` only finds a reaction with the id `ATPM` (not `R_ATPM` of the
  `iAF1260`/`iMM904` data sets).
- `test_moment.ipynb` and `test_ccFBA_mr.ipynb` use Windows paths with lower case file names
  (`data\ecoli_kcat.csv`). The files are `data/Ecoli_kcat.csv` and `data/Ecoli_MW.csv`, which
  matters on case sensitive file systems (Linux).

## Errors in the R help pages

These are not repeated in the pages here:

- The unit of `MW` is given as g/mol; the values are in g/mmol (kDa).
- `getccFBA_sol.Rd` has the title and description of `getccFBA_mat`.
- `yst_kcat.Rd` says "a data frame"; the data set is a numeric matrix.
- `sybilccFBA-package.Rd` gives version 0.0.1 (2013); the package is version 3.0.2.
- The example of `addGlcTrns.Rd` is the body of the function, and the R function cannot run (see
  [addGlcTrns](reference/addGlcTrns.md)).
- The example of `cfba_moment_pw.Rd` calls `cfba_moment`, not `cfba_moment_pw`.
