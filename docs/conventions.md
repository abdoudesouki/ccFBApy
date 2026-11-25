# Inputs, units and conventions

The R help pages describe the same arguments (`model`, `mod2`, `Kcat`, `MW`, `RHS`, `medval`,
`cpx_stoich`, ...) on every page. This page describes their Python formats once. The function
pages link here.

## Units

| Quantity | Symbol | Unit |
|---|---|---|
| flux | `v` | mmol/gDW/h |
| turnover number, input | `kcat` | 1/s (converted to 1/h inside the functions: × 3600) |
| turnover number for missing values | `medval` | **1/h** (not converted, e.g. `3600*22.6`) |
| enzyme (gene) concentration | `g` | mmol/gDW |
| molecular weight | `MW` | g/mmol (= kDa), e.g. 89.12 for b0002 |
| budget of metabolic protein | `C`, `RHS` | g protein/gDW (0.27 for *E. coli*) |
| protein cost (molecular crowding) | `MC = MW * g` | g protein/gDW |

The R help pages give the unit of `MW` as g/mol. The values in the data sets (`mw`,
`mw_iML1515`, `yst_mw`, `data/*MW.csv`) are in g/mmol (kDa), which is what the budget
constraint needs (`MW * g` in g/gDW). [`calc_MW`](reference/calc_MW.md) returns g/mol,
so divide its output by 1000: `pd.DataFrame({'gene': geneCnt['Synonym'], 'mw': geneCnt['mw'] / 1000})`.

## Metabolic model (`model`)

A `cobra.Model`. There are two ways to get the models used in the examples:

| Source | Reaction ids | Example |
|---|---|---|
| `datasets.load_model('iAF1260')` (the R data set `iAF1260`) | as in sybil | `R_EX_glc_e_`, `R_ATPM` |
| `datasets.load_model('iML1515')` (the R data set `iML1515`) | as in sybil | `EX_glc__D(e)`, `ATPM` |
| `cobra.io.read_sbml_model('models/iAF1260.xml')` or `datasets.load_model('iAF1260', source='bigg')` | BiGG ids, no `R_` prefix | `EX_glc__D_e`, `ATPM` |

`datasets` is `ccfbapy.datasets`. The models are SBML files; see [Data sets](data.md) and
[Two versions of the models](data.md#two-versions-of-the-models).

A reaction is treated as reversible when its `lower_bound < 0`.

## Irreversible model (`mod2`)

In R, `mod2` is the irreversible model of `model` (sybil's `mod2irrev`). In Python:

- `cfba_moment` and `ccFBA_mr` have **no** `mod2` argument. They work on a copy of the cobra model
  and use its bounds. A reversible reaction gets one kcat constraint on its forward variable and one
  on its reverse variable.
- `getccFBA_mat`, `getccFBA_sol`, `simulate_EColi` and `simGeneOverExpr` (with `getccFBA_sol`) use
  an [`IrrevModel`](reference/mod2irrev.md) made by `mod2irrev(model)`. Reversible reactions are
  split in `<id>_f` and `<id>_b`. The bounds of the fluxes are taken from `mod2.lowbnd` and
  `mod2.uppbnd`, so you can change them there before solving, as in R.

## Kinetic data (`Kcat`)

A `pandas.DataFrame` with one row per reaction and direction:

| Column | Content |
|---|---|
| `rxn_id` | reaction id |
| `val` | kcat in 1/s |
| `dirxn` | `1` forward, `-1` reverse direction |
| `src` | optional: source of the value (used in `getccFBA_mat` outputs, `kcat_src_ir`) |

Reactions without a value for a direction use `medval` (in 1/h). The default is the median of `val`
(after conversion to 1/h).

The functions read the table in different ways:

| Function | Column lookup | Matching of `rxn_id` |
|---|---|---|
| `cfba_moment`, `ccFBA_mr` | `val` must be the **second** column (it is converted by position) | `rxn_id` must be `'R_' + reaction id`. The prefix is added to the reaction id unless it already starts with `R_` |
| `getccFBA_mat`, `getccFBA_sol`, `simulate_EColi` | by name | reaction id with or without the `R_` prefix |

`data/Ecoli_kcat.csv` and the R data set `kcat` already have the `R_` prefix. The kcat table
`allKcats_upd34_dd_h.csv` (column `ijo_id`) does not. To use it with `cfba_moment` or
`ccFBA_mr`, add the prefix: `kcat['rxn_id'] = 'R_' + kcat['rxn_id']`. Otherwise all
reactions silently get `medval`.

## Molecular weights (`MW`)

A `pandas.DataFrame`. The **first column** is the gene id and the **second column** is the molecular
weight in g/mmol. Column names do not matter (`Synonym`/`mw`, `gene`/`mw`, ...). If a gene appears
more than once, the first row is used.

Genes without a molecular weight get no cost: they are left out of the budget constraint.

## Spontaneous reactions (`s0001`)

*E. coli* models use the pseudo gene `s0001` for spontaneous reactions. It has no molecular weight
in `mw`.

- `cfba_moment`, `ccFBA_mr`: reactions whose whole GPR is `s0001` are not constrained. Other
  reactions that contain `s0001` (e.g. `b0875 or s0001`) use a variable `g_s0001` that has no
  cost, so in practice these reactions are free. The R version removes `s0001` from these rules
  (e.g. `b0875 or s0001` becomes `b0875`) for the rule strings it knows. The Python version
  has the same list, but cobra writes the rules without parentheses (`b0875 or s0001`), so
  they do not match and are not changed.
- `cfba_moment_pw`: `s0001` is removed from the OR rules (e.g. `b0875 or s0001` becomes `b0875`),
  as R does.
- `getccFBA_mat`: `getGpr1iso` prefers terms with unknown MW, so `s0001` is chosen for
  `b0875 or s0001` and the reaction has no cost (as in R).
- `ccFBA_iso`: `s0001` is a gene without MW (no cost), unless it is removed with `rm_genes=['s0001']`.
- The R examples add `s0001` with a very small MW (0.001) to the MW table.

## Protein complex stoichiometry (`cpx_stoich`)

Optional `pandas.DataFrame` with at least the columns:

| Column | Content |
|---|---|
| `genes` | comma separated gene ids of the complex, **sorted**, e.g. `b2221,b2222` |
| `stoich` | comma separated number of subunits of each gene, same order, e.g. `2,2` |

Complexes that are not found in the table get one subunit per gene. An example is
`R/sybilccFBA/inst/extdata/cpx_stoich_me.csv`.

## Budget (`RHS`) and objective

- `RHS` is the budget `C`: `sum(MW_g * g_g) <= RHS`.
- `cfba_moment`, `ccFBA_mr`: when `RHS` is `None` no budget constraint is added and `objVal`
  is ignored, so the result is an FBA with the kcat constraints only.
- `getccFBA_mat`: when `RHS` is `None` the problem is to minimize the budget needed to reach the
  objective value `objVal` (default: the FBA optimum).

## Gene concentration bounds (`geneConcLB`, `geneConcUB`)

Lower and upper bounds of the gene concentrations, in mmol/gDW. The defaults are 0 and 1000. You
can give a `dict` `{gene_id: value}` to all functions. `getccFBA_mat` and `getccFBA_sol` also
accept a list in the order of `model.genes`. For example, you can simulate over expression with a
lower bound (see [`simGeneOverExpr`](reference/simGeneOverExpr.md)).

## Solver

The R functions have a `solver` argument (`glpkAPI`, `cplexAPI`, ...). In Python the argument is
accepted but **not used**:

- `cfba_moment` and `ccFBA_mr` use the solver of the model (`model.solver`). To change it:
  `model.solver = 'glpk'`.
- `getccFBA_mat` / `solve_ccFBA_mat` use cobra's default solver
  (`cobra.Configuration().solver`), or the `interface` argument of `solve_ccFBA_mat`.

## `verboseMode`

`0` prints almost nothing, `1` prints status messages, `2` prints more detail, `>2` prints the
processing of each reaction (`cfba_moment`, `ccFBA_mr`). The R option `3` also writes the LP
problem to a file. There is no such option in Python. Use the optlang problem instead, e.g.
`open('prob.lp', 'w').write(res['prob'].to_lp())` (GLPK interface).

## Results

All functions return a `dict`. The ccFBA functions share these keys:

| Key | Content |
|---|---|
| `obj` | objective value (e.g. growth rate, 1/h) |
| `status` | solver status, `'optimal'` when solved |
| `rxnFlx` | `dict` reaction id → net flux (forward − reverse) of the reactions of `model` |
| `geneConc` | `dict` gene id → concentration (mmol/gDW) |
| `prob` | the optlang problem |

The R functions return sybil objects (`sol`, `prob`, `geneCol`). Row and column indices
in the Python results are **0-based**, while they are 1-based in R.
