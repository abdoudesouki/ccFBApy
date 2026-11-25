# cfba_moment_pw: MOMENT with pairwise OR

*R help page: `man/cfba_moment_pw.Rd`. Python: `ccfbapy/cfba_moment_pw.py`, function `cfba_moment_pw()`.*

## Description

Computes fluxes with the MOMENT method (GPR, kcat and molecular weights), written like the original
MATLAB implementation of MOMENT. It works on the irreversible model, and a GPR rule with more than two
isoenzymes (OR terms) is written as a chain of pairwise sums.

## Usage

```python
from ccfbapy import cfba_moment_pw

cfba_moment_pw(model, mod2=None, Kcat=None, MW=None, selected_rxns=None, verboseMode=2, objVal=None,
               RHS=None, solver=None, medval=None, pairwise=True, interface=None)
```

R:

```r
cfba_moment_pw(model, mod2=NULL, Kcat, MW=NULL, selected_rxns=NULL, verboseMode=2, objVal=NULL,
               RHS=NULL, solver=SYBIL_SETTINGS("SOLVER"), medval=NULL)
```

## Arguments

| Argument | Description |
|---|---|
| `model` | `cobra.Model` |
| `mod2` | irreversible model, `mod2irrev(model)` (see [mod2irrev](mod2irrev.md)). Its bounds are used. Default: `mod2irrev(model)`. |
| `Kcat` | **required**. kcat values in 1/s: DataFrame `rxn_id`, `val`, `dirxn`. Reaction ids of `model`, with or without the `R_` prefix. |
| `MW` | molecular weights in g/mmol: first column gene, second column MW. When `None`, there is no budget constraint. |
| `selected_rxns` | reactions to constrain: ids of the **irreversible** model (e.g. `R_PGI_f`). The id of a split reaction of `model` (e.g. `R_PGI`) selects both directions. Default: all reactions of `mod2`. |
| `verboseMode` | `0` silent, `1` status, `>2` each reaction and its rule. Default `2`. |
| `objVal` | when given, the problem is to find the minimum budget that gives this objective value |
| `RHS` | the budget C in g/gDW, 0.27 for *E. coli* |
| `solver` | accepted for compatibility with R, not used (see `interface`) |
| `medval` | kcat in **1/h** for reactions without kcat. Default: median of `Kcat['val']` (in 1/h). |
| `pairwise` | `True` (default): pairwise chain of ORs, as in R. `False`: one row `v <= kcat * (sum of isoenzymes)`, as in [`cfba_moment`](cfba_moment.md); the solution is the same. |
| `interface` | optlang interface used to solve. Default: cobra's default solver. |

## Details

The function works on the irreversible model. The kcat of a backward reaction `<id>_b` is the value of
`<id>` with `dirxn = -1`; the other reactions use `dirxn = 1`. Reactions without a value use
`medval`.

Main steps:

1. Add a variable `g_<gene>` for every gene except `s0001`.
2. For each selected reaction, parse the GPR rule of its original reaction. Reactions without a GPR,
   or with the rule `s0001`, are not constrained. The isoenzyme `s0001` is removed from OR rules
   (e.g. `b0875 or s0001` becomes `b0875`). R does the same for the rule strings of iAF1260 and
   iJO1366 that it knows.
3. Add variables and constraints according to the rule (table below).
4. Add the budget ("solvent") constraint `MW`: `sum(MW_g * g_g) <= RHS`. When `objVal` is given,
   add `ObjC_`: `objective >= objVal`, and minimize `sum(MW_g * g_g)` instead.

Constraints for each type of GPR rule:

| GPR rule | Constraints |
|---|---|
| one gene `A` | `rkcat_<rxn>`: `v <= kcat * g_A` |
| complex `A and B` | `rkcat_<rxn>`: `v <= kcat * cpx_<rxn>`; `min_<gene>_<rxn>`: `cpx <= g` for each gene |
| two isoenzymes `A or B` | `rkcat_<rxn>`: `v <= kcat * (g_A + g_B)` |
| more isoenzymes `A or B or C or D` | `AuxOR_<rxn>_1`: `auxOR_1 <= g_A + g_B`; `AuxOR_<rxn>_2`: `auxOR_2 <= g_C + auxOR_1`; `rkcat_<rxn>`: `v <= kcat * (g_D + auxOR_2)` |

An AND term of an OR rule gets a variable `cpx_<rxn>_<p>` (term `p`, 0-based), with the rows
`min_<gene>_<rxn>_<p>`. The terms are chained in the order of their columns (genes first, in the order of the
model, then the complexes), as in R.

The chain is equivalent to `v <= kcat * (g_A + g_B + g_C + g_D)`, which [`cfba_moment`](cfba_moment.md)
writes in one row. `pairwise=False` writes the single row and gives the same optimum.

## Value

A `dict`:

| Key | Content |
|---|---|
| `obj`, `status` | objective value (or the minimum budget, with `objVal`) and solver status |
| `rxnFlx` | `dict` reaction id of `model` → net flux |
| `flx_ir` | fluxes of the reactions of `mod2` |
| `geneConc` | `dict` gene → concentration (mmol/gDW) |
| `x` | values of all columns |
| `geneCol` | DataFrame `gene`, `Col`: column of each gene (R: `geneCol`) |
| `prob` | the optlang problem |
| `mat` | the constraint matrix, in the format of [`getccFBA_mat`](getccFBA_mat.md#value) (`LHS`, `rlb`, `rub`, `rtype`, `rnames`, `clb`, `cub`, `cnames`, `obj_cf`, `obj_dir`, `geneCol`, `model`, `cc_ind`, `mod_objc_ind`) and `rxnkcatRow` (constrained reactions: `vind`, `rxn_id`, `dirxn`, `kval` (1/h), `kcat_src`, `rowind`, `nterms`). Change the bounds of `mat['model']` and call [`solve_ccFBA_mat(mat)`](getccFBA_mat.md#solve_ccfba_mat) to solve again. |
| `sol` | result of `solve_ccFBA_mat` |

R returns an unnamed list: the solution, `geneCol` and the problem object.

## Differences from R

- The problem is built as a matrix and solved with `solve_ccFBA_mat`, so the result can be solved
  again with other bounds.
- The kcat direction comes from the irreversible model (`is_bwd`), not from the `_f`/`_b` suffix of the id.
  Kcat ids are matched with or without the `R_` prefix.
- `s0001` is removed from all OR rules, not only from the rule strings that R knows.
- The rules are parsed with cobra (disjunctive normal form), so nested rules like `A and (B or C)` are
  expanded. R splits the rule text on `or`/`and`, which assumes rules that are already in this form.
- When neither `RHS` nor `objVal` is given, there is no budget constraint (R stops with an error). The
  default `selected_rxns` also works when `mod2` is not given (R uses `mod2` before making it).
- The LP file of `verboseMode > 2` is not written; use `res['prob'].to_lp()`.

## Check against the MATLAB model

iAF1260 from the R data set, data sets `kcat` and `mw`, glucose, `RHS = 0.27`, `medval = 3600*22.6`:
growth 0.508335. The MATLAB MOMENT model solved with [`readmodel`](readmodel.md) gives 0.508309.
[`cfba_moment`](cfba_moment.md) gives 0.524 for the same inputs, because it keeps `s0001` as a free
isoenzyme (see [conventions](../conventions.md#spontaneous-reactions-s0001)).

## References

Adadi R, et al. (2012), see [`cfba_moment`](cfba_moment.md#references).

## See also

[`cfba_moment`](cfba_moment.md), [`readmodel`](readmodel.md), [`getccFBA_mat`](getccFBA_mat.md).

## Examples

Python version of the example of the R help page (which calls `cfba_moment` instead of `cfba_moment_pw`):

```python
from ccfbapy import datasets
from ccfbapy import mod2irrev
from ccfbapy import cfba_moment_pw

model = datasets.load_model('iAF1260')       # R: data(iAF1260)
mw = datasets.load_dataset('mw')             # R: data(mw)
kcat = datasets.load_dataset('kcat')         # R: data(kcat)

model.reactions.get_by_id('R_EX_glc_e_').lower_bound = -1000   # uppbnd(mod2)[R_EX_glc_e__b] = 1000
model.reactions.get_by_id('R_EX_glyc_e_').lower_bound = 0
model.reactions.get_by_id('R_EX_ac_e_').lower_bound = 0
model.reactions.get_by_id('R_EX_o2_e_').lower_bound = -1000
model.reactions.get_by_id('R_ATPM').lower_bound = 0
mod2 = mod2irrev(model)

sol = cfba_moment_pw(model, mod2, kcat, MW=mw, verboseMode=1, RHS=0.27, medval=3600*22.6)
print(sol['status'], sol['obj'])
rk = sol['mat']['rxnkcatRow']
print(rk[rk['nterms'] > 2].head())                   # rules written as a pairwise chain
```

The chain gives the same solution as one row per rule, and the budget used is C:

```python
sol1 = cfba_moment_pw(model, mod2, kcat, MW=mw, verboseMode=0, RHS=0.27, medval=3600*22.6, pairwise=False)
print('pairwise:', sol['obj'], 'one row:', sol1['obj'])

# minimum budget for this growth rate (objVal)
sol2 = cfba_moment_pw(model, mod2, kcat, MW=mw, verboseMode=0, objVal=sol['obj'], medval=3600*22.6)
print('minimum budget:', sol2['obj'])               # about 0.27
```
