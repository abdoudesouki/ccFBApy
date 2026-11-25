# getGpr1iso: choose the cheapest isoenzyme

*R help page: `man/getGpr1iso.Rd`. Python: `ccfbapy/gpr1iso.py`, function `getGpr1iso()`.*

## Description

Preprocessing of the model for [`getccFBA_mat`](getccFBA_mat.md). The function computes the
molecular weight of every enzyme (gene or complex) in each GPR rule. For each reaction it keeps
only the isoenzyme with the smallest molecular weight, so the ORs are removed from the GPR
rules.

## Usage

```python
from ccfbapy import getGpr1iso

getGpr1iso(model, MW, cpx_stoich=None)
```

R (note the different order of the arguments):

```r
getGpr1iso(model, cpx_stoich, MW)
```

## Arguments

| Argument | Description |
|---|---|
| `model` | `cobra.Model` |
| `MW` | molecular weights in g/mmol: DataFrame, first column gene id, second column MW |
| `cpx_stoich` | optional stoichiometry of the complexes: DataFrame with the columns `genes` (sorted, comma separated gene ids) and `stoich` (comma separated number of subunits, same order). See [cpx_stoich](../conventions.md#protein-complex-stoichiometry-cpx_stoich). Default `None`: one subunit per gene. |

## Details

Each GPR rule is written in disjunctive normal form (OR of AND terms). The rule gets one of four
types:

| `gprtype` | Rule | Example |
|---|---|---|
| `r1g` | one gene | `b0001` |
| `r1cpx` | one complex | `b0001 and b0002` |
| `rmg` | isoenzymes, each one gene | `b0001 or b0002` |
| `rmcpx` | isoenzymes, at least one complex | `(b0001 and b0002) or b0003` |

The molecular weight of a term is `sum(MW_gene * stoich_gene)`. It is `sum(MW_gene)` when the
complex is not in `cpx_stoich`. For each reaction with more than one term, the term with the
smallest MW is kept. If two terms have the same MW, the first one is kept. As in R, a term with
an unknown MW (a gene that is not in `MW`, e.g. `s0001`) is preferred, because R's `min()` of
a group with `NA` is `NA`. The function warns about genes without a molecular weight.

A single-gene term gets the stoichiometry `1` unless `cpx_stoich` gives another value, e.g. `2` for a
homodimer.

## Value

A `dict`:

| Key | Content |
|---|---|
| `gpr1iso` | DataFrame, one row per reaction with a GPR: the chosen term |
| `mod_cpx_mw` | DataFrame, all terms of all GPR rules with their molecular weight (same columns) |
| `rgst` | DataFrame reactions × genes (index `model.reactions` ids, columns `model.genes` ids): stoichiometry of the genes of the chosen term, 0 otherwise |

Columns of `gpr1iso` and `mod_cpx_mw`:

| Column | Content |
|---|---|
| `rsn` | index of the reaction in `model.reactions` (0-based) |
| `rxn` | reaction id |
| `bool_expr` | the GPR rule (`rmcpx`: the term) |
| `genes` | sorted, comma separated genes of the term |
| `termNo` | number of the term in the rule (1-based) |
| `nOR`, `nAND` | number of ORs and ANDs in the rule |
| `gprtype` | `r1g`, `r1cpx`, `rmg` or `rmcpx` |
| `stoich` | comma separated subunit numbers of the genes |
| `termNG` | number of genes in the term |
| `mod_cpx_MW1s` | MW of the term with one subunit per gene |
| `mod_st_MW` | MW of the term with the stoichiometry of `cpx_stoich` (NaN when not available) |
| `mod_cpx_MWtxt` | MW of each gene, comma separated |
| `MW` | MW used: `mod_st_MW`, or `mod_cpx_MW1s` when it is NaN |

## See also

[`getccFBA_mat`](getccFBA_mat.md). [`ccFBA_iso`](ccFBA_iso.md) solves the same problem without choosing
one isoenzyme beforehand: the LP chooses, and knocked out or limited genes are handled.

## Examples

```python
import cobra
import pandas as pd
from ccfbapy import getGpr1iso

model = cobra.io.read_sbml_model('models/iAF1260.xml')
mw = pd.read_csv('data/Ecoli_MW.csv')
cpx_stoich = pd.read_csv('R/sybilccFBA/inst/extdata/cpx_stoich_me.csv')

r1i = getGpr1iso(model, MW=mw, cpx_stoich=cpx_stoich)
gpr1iso, mod_cpx_mw, rgst = r1i['gpr1iso'], r1i['mod_cpx_mw'], r1i['rgst']
print(gpr1iso['gprtype'].value_counts())

# aconitase: two isoenzymes, the lighter one (b0118) is chosen
print(mod_cpx_mw.loc[mod_cpx_mw['rxn'] == 'ACONTa', ['genes', 'stoich', 'MW']])
print(gpr1iso.loc[gpr1iso['rxn'] == 'ACONTa', ['genes', 'stoich', 'MW']])
print(rgst.loc['ACONTa', rgst.loc['ACONTa'] != 0])
```
