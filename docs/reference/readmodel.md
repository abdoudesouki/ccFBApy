# readmodel: read the MOMENT model of MATLAB

*R help page: `man/readmodel.Rd`. Python: `ccfbapy/readmodel.py`, functions `readmodel()` and
`read_moment_files()`.*

## Description

Creates a linear problem from the lists exported from the MATLAB implementation of MOMENT, solves it,
and returns the fluxes. It is used to compare the R and Python implementations with the original
MOMENT model.

## Usage

```python
from ccfbapy import readmodel, read_moment_files

files = read_moment_files(path='R/sybilccFBA/inst/extdata')
readmodel(mat, mets, rxns, rbnds, cbnds, solver="glpk", C=0.27, lp_fname=None, verboseMode=1)
```

R:

```r
readmodel(mat, mets, rxns, rbnds, cbnds, solver = "glpkAPI")
```

## Arguments

| Argument | Description |
|---|---|
| `mat` | the constraint matrix as triplets: row, column, value (1-based indices, as in MATLAB and R) |
| `mets` | metabolites (first column: ids). There is one equality row for each, and they are the first rows. |
| `rxns` | reactions (first column: ids) and their bounds; column `ocf` holds the objective coefficients. There is one flux column for each, and they are the first columns. |
| `rbnds` | bounds of the rows (`lb`, `ub`) |
| `cbnds` | bounds of the columns (`lb`, `ub`) |
| `solver` | optlang interface: `'glpk'` (default), `'cplex'`, `'gurobi'`, ..., or an interface module. The R names `"glpkAPI"` and `"cplexAPI"` are accepted. |
| `C` | budget: the upper bound of the last row (the enzyme capacity constraint). Default 0.27, as fixed in R. `None` keeps the value of `rbnds` (0.56 in the example files). |
| `lp_fname` | when given, the problem is written to this file in LP format. R always writes an LP file (`glpk_eFBA_<date>.lp` or `Cplex_moment_<date>.lp`) to the working folder. |
| `verboseMode` | `0` silent; `1` prints the size, status, objective and the fluxes of some exchange reactions (as R) |

`read_moment_files(path)` reads the five files of the MATLAB model from `path` (default
`R/sybilccFBA/inst/extdata`): `mat.csv`, `mets.csv`, `rxns.csv`, `rbnds.csv`, `cbnds.csv`
(separator `;`). It returns a `dict` with the keys `mat`, `mets`, `rxns`, `rbnds`, `cbnds`, so you can call
`readmodel(**files)`.

## Details

- The rows `1..len(mets)` are `S v = b`, fixed to `rbnds.lb`. The other rows are `<= rbnds.ub`,
  except the last row (the enzyme capacity), which is `<= C`.
- The columns have the bounds of `cbnds`. The first `len(rxns)` columns are the fluxes; the others are
  the enzyme and auxiliary variables of MOMENT.
- The objective `rxns$ocf` (the biomass reaction) is maximized.

The example files are the *E. coli* iAF1260 MOMENT model: 6705 rows (1674 metabolites) and 4991
columns (3234 fluxes), with glucose uptake (`REV_EX_glc(e)`) open. With `C = 0.27` the growth rate is
0.508 1/h and the glucose uptake 12.6 mmol/gDW/h.

## Value

A `dict`:

| Key | Content |
|---|---|
| `x` | values of all columns (numpy array). R returns this vector. |
| `obj`, `status` | objective value (growth rate) and solver status |
| `rxnFlx` | Series: values of the flux columns, with the reaction ids of `rxns` as index |
| `prob` | the optlang problem. Columns are named with the reaction ids and `x<j>` (other columns, `j` 1-based); rows with the metabolite ids, `r<i>`, and `MW` (last row). |

## Differences from R

- R fixes the size of the problem in the code (6705 rows, 1674 metabolites, 4991 columns, 3234 fluxes).
  Python takes it from the inputs, and these are the sizes of the example files.
- The file `mat.csv` has no header line. The comments of the R source read it with
  `read.table(..., header=T)`, which takes the first entry of the matrix as the header and so drops it.
  `read_moment_files` reads it as data.
- Six lines of `mets.csv` (alcohols added at the end, e.g. `isobutanol`) have one field more.
  `read_moment_files` drops the extra field.
- With `cplexAPI`, R also sets the upper bounds of five exchange columns (acetate, pyruvate, fructose,
  lactate: 0; glucose: 1000). These are the values already in `cbnds.csv`, so both R solvers, and
  Python, solve the same problem.
- The budget can be changed (`C`), and the LP file is only written when `lp_fname` is given.

## References

Adadi R, et al. (2012), see [`cfba_moment`](cfba_moment.md#references).

## See also

[`cfba_moment`](cfba_moment.md)

## Examples

```python
from ccfbapy import read_moment_files, readmodel

files = read_moment_files()                      # R/sybilccFBA/inst/extdata
print({k: v.shape for k, v in files.items()})
res = readmodel(**files)                         # budget C = 0.27, as in R
print(res['status'], res['obj'])
print(res['rxnFlx'][['REV_EX_glc(e)', 'EX_o2(e)', 'REV_EX_o2(e)']])

# the same model with a smaller budget
res2 = readmodel(**files, C=0.2, verboseMode=0)
print('C=0.2:', res2['obj'])
```
