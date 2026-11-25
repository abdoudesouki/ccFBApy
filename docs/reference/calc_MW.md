# calc_MW: molecular weights from the genome

*R help page: `man/calc_MW.Rd`. Python: `ccfbapy/calc_MW.py`, function `calc_MW()`.*

## Description

Calculates the molecular weights of the proteins of a genome from their sequences (`.faa` file).

## Usage

```python
from ccfbapy import calc_MW

calc_MW(aa_fname="aa.txt", ptt_fname="test2.ptt", faa_fname="NC_000913.faa", nchrm=None)
```

R:

```r
calc_MW(aa_fname = "aa.txt", ptt_fname = "test2.ptt", faa_fname = "NC_000913.faa", nchrm = 1)
```

## Arguments

| Argument | Description |
|---|---|
| `aa_fname` | file with the list of amino acids and their molecular weights: tab separated, columns `AA_1` (one letter code) and `MW` (g/mol) |
| `ptt_fname` | protein table (`.ptt`, tab separated) with the protein ids (`PID`) and the gene names (e.g. `Synonym`: b-numbers) |
| `faa_fname` | protein sequences in FASTA format; the protein id (`PID`) is the second field of the header line (fields separated by vertical bars). One file per chromosome: a file name or a list of file names. |
| `nchrm` | number of chromosomes (number of `.faa` files used). Default: all the files in `faa_fname`. In R the default is 1. |

If a file name does not exist, the function looks for the file in `R/sybilccFBA/inst/extdata`
(as `system.file("extdata", ...)` in R). The defaults are the *E. coli* K-12 MG1655 files there.

## Details

Molecular weight = sum of the weights of the residues − water (18.015 g/mol) × (number of residues − 1).
Only the 20 standard amino acids are counted (other letters, e.g. `X`, are ignored, as in R).

## Value

A DataFrame with one row per protein found in both the `.faa` file and the `.ptt` table (merged on
`PID`, sorted by `PID`):

| Column | Content |
|---|---|
| `PID` and the other columns of the `.ptt` file | e.g. `LOCATION`, `STRAND`, `LENGTH`, `GENE`, `Synonym`, `COG`, `PRODUCT` |
| `gsn` | serial number of the protein in its `.faa` file |
| `mw` | molecular weight in **g/mol** |
| `A`, `C`, `D`, ..., `Y` | number of each amino acid |

The ccFBA functions need the molecular weight in g/mmol (kDa), so divide `mw` by 1000. For
example, b0002 is 89119.9 g/mol, and 89.12 in the data sets.

## Differences from R

- R returns all columns as text (the counts are built with `cbind`); Python returns numbers.
- R does not record the last protein of a `.faa` file that has only one protein. Python does.
- `nchrm` defaults to the number of files given.
- Missing files are looked up in the `extdata` folder.

## Check against the data sets

For the 1514 genes of `mw_iML1515` that are in `test2.ptt`, `calc_MW()['mw'] / 1000` equals the data
set (the difference is below 1e-13). `b4104` is not in `test2.ptt`. The older data set `mw` (1260
genes, used by MOMENT) was made with slightly different values and differs by up to about 3.4 g/mol
(0.0034 kDa).

## See also

[Data sets](../data.md) `mw`, `mw_iML1515`, `yst_mw`; [units](../conventions.md#units).

## Examples

Python version of the example of the R help page:

```python
from ccfbapy import calc_MW

aa_fname = 'R/sybilccFBA/inst/extdata/aa.txt'
ptt_fname = 'R/sybilccFBA/inst/extdata/test2.ptt'
faa_fname = 'R/sybilccFBA/inst/extdata/NC_000913.faa'

geneCnt = calc_MW(aa_fname, ptt_fname, faa_fname)
print(geneCnt[['PID', 'Synonym', 'GENE', 'LENGTH', 'gsn', 'mw']].head())
# geneCnt.to_csv('geneCnt.csv')                  # R: write.csv(file="geneCnt.csv", geneCnt)
```

Molecular weights in the format of the `MW` argument of the ccFBA functions (g/mmol), and a check
against the data set `mw_iML1515`:

```python
import pandas as pd
from ccfbapy import datasets

mw = pd.DataFrame({'gene': geneCnt['Synonym'], 'mw': geneCnt['mw'] / 1000})
mw_iML1515 = datasets.load_dataset('mw_iML1515')
cmp = mw_iML1515.merge(mw, left_on='Synonym', right_on='gene', suffixes=('_data', '_calc'))
print(len(cmp), 'genes, max difference:', (cmp['mw_data'] - cmp['mw_calc']).abs().max())
```
