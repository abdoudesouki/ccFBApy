# Data sets

The R package has nine data sets in `R/sybilccFBA/data/*.RData`, each with a help page in `man/`. ccFBApy
ships them as **SBML** (models) and **CSV** (tables) files, read with
[`ccfbapy.datasets`](#data-sets-in-the-package). The package does not read `.RData` files and does not need
R or the package `rdata`: the files were exported once from the `.RData` files (see
[How the data sets were made](#how-the-data-sets-were-made)).

```python
from ccfbapy import datasets
iAF1260 = datasets.load_model('iAF1260')                              # R: data(iAF1260)
print(iAF1260.id, len(iAF1260.reactions), 'reactions')
```

| R data set | Content | Python | File in the package (`ccfbapy/data/`) |
|---|---|---|---|
| [`iAF1260`](#iaf1260) | *E. coli* model | `datasets.load_model('iAF1260')` | `sybilccFBA/iAF1260.xml.gz` |
| [`iJO1366`](#ijo1366) | *E. coli* model | `datasets.load_model('iJO1366')` | `sybilccFBA/iJO1366.xml.gz` |
| [`iML1515`](#iml1515) | *E. coli* model | `datasets.load_model('iML1515')` | `sybilccFBA/iML1515.xml.gz` |
| [`iMM904`](#imm904) | *S. cerevisiae* model | `datasets.load_model('iMM904')` | `sybilccFBA/iMM904.xml.gz` |
| [`kcat`](#kcat) | *E. coli* kcat values | `datasets.load_dataset('kcat')` | `sybilccFBA/kcat.csv` |
| [`mw`](#mw) | *E. coli* molecular weights | `datasets.load_dataset('mw')` | `sybilccFBA/mw.csv` |
| [`mw_iML1515`](#mw_iml1515) | molecular weights of the iML1515 genes | `datasets.load_dataset('mw_iML1515')` | `sybilccFBA/mw_iML1515.csv` |
| [`yst_kcat`](#yst_kcat) | *S. cerevisiae* kcat values | `datasets.load_dataset('yst_kcat')` | `sybilccFBA/yst_kcat.csv` |
| [`yst_mw`](#yst_mw) | *S. cerevisiae* molecular weights | `datasets.load_dataset('yst_mw')` | `sybilccFBA/yst_mw.csv` |

## Two versions of the models

There are two SBML versions of each model:

- `datasets.load_model(name)` (default, `source='sybilccFBA'`): **the data set of the R package**, with the ids of
  sybil, its bounds and GPR rules. The results are the same as with the R data set (e.g. the tests of the
  package: `cfba_moment_pw` on iAF1260 gives 0.508335 with both).
- `datasets.load_model(name, source='bigg')`, or `cobra.io.read_sbml_model('models/<name>.xml')` in the
  repository: the SBML models of the folder `models/`, with BiGG ids. They are another version of the models:
  with the default bounds the FBA optimum is the same, but some ids and bounds differ (e.g. upper bounds 999999
  instead of 1000, ATPM). The examples of this manual that use `models/` are written for this version.

| Model | R data set (`source='sybilccFBA'`): reactions / metabolites / genes | example ids | BiGG (`source='bigg'`, `models/*.xml`) | example ids |
|---|---|---|---|---|
| iAF1260 | 2382 / 1668 / 1261 | `R_EX_glc_e_`, `R_ATPM` | 2382 / 1668 / 1261 | `EX_glc__D_e`, `ATPM` |
| iJO1366 | 2583 / 1805 / 1367 | `EX_glc(e)`, `ATPM` | 2583 / 1805 / 1367 | `EX_glc__D_e` |
| iML1515 | 2712 / 1877 / 1516 | `EX_glc__D(e)`, `ATPM` | 2712 / 1877 / 1516 | `EX_glc__D_e` |
| iMM904 | 1577 / 1228 / 904 | `R_EX_glc_e_`, `R_ATPM` | 1577 / 1226 / 905 | `EX_glc__D_e` |

The gene counts of the *E. coli* models include the spontaneous pseudo gene `s0001`. The kcat data set `kcat`
(ids of the R iAF1260) matches 508 reactions of the R iAF1260 and 501 of the BiGG version. With glucose and
oxygen open, `cfba_moment_pw` gives 0.508335 with the R version (as the MATLAB MOMENT model, 0.508309) and
0.511998 with the BiGG version.

## iAF1260

*Escherichia coli* metabolic model iAF1260 (`man/iAF1260.Rd`).

A genome scale metabolic network of *E. coli*: 2077 internal reactions, 304 exchange reactions and
a biomass objective function.

- **Usage** R: `data(iAF1260)`. Python: `datasets.load_model('iAF1260')`. BiGG version:
  `datasets.load_model('iAF1260', source='bigg')` or `cobra.io.read_sbml_model('models/iAF1260.xml')`.
- **Format** R: an object of class `modelorg`. Python: `cobra.Model`.
- **Reference** Feist AM, Henry CS, Reed JL, Krummenacker M, Joyce AR, Karp PD, Broadbelt LJ,
  Hatzimanikatis V, Palsson BØ (2007). A genome-scale metabolic reconstruction for *Escherichia
  coli* K-12 MG1655 that accounts for 1260 ORFs and thermodynamic information. *Mol Syst Biol* 3: 121.

## iJO1366

*Escherichia coli* metabolic model iJO1366 (`man/iJO1366.Rd`).

A genome scale metabolic network of *E. coli*: 2253 internal reactions, 330 exchange reactions and
a biomass objective function.

- **Usage** R: `data(iJO1366)`. Python: `datasets.load_model('iJO1366')`. BiGG version:
  `datasets.load_model('iJO1366', source='bigg')` or `cobra.io.read_sbml_model('models/iJO1366.xml')`.
- **Format** R: `modelorg`. Python: `cobra.Model`.
- **Reference** Orth JD, Conrad TM, Na J, Lerman JA, Nam H, Feist AM, Palsson BØ (2011). A
  comprehensive genome-scale reconstruction of *Escherichia coli* metabolism—2011. *Mol Syst Biol* 7(1).

## iML1515

*Escherichia coli* metabolic model iML1515 (`man/iML1515.Rd`).

A genome scale metabolic network of *E. coli*: 2381 internal reactions, 331 exchange reactions and
a biomass objective function. It has 1515 genes. Used in the example of
[`simulate_EColi`](reference/simulate_EColi.md).

- **Usage** R: `data(iML1515)`. Python: `datasets.load_model('iML1515')`. BiGG version:
  `datasets.load_model('iML1515', source='bigg')` or `cobra.io.read_sbml_model('models/iML1515.xml')`.
- **Format** R: `modelorg`. Python: `cobra.Model`.
- **Reference** Monk JM, et al. (2017). iML1515, a knowledgebase that computes *Escherichia coli*
  traits. *Nature Biotechnology* 35(10): 904.

## iMM904

*Saccharomyces cerevisiae* metabolic model (`man/iMM904.Rd`).

A genome scale metabolic network of *S. cerevisiae*: 1412 internal reactions, 164 exchange
reactions and a biomass objective function.

- **Usage** R: `data(iMM904)`. Python: `datasets.load_model('iMM904')`. BiGG version:
  `datasets.load_model('iMM904', source='bigg')` or `cobra.io.read_sbml_model('models/iMM904.xml')`.
- **Format** R: `modelorg`. Python: `cobra.Model`.
- **Reference** Mo ML, Palsson BØ, Herrgård MJ (2009). Connecting extracellular metabolomic
  measurements to intracellular flux states in yeast. *BMC Syst Biol* 3: 37.

## kcat

*Escherichia coli* kcat values used in the MOMENT method (`man/kcat.Rd`).

A list of kcat values used in MOMENT, in 1/s.

- **Usage** R: `data(kcat)`. Python: `datasets.load_dataset('kcat')`. The same values are in
  `data/Ecoli_kcat.csv` (`datasets.load_kcat()`).
- **Format** a data frame with three columns and 508 rows: `rxn_id` (iAF1260 reaction id with the
  `R_` prefix, e.g. `R_ASPTA`), `val` (kcat, 1/s), `dirxn` (`1`: 452 rows, `-1`: 56 rows). It can be
  given as `Kcat` to all the functions (see [Kcat](conventions.md#kinetic-data-kcat)).
- **Reference** Adadi R, Volkmer B, Milo R, Heinemann M, Shlomi T (2012). Prediction of microbial
  growth rate versus biomass yield by a metabolic network with kinetic parameters. *PLoS Comput
  Biol* 8(7). doi:10.1371/journal.pcbi.1002575

## mw

*Escherichia coli* molecular weights used in the MOMENT method (`man/mw.Rd`).

- **Usage** R: `data(mw)`. Python: `datasets.load_dataset('mw')`.
- **Format** a data frame with two columns and 1260 rows: `Synonym` (b-number) and `mw`. The help
  page gives the unit as g/mol, but the values are in g/mmol (kDa), e.g. b0002: 89.12 (see
  [units](conventions.md#units)). `s0001` is not included.
- `data/Ecoli_MW.csv` (`datasets.load_mw()`, columns `Synonym`, `mw`) has 1548 genes: the 1260 genes of `mw`, and 288
  more. For 7 genes the values differ by up to 2%.
- **Reference** Adadi et al. (2012), see [kcat](#kcat).

## mw_iML1515

Molecular weights of the *E. coli* metabolic genes of iML1515, used in MOMENT/ccFBA (`man/mw_iML1515.Rd`).

Calculated with [`calc_MW`](reference/calc_MW.md) from the *E. coli* genome sequence available at NCBI.

- **Usage** R: `data(mw_iML1515)`. Python: `datasets.load_dataset('mw_iML1515')`.
- **Format** a data frame with two columns and 1515 rows: `Synonym` and `mw` (g/mmol). The
  example of [`simulate_EColi`](reference/simulate_EColi.md) adds `s0001` with MW 0.001 and
  renames the first column to `gene`.
- **Reference** Adadi et al. (2012), see [kcat](#kcat).

## yst_kcat

*Saccharomyces cerevisiae* kcat values used in the MOMENT method (`man/yst_kcat.Rd`).

A list of kcat values for yeast, in 1/s.

- **Usage** R: `data(yst_kcat)`. Python: `datasets.load_dataset('yst_kcat')`. The same values are in
  `data/yeast_imm_kcat.csv` (`datasets.load_yeast_kcat()`).

    ```python
    from ccfbapy import datasets
    yst_kcat = datasets.load_dataset('yst_kcat')
    print(yst_kcat.head())
    ```

- **Format** the help page says "a data frame with three columns". In R it is a numeric matrix with
  938 rows and the columns `rxn_id`, `val`, `dirxn` (a DataFrame in Python). `rxn_id` is a number (1–469, each with
  `dirxn` 1 and −1), not a reaction id. The help page does not say which reactions the numbers
  refer to, so check this before you use the values with `iMM904`.
- **References** Richter C (2011). *Kosten und Effizienz von Enzymen als zusätzliche Bedingung für
  die Flussverteilung in metabolischen Netzwerken.* Diploma thesis, Humboldt University Berlin.
  Adadi et al. (2012), see [kcat](#kcat).

## yst_mw

*Saccharomyces cerevisiae* molecular weights used in the MOMENT method (`man/yst_mw.Rd`).

Calculated with [`calc_MW`](reference/calc_MW.md) from the yeast genome sequence available at NCBI
(*S. cerevisiae* S288c).

- **Usage** R: `data(yst_mw)`. Python: `datasets.load_dataset('yst_mw')`. The same values are in
  `data/yeast_mw.csv` (`datasets.load_yeast_mw()`).
- **Format** a data frame with 5899 rows and the columns `gene` (systematic name, e.g. `YAL005C`),
  `mw` (g/mmol) and `LENGTH` (number of amino acids).
- **Reference** Adadi et al. (2012), see [kcat](#kcat).

## Data sets in the package

The module `ccfbapy.datasets` gives the data of the package, so they can be used without a copy of
the repository. Only SBML and CSV files are read (cobra and pandas).

| Function | Returns |
|---|---|
| `load_model(name, source='sybilccFBA')` | a model as `cobra.Model`: `iAF1260`, `iJO1366`, `iML1515`, `iMM904`. `source='sybilccFBA'`: the R data set (sybil ids); `source='bigg'`: the model of `models/` (BiGG ids), see [Two versions of the models](#two-versions-of-the-models) |
| `load_dataset(name)` | an R data set: the models (as `load_model(name)`) and the tables `kcat`, `mw`, `mw_iML1515`, `yst_kcat`, `yst_mw` (DataFrame) |
| `list_datasets()` | the names of the data sets |
| `model_path(name, source='sybilccFBA')` | the path of the SBML file (gzip) of a model |
| `load_kcat()`, `load_mw()` | *E. coli* kcat values (1/s) and molecular weights (g/mmol): `data/Ecoli_kcat.csv`, `data/Ecoli_MW.csv` |
| `load_yeast_kcat()`, `load_yeast_mw()` | the same for *S. cerevisiae*: `data/yeast_imm_kcat.csv`, `data/yeast_mw.csv` |
| `load_cpx_stoich()` | protein complex stoichiometry, `cpx_stoich_me.csv` (for `cpx_stoich`) |
| `load_kcat_list(model=None, column='kcat_max')` | `allKcats_upd34_dd_h.csv` as a Kcat table (`rxn_id`, `val`, `dirxn`, `src`), as in the example of `simulate_EColi`; only the reactions of `model` when it is given |
| `extdata_path(fname=None)` | the folder of the example files of sybilccFBA (`inst/extdata`, see below), or the path of one of them. [`calc_MW`](reference/calc_MW.md) and [`readmodel`](reference/readmodel.md) use these files by default. |
| `data_path(fname=None)` | the data folder of the package |

```python
from ccfbapy import datasets

model = datasets.load_model('iML1515')
kcat = datasets.load_kcat_list(model)            # kcat_max, reactions of iML1515
cpx_stoich = datasets.load_cpx_stoich()
print(model.id, len(model.reactions), kcat.shape, cpx_stoich.shape)
print(datasets.load_kcat().head(3))
bigg = datasets.load_model('iML1515', source='bigg')
print(bigg.reactions.get_by_id('EX_glc__D_e').name)
```

## How the data sets were made

The development tool `tools/make_package_data.py` of the repository prepares `ccfbapy/data`. It copies the CSV
tables of `data/` and the files of `R/sybilccFBA/inst/extdata`, compresses the SBML models of `models/`
(`bigg/*.xml.gz`), and exports the data sets of the R package to `sybilccFBA/`: the models with
`cobra.io.write_sbml_model` (gzip), the tables as CSV. Only this export reads `.RData` files: it needs
`pip install rdata` and `tools/sybil2cobra.py`, not the package.

Each export is read again and compared with the R data set:

- models: reaction, metabolite and gene ids and their order, bounds, stoichiometry, GPR rules, objective,
  names and formulas, and the FBA optimum;
- tables: all values, exactly (the package reads the CSV files with `float_precision='round_trip'`).

Three things SBML cannot store as in R:

- in the R data sets of iAF1260 and iMM904 the metabolites have no compartment (`met_comp` is NA, ids such as
  `M_glc_D_c[NA]`). SBML needs one: it is taken from the id (the letter before `[`, e.g. `c`); the ids are kept;
- SBML writes the GPR rules without redundant parentheses (`(A and B) and C` becomes `A and B and C`, the
  same rule);
- a missing metabolite charge (iAF1260) is read as 0.

None of them changes the results.

## Other data files

Files in `R/sybilccFBA/inst/extdata` (R: `system.file("extdata", ..., package="sybilccFBA")`). They are also in
the package: `datasets.extdata_path(fname)`.

| File | Content | Used by |
|---|---|---|
| `allKcats_upd34_dd_h.csv` | *E. coli* kcat values from several sources: `ijo_id` (reaction id without `R_`), `iaf_id`, `dirxn`, `src`, `kcat_max`, ... | examples of `simulate_EColi`, package help page |
| `cpx_stoich_me.csv` | protein complex stoichiometry: `genes`, `stoich`, `Enzyme.ID`, `Gene.composition`, ... | `cpx_stoich` argument |
| `Ec_core_genedef.csv`, `Ec_iAF1260_genedef.csv`, `iJO1366_genedef.csv` | gene definitions of *E. coli* models: gene, molecular weight (`mw`/`MW`, g/mmol), gene name, product, ... | examples of the R help pages (`Ec_core`), gene names in `test_simulate_EColi.ipynb` |
| `yst_genedef.csv` | yeast gene definitions (`mw` in g/mol) | |
| `aa.txt`, `test2.ptt`, `NC_000913.faa`, `ALL_PTT.TSV` | amino acid weights, protein tables (*E. coli* `test2.ptt`, yeast `ALL_PTT.TSV`) and *E. coli* protein sequences | [`calc_MW`](reference/calc_MW.md) |
| `mat.csv`, `mets.csv`, `rxns.csv`, `rbnds.csv`, `cbnds.csv` | MOMENT model exported from MATLAB (iAF1260) | [`readmodel`](reference/readmodel.md), read with `read_moment_files()` |
| `iJO1366_mod2_mppv.csv`, `iJO_M_prot_pool.csv` | protein pool values of the iJO1366 reactions (`MppVal`) | |

`R/OverExpr/kcat2` is the kcat file of `test_simGeneOverExpr.R` (columns `rxn_id`, `val`,
`dirxn`, tab separated, no header).
