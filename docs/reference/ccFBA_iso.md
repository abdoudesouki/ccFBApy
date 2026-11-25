# ccFBA_iso: isoenzymes without getGpr1iso (gene allocation)

*No R help page: new in the Python version. Python: `ccfbapy/ccFBA_iso.py`, functions `ccFBA_iso()` and
`getccFBA_iso_mat()`.*

## Description

Cost constrained FBA in which the concentration of each gene is **allocated** to the enzymes that use
it. This solves the under estimation of the cost of rules like `(A and B) or (A and C)` without
[`getGpr1iso`](getGpr1iso.md). The isoenzymes are kept, and the LP chooses between them (or uses a
mixture). The function also includes multi-functional enzymes and the stoichiometry of protein
complexes, like [`getccFBA_mat`](getccFBA_mat.md).

## The problem

In MOMENT ([`cfba_moment`](cfba_moment.md)) each term of a rule gets its own limit. For
`(A and B) or (A and C)`:

```
v <= kcat * (cpx1 + cpx2)
cpx1 <= g_A,  cpx1 <= g_B          (complex A-B)
cpx2 <= g_A,  cpx2 <= g_C          (complex A-C)
```

Both complexes can use the whole concentration of A, so A is counted once but used twice. The enzyme
is cheaper than it really is, and the flux is over estimated. [`ccFBA_mr`](cfba_moment_mr.md) shares
a gene between *reactions*, but not between the complexes of the same reaction, so it has the same
problem.

[`getccFBA_mat`](getccFBA_mat.md) avoids it with [`getGpr1iso`](getGpr1iso.md): before the
problem is built, each reaction keeps only its cheapest isoenzyme (smallest molecular weight). This is a
preprocessing heuristic:

- the other isoenzymes are removed from the model. If a gene of the chosen isoenzyme is knocked out or
  limited (`geneConcUB`), the reaction is blocked or limited, although the other isoenzymes could
  catalyze it;
- the choice is made from the molecular weights only, before the problem is solved. It cannot
  depend on the condition, and mixtures of isoenzymes are not possible.

## The formulation

Every term `t` of the rule of reaction `r` (a gene or a complex) is an enzyme with its own
concentration `e_rt`. The concentration of each gene is the sum of its uses:

```
v_r  <= kcat_r * sum_t e_rt                                  for each direction of r (own kcat)
sum_{(r,t) that contain g} st_gt * e_rt  =  g_g              for each gene g
sum_g MW_g * g_g  <=  C
```

`st_gt` is the number of subunits of gene `g` in enzyme `t` (from `cpx_stoich`; 1 by default).

- A subunit shared by several complexes, of the same or of different reactions, is counted once for
  each complex. The cost is never under estimated.
- An enzyme that catalyzes several reactions gets one `e_rt` per reaction, and they share the
  gene: this is the multi-functional enzyme constraint.
- The problem is still an LP. The isoenzymes are chosen by the optimization, for each condition.

Every solution of `getccFBA_mat` is also a solution of this problem (set the other isoenzymes to 0).
Every solution of this problem is a solution of `ccFBA_mr` and MOMENT (set each complex variable to its
enzyme). So, with the same inputs, the growth rates are ordered:

```
getccFBA_mat  <=  ccFBA_iso  <=  ccFBA_mr  <=  cfba_moment
```

## Usage

```python
from ccfbapy import ccFBA_iso, getccFBA_iso_mat

getccFBA_iso_mat(model, mod2=None, Kcat=None, MW=None, selected_rxns=None, verboseMode=1, objVal=None,
                 RHS=None, medval=None, geneConcLB=None, geneConcUB=None, cpx_stoich=None, C_mu_coef=0,
                 gene_alloc='E', rm_genes=None)

ccFBA_iso(model, mod2=None, Kcat=None, MW=None, selected_rxns=None, verboseMode=1, objVal=None, RHS=None,
          medval=None, geneConcLB=None, geneConcUB=None, cpx_stoich=None, C_mu_coef=0, gene_alloc='E',
          rm_genes=None, mat=None, interface=None)
```

## Arguments

| Argument | Description |
|---|---|
| `model` | `cobra.Model` |
| `mod2` | irreversible model, `mod2irrev(model)` (default). The flux bounds are taken from it. See [mod2irrev](mod2irrev.md). |
| `Kcat` | **required**. kcat values in 1/s: DataFrame `rxn_id`, `val`, `dirxn`, optionally `src`; ids with or without `R_` (as in `getccFBA_mat`) |
| `MW` | molecular weights in g/mmol: first column gene, second column MW. Genes without MW have no cost. |
| `selected_rxns` | reactions of `model` to constrain (they must have a GPR). Default: all reactions with a GPR. |
| `verboseMode` | `0` silent, `1` status |
| `objVal`, `RHS` | `RHS` is the budget C (g/gDW). With `RHS=None`: minimize the budget needed for the objective value `objVal` (default: FBA optimum). |
| `medval` | kcat (1/h) for reactions without kcat. Default: median of `Kcat['val']`. |
| `geneConcLB`, `geneConcUB` | gene concentration bounds (mmol/gDW): `dict` gene → value, or a list in the order of `model.genes`. `geneConcUB={g: 0}` knocks out gene `g`. |
| `cpx_stoich` | stoichiometry of the complexes: DataFrame `genes`, `stoich` (see [conventions](../conventions.md#protein-complex-stoichiometry-cpx_stoich)). It is used for all the terms, also for single genes (e.g. a homodimer: `2`). |
| `C_mu_coef` | budget as a function of growth: `C = RHS + C_mu_coef * mu` |
| `gene_alloc` | `'E'` (default): each gene concentration equals the sum of its uses. `'U'`: the sum of its uses is at most the gene concentration, so an over expressed protein can stay unused. |
| `rm_genes` | genes to remove from the OR rules, e.g. `['s0001']` as in R's `cfba_moment`. A term with one of these genes is dropped when the rule has other terms. Default `None`: keep them; `s0001` has no MW and so no cost. |
| `mat` | (`ccFBA_iso`) a matrix from `getccFBA_iso_mat` to solve again. The flux bounds come from `mod2` (default `mat['model']`), and `geneConcLB`/`geneConcUB` are applied to it. |
| `interface` | (`ccFBA_iso`) optlang interface. Default: cobra's default solver. |

## Value

`getccFBA_iso_mat` returns a matrix in the format of [`getccFBA_mat`](getccFBA_mat.md#value), so it can be
solved with [`solve_ccFBA_mat`](getccFBA_mat.md#solve_ccfba_mat) (e.g. in a loop over conditions):

| Key | Content |
|---|---|
| `LHS`, `rlb`, `rub`, `rtype`, `rnames`, `clb`, `cub`, `cnames`, `obj_cf`, `obj_dir` | the LP. Columns: `R_<rxn>` (fluxes of `mod2`), `g_<gene>`, `enz_<rxn>_<t>` (term `t` of the rule, 1-based). Rows: `S_<met>`, `rkcat_<rxn of mod2>`, `galloc_<gene>`, `cc` or `ObjC_`. |
| `model` | the irreversible model `mod2` |
| `geneCol` | `gene`, `Col`, `MW`, `nenz` (number of enzymes that use the gene) |
| `enzCol` | one row per enzyme: `rsn`, `rxn`, `term`, `genes`, `stoich`, `nterms` (number of isoenzymes of the reaction), `MW` (of the enzyme, with the stoichiometry), `Col` |
| `rxnkcatRow` | kcat constraints: `vind`, `rxn_id`, `dirxn`, `kval` (1/h), `kcat_src`, `rowind` |
| `allocRow` | gene allocation rows: `gene`, `rowind`, `nuses` |
| `kcat_ir`, `kcat_src_ir`, `geneMW`, `fluxind`, `cc_ind`, `mod_objc_ind` | as in `getccFBA_mat` |

`ccFBA_iso` returns a `dict`:

| Key | Content |
|---|---|
| `obj`, `status` | objective value and solver status |
| `rxnFlx` | `dict` reaction id of `model` → net flux |
| `flx_ir` | fluxes of `mod2` |
| `geneConc` | `dict` gene → concentration (mmol/gDW) |
| `enzConc` | `enzCol` with `conc` (enzyme concentration), `MC` (= `MW * conc`, g/gDW) and `flx` (net flux of the reaction) |
| `rxnMC` | Series: protein cost of each reaction (sum of `MC` of its enzymes), largest first. With `gene_alloc='E'` the costs add up to the budget used. |
| `x`, `prob`, `mat` | values of all columns, the optlang problem, the matrix |

## Results on iAF1260

`models/iAF1260.xml`, `data/Ecoli_kcat.csv`, `data/Ecoli_MW.csv`, oxygen uptake open, ATPM = 0,
glucose, `RHS = 0.27`, `medval = 3600*22.6`:

| Method | Growth rate (1/h) |
|---|---|
| `cfba_moment` | 0.50581 |
| `ccFBA_mr` | 0.27393 |
| `ccFBA_iso` | 0.26254 |
| `getccFBA_mat` | 0.26254 |
| `ccFBA_iso`, complex stoichiometry (`cpx_stoich_me.csv`) | 0.12160 |
| `getccFBA_mat`, complex stoichiometry | 0.12160 |

Without gene limits, the cheapest isoenzyme is the best choice, so `ccFBA_iso` gives the same growth
as `getccFBA_mat`. This confirms that `getGpr1iso` is optimal in this case. The results differ
when genes are knocked out or limited.

**Knock-out scan.** For each of the 68 genes that `getGpr1iso` chose for an active reaction with
isoenzymes, the gene was knocked out (`geneConcUB = 0`):

| | `getccFBA_mat` | `ccFBA_iso` | FBA (cobra, GPR rules) |
|---|---|---|---|
| knock-outs without growth | 42 | 5 | 5 |
| without growth, but FBA grows | 37 | 0 | |
| growth higher than `getccFBA_mat` | | 61 | |

`getccFBA_mat` blocks a reaction when a gene of its chosen isoenzyme is knocked out, so it predicts 37
knock-outs as lethal although other isoenzymes exist. `ccFBA_iso` predicts the same lethal genes as FBA
(b2323, b2600, b2762, b3040, b3770). For half of the genes the growth stays within 0.01% of the wild
type, because the other isoenzyme costs about the same. For acnB (b0118) the growth is 0.26252 instead
of 0, because acnA (b1276) takes over.

## See also

[`getGpr1iso`](getGpr1iso.md), [`getccFBA_mat`](getccFBA_mat.md), [`cfba_moment_mr`](cfba_moment_mr.md),
[`cfba_moment`](cfba_moment.md), [`simGeneOverExpr`](simGeneOverExpr.md).

## Examples

A toy model with one reaction and the rule `(A and B) or (A and C)`: kcat 1/s, all MW 1, budget 0.1.
The enzyme costs 2 per unit (A + B), so the right flux is `3600 * 0.1 / 2 = 180`:

```python
import cobra
import pandas as pd
from ccfbapy import cfba_moment
from ccfbapy import ccFBA_mr
from ccfbapy import getccFBA_mat, solve_ccFBA_mat
from ccfbapy import ccFBA_iso

toy = cobra.Model('toy')
s, p = cobra.Metabolite('s', compartment='c'), cobra.Metabolite('p', compartment='c')
ex_s = cobra.Reaction('EX_s', lower_bound=-1000, upper_bound=1000)
ex_s.add_metabolites({s: -1})
r1 = cobra.Reaction('R1', lower_bound=0, upper_bound=1000)
r1.add_metabolites({s: -1, p: 1})
r1.gene_reaction_rule = '(A and B) or (A and C)'
ex_p = cobra.Reaction('EX_p', lower_bound=0, upper_bound=1000)
ex_p.add_metabolites({p: -1})
toy.add_reactions([ex_s, r1, ex_p])
toy.objective = 'EX_p'

kcat = pd.DataFrame({'rxn_id': ['R_R1'], 'val': [1.0], 'dirxn': [1]})
mw = pd.DataFrame({'gene': ['A', 'B', 'C'], 'mw': [1.0, 1.0, 1.0]})
kw = dict(Kcat=kcat, MW=mw, RHS=0.1, verboseMode=0)
print('cfba_moment ', cfba_moment(toy, **kw)['obj'])      # 240: A is used by both complexes
print('ccFBA_mr    ', ccFBA_mr(toy, **kw)['obj'])         # 240
print('getccFBA_mat', solve_ccFBA_mat(getccFBA_mat(toy, **kw), verboseMode=0)['obj'])   # 180
res = ccFBA_iso(toy, **kw)
print('ccFBA_iso   ', res['obj'])                          # 180
print(res['enzConc'][['rxn', 'genes', 'MW', 'conc', 'MC']])
```

Knock out B. `getccFBA_mat` only has the complex A-B (the first of two equal isoenzymes), while
`ccFBA_iso` switches to A-C:

```python
a = solve_ccFBA_mat(getccFBA_mat(toy, geneConcUB={'B': 0}, **kw), verboseMode=0)['obj']
res = ccFBA_iso(toy, geneConcUB={'B': 0}, **kw)
print('B knocked out: getccFBA_mat %.1f, ccFBA_iso %.1f' % (a, res['obj']))    # 0 and 180
print(res['enzConc'][['genes', 'conc']].to_string(index=False))
```

Make C heavier (MW 1.5) and limit B to 0.02. `ccFBA_iso` uses A-B up to the limit, then A-C for the rest
of the budget: `3600 * (0.02 + 0.06 / 2.5) = 158.4`. `getccFBA_mat` is limited by B: `3600 * 0.02 = 72`.

```python
kw2 = dict(kw, MW=pd.DataFrame({'gene': ['A', 'B', 'C'], 'mw': [1.0, 1.0, 1.5]}))
a = solve_ccFBA_mat(getccFBA_mat(toy, geneConcUB={'B': 0.02}, **kw2), verboseMode=0)['obj']
res = ccFBA_iso(toy, geneConcUB={'B': 0.02}, **kw2)
print('B <= 0.02: getccFBA_mat %.1f, ccFBA_iso %.1f' % (a, res['obj']))        # 72 and 158.4
print(res['enzConc'][['genes', 'MW', 'conc']].to_string(index=False))
```

iAF1260: the wild type, and a knock out of the isoenzyme chosen by `getGpr1iso` for aconitase
(`b0118`, acnB). The matrix is built once and solved with other gene bounds (`mat=`):

```python
from ccfbapy import mod2irrev
from ccfbapy import getccFBA_iso_mat

model = cobra.io.read_sbml_model('models/iAF1260.xml')
mw = pd.read_csv('data/Ecoli_MW.csv')
kcat = pd.read_csv('data/Ecoli_kcat.csv')
model.reactions.get_by_id('EX_o2_e').lower_bound = -1000
model.reactions.get_by_id('ATPM').bounds = (0, 0)
mod2 = mod2irrev(model)
kw = dict(Kcat=kcat, MW=mw, RHS=0.27, medval=3600*22.6, verboseMode=0)

mat = getccFBA_mat(model, mod2=mod2, **kw)
imat = getccFBA_iso_mat(model, mod2=mod2, **kw)
wt = ccFBA_iso(model, mat=imat, verboseMode=0)
print('wild type: getccFBA_mat %.5f, ccFBA_iso %.5f' % (solve_ccFBA_mat(mat, verboseMode=0)['obj'], wt['obj']))
print(wt['rxnMC'].head(), 'budget used:', wt['rxnMC'].sum())

ko = {'b0118': 0}                                   # acnB
mat_ko = dict(mat, cub=mat['cub'].copy())
mat_ko['cub'][int(mat['geneind']['b0118'])] = 0
res_ko = ccFBA_iso(model, mat=imat, geneConcUB=ko, verboseMode=0)
print('acnB KO: getccFBA_mat %.5f, ccFBA_iso %.5f' % (solve_ccFBA_mat(mat_ko, verboseMode=0)['obj'], res_ko['obj']))
e = res_ko['enzConc']
print(e.loc[e['rxn'].isin(['ACONTa', 'ACONTb']), ['rxn', 'genes', 'MW', 'conc', 'flx']])
```
