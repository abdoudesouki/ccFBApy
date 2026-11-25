# mod2irrev: irreversible model

*No R help page: this is sybil's `mod2irrev()` / class `modelorg_irrev`. Python: `ccfbapy/getccFBA_mat.py`,
function `mod2irrev()` and class `IrrevModel`.*

## Description

Makes an irreversible version of a cobra model, as sybil's `mod2irrev`. A reaction that can carry
a negative flux (`lower_bound < 0`) is split into a forward reaction `<id>_f` and a backward
reaction `<id>_b`, placed next to each other. The other reactions keep their id.
[`getccFBA_mat`](getccFBA_mat.md), [`getccFBA_sol`](getccFBA_sol.md) and
[`simulate_EColi`](simulate_EColi.md) use it. It plays the role of the `mod2` argument of the R
functions.

## Usage

```python
from ccfbapy import mod2irrev, IrrevModel

mod2 = mod2irrev(model, rev=None)      # same as IrrevModel(model, rev=None)
```

| Argument | Description |
|---|---|
| `model` | `cobra.Model` |
| `rev` | optional list of reaction ids to split even though their lower bound is ≥ 0. In R you would set `react_rev(model)[...] = TRUE` before `mod2irrev`. Their backward reaction gets the bounds `(0, -lower_bound)`, e.g. `(0, 0)`, and can be opened later with `uppbnd`. |

## Value / attributes

`n` is the number of reactions of the irreversible model.

| Attribute | Content |
|---|---|
| `react_id`, `react_name` | ids and names (length `n`) |
| `lowbnd`, `uppbnd` | bounds (numpy arrays, length `n`). You can change them before solving, as `uppbnd(mod2)` in R. |
| `obj_coef` | objective coefficients (the backward reaction gets `-c`) |
| `obj_dir` | `'max'` or `'min'` |
| `S` | stoichiometric matrix (`scipy.sparse`, metabolites × `n`) |
| `met_id` | metabolite ids |
| `react_num`, `met_num` | number of reactions and metabolites |
| `irrev2rev` | index of the original reaction of each reaction |
| `rev2irrev` | array (reactions of `model` × 2): forward and backward index (both equal to the index of the reaction when it is not split) |
| `matchrev` | index of the opposite direction, `-1` when the reaction is not split |
| `is_bwd` | `True` for the backward (`_b`) reactions |
| `react_rev` | per reaction of `model`: `True` when it was split |
| `allGenes`, `rxnGeneMat`, `gpr` | genes, reaction × gene matrix (`n` × genes), GPR rules of the original reactions |
| `rxn_id_orig`, `mod_id` | reaction ids and id of the original model |

Method `net_flux(x)`: fluxes of the irreversible reactions (the first `n` values of `x`) →
net flux (forward − backward) of each reaction of the original model.

For a reaction `r` with bounds `(lb, ub)`:

| | lower bound | upper bound |
|---|---|---|
| `r_f` | `max(0, lb)` | `max(0, ub)` |
| `r_b` | `max(0, -ub)` | `max(0, -lb)` |

## Details

Bounds and objective coefficients can be changed before solving. A reaction that was not split
cannot become reversible later. To use a reaction in both directions, set its lower bound < 0
before you call `mod2irrev`, or give it in `rev`.

## Examples

```python
import cobra
from ccfbapy import mod2irrev

model = cobra.io.read_sbml_model('models/iAF1260.xml')
mod2 = mod2irrev(model, rev=['EX_glyc_e'])     # split glycerol exchange too (lower bound is 0)
print(mod2)                                   # <IrrevModel iAF1260: ... reactions (... split), 1668 metabolites>

i = mod2.react_id.index('EX_glyc_e_b')        # glycerol uptake
mod2.uppbnd[i] = 10
mod2.uppbnd[mod2.react_id.index('EX_glc__D_e_b')] = 0   # no glucose uptake
```

## See also

[`getccFBA_mat`](getccFBA_mat.md), [`getRevFlux`](getRevFlux.md)
