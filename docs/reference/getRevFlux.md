# getRevFlux: fluxes of the reversible model

*R help page: `man/getRevFlux.Rd`. R source: `R/sybilccFBA/R/getRevFlux.R`.*

!!! note "Not needed in Python"
    The Python functions already return the net fluxes of the reactions of the original model
    (`rxnFlx`). For a flux vector of an [irreversible model](mod2irrev.md), use
    `IrrevModel.net_flux(x)`.

## Description

Gets the fluxes of the reversible model from the fluxes of the irreversible model.

## Usage (R)

```r
getRevFlux(model, irrxns, fdirrev)
```

| Argument | Description |
|---|---|
| `model` | metabolic model of class `modelorg` |
| `irrxns` | ids of the reactions of the irreversible model |
| `fdirrev` | fluxes of the irreversible model |

## Value

A data frame with the reaction id of `model` (`rxn`), the forward (`fwd`) and backward (`bwd`)
flux. The net flux is `fwd - bwd`.

## Python equivalent

```python
import cobra
import numpy as np
import pandas as pd
from ccfbapy import getccFBA_mat, solve_ccFBA_mat

model = cobra.io.read_sbml_model('models/iAF1260.xml')
mw = pd.read_csv('data/Ecoli_MW.csv')
kcat = pd.read_csv('data/Ecoli_kcat.csv')
mat = getccFBA_mat(model, Kcat=kcat, MW=mw, RHS=0.27, medval=3600*22.6, verboseMode=0)
Res = solve_ccFBA_mat(mat, verboseMode=0)

mod2 = mat['model']
flx_ir = Res['flx_ir']                     # fluxes of the irreversible model (R: fdirrev)
net = mod2.net_flux(flx_ir)                # fwd - bwd, in the order of model.reactions
fluxes = pd.DataFrame({'rxn': mod2.rxn_id_orig,
                       'fwd': flx_ir[mod2.rev2irrev[:, 0]],
                       'bwd': np.where(mod2.react_rev, flx_ir[mod2.rev2irrev[:, 1]], 0)})
print(fluxes[fluxes['fwd'] + fluxes['bwd'] > 1e-6].head())
```

`Res['rxnFlx']` already holds `dict(zip(mod2.rxn_id_orig, net))`.

## See also

[`mod2irrev`](mod2irrev.md), [`getccFBA_sol`](getccFBA_sol.md) (`metflx`, `metflx_ir`)
