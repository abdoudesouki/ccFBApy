# 2/12/2018
# by Abdelmoneim Amer Desouki
# part of ccFBApy
# Python version of simulate_EColi (sybilccFBA, R/simulate_EColi.R)
import copy

import numpy as np
import pandas as pd

from getccFBA_mat import getccFBA_mat, solve_ccFBA_mat
from getccFBA_sol import flx_MC_tables


def simulate_EColi(model, mod2, kcat, mw, budget_C, CSList, atpmz=True, trns_rxns=None,
                   cpx_stoich=None, solver=None, verboseMode=1):
    """Simulate growth on a list of carbon sources with ccFBA (getccFBA_mat). The matrix is built once;
    for each carbon source the uptake of all carbon sources is closed, except the current one.

    model:     cobra Model
    mod2:      irreversible model (getccFBA_mat.mod2irrev(model, rev=...)): the uptake reactions in CSList
               must exist in it (i.e. the exchange reactions of the carbon sources must be split).
               It is not modified (a copy is used).
    kcat:      DataFrame with columns 'rxn_id', 'val' (1/s), 'dirxn', optionally 'src'
    mw:        DataFrame, first column gene, second column molecular weight (g/mmol)
    budget_C:  the budget C, for E. coli 0.27
    CSList:    carbon sources to be simulated: ids of reactions in mod2 (e.g. 'EX_glc__D(e)_b')
    atpmz:     set the ATPM reaction to zero, default True
    trns_rxns: transport reactions to be excluded from the budget (may contain other reactions if
               necessary): boolean sequence in the order of model.reactions, or list of reaction ids
    cpx_stoich: complex stoichiometry, DataFrame with columns 'genes' and 'stoich' (see getGpr1iso)
    solver:    accepted for compatibility with R, not used (cobra's default solver is used)

    Returns a dictionary:
        gr_data:    rows of all_flx of the objective reaction (growth rate in column flx) of all carbon sources
        uptake:     rows of all_flx of the uptake reaction of each carbon source
        all_flx:    cs, csname, rxn_id, flux (flx), kcat, kcat_src, enzyme concentration (enzConc),
                    obj, ub, ubR for all reactions of mod2 and all carbon sources
        rg_MC:      reaction gene molecular crowding (MC): all_flx_MC sorted in decreasing order of rxnMC
                    (rxnMCorder) for each carbon source; rgMCorder is the order by rgMC
        all_flx_MC: reaction-gene pairs (mrMat['rxnGeneCol']) with flx, grConc (concentration of the gene
                    used by the reaction), enzConc (concentration of the enzyme of the reaction, NaN for genes
                    of isoenzymes that were not chosen), geneConc, rgMC = gMW*grConc: MC of the gene part,
                    rxnMC = MW of enzyme * enzConc: MC of the reaction, not sorted
        mrMat:      the matrix from getccFBA_mat (not returned in R)
    cs is the number of the carbon source in CSList, starting at 1 as in R.
    """
    mod2 = copy.deepcopy(mod2)
    if atpmz:
        for i, rid in enumerate(mod2.react_id):
            if rid in ('ATPM', 'ATPM_f', 'ATPM_b'):
                mod2.lowbnd[i] = 0
                mod2.uppbnd[i] = 0
    bad = [c for c in CSList if c not in mod2.react_id]
    if bad:
        raise ValueError("CSList: reactions not in mod2: %s" % ', '.join(bad))
    nr = mod2.react_num
    medianVal = float(np.median(kcat['val'].astype(float)))
    mod2R = copy.deepcopy(mod2)
    react_id = [r.id for r in model.reactions]
    has_gpr = np.array([str(r.gpr).strip() != "" for r in model.reactions])
    if trns_rxns is not None:
        if len(trns_rxns) == len(react_id) and all(isinstance(x, (bool, np.bool_)) for x in trns_rxns):
            trns = np.asarray(trns_rxns, dtype=bool)
        else:
            trns = np.isin(react_id, [getattr(r, 'id', r) for r in trns_rxns])
        selected_rxns = [r for r, g, t in zip(react_id, has_gpr, trns) if g and not t]
        print('excluded transport reactions:%d' % trns.sum())
    else:
        selected_rxns = [r for r, g in zip(react_id, has_gpr) if g]

    mrMat = getccFBA_mat(model, mod2, kcat, MW=mw, verboseMode=2 if verboseMode > 0 else 0, RHS=budget_C,
                         cpx_stoich=cpx_stoich, medval=3600 * medianVal, selected_rxns=selected_rxns)
    rgm = mrMat['rxnGeneCol']
    # enzyme column of each flux column (first kcat constraint of the reaction)
    rec = mrMat['rxnExprCol'].drop_duplicates('flxCol').set_index('flxCol')['exprCol']
    flx_expr = rec.reindex(np.arange(nr))
    cs_ix = [mod2.react_id.index(c) for c in CSList]

    def take(x, cols, na=np.nan):
        cols = pd.array(cols, dtype='Int64')
        out = np.full(len(cols), na, dtype=float)
        ok = ~np.asarray(cols.isna())
        out[ok] = x[np.asarray(cols[ok], dtype=int)]
        return out

    all_flx, all_flx_MC, all_rg_MC = [], [], []
    prob = None
    for k, csname in enumerate(CSList):
        cs = k + 1
        print('%d %s' % (cs, csname))
        mod2.lowbnd[:] = mod2R.lowbnd
        mod2.uppbnd[:] = mod2R.uppbnd
        mod2.uppbnd[cs_ix] = 0
        mod2.uppbnd[cs_ix[k]] = 1000
        sol = solve_ccFBA_mat(mrMat, mod2=mod2, verboseMode=0, prob=prob)
        prob = sol['prob']
        print(sol['status'])
        x = sol['x']
        ###
        all_flx.append(pd.DataFrame({'cs': cs, 'csname': csname, 'rxn_id': mod2.react_id, 'flx': x[:nr],
                                     'kcat': mrMat['kcat_ir'], 'kcat_src': mrMat['kcat_src_ir'],
                                     'enzConc': take(x, flx_expr.values), 'obj': mod2.obj_coef,
                                     'ub': mod2.uppbnd.copy(), 'ubR': mod2R.uppbnd.copy()}))
        flx_MC, rg_MC = flx_MC_tables(mrMat, x)
        flx_MC.insert(0, 'csname', csname)
        flx_MC.insert(0, 'cs', cs)
        rg_MC.insert(0, 'csname', csname)
        rg_MC.insert(0, 'cs', cs)
        all_flx_MC.append(flx_MC)
        all_rg_MC.append(rg_MC)

    all_flx = pd.concat(all_flx, ignore_index=True)
    upt = all_flx[all_flx['csname'] == all_flx['rxn_id']]
    gr_data = all_flx[all_flx['obj'] != 0]
    return {'gr_data': gr_data, 'uptake': upt, 'all_flx': all_flx,
            'rg_MC': pd.concat(all_rg_MC, ignore_index=True), 'all_flx_MC': pd.concat(all_flx_MC, ignore_index=True),
            'mrMat': mrMat}
