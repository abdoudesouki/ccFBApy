# 14/4/2021
# implementation of over expression of one gene using ccFBA (Tin Y Pang).
# By Abdelmonem Amer Desouki
# part of ccFBApy, python version of simGeneOverExpr.R
import inspect

import numpy as np
import pandas as pd


def simGeneOverExpr(gx, y_pct, model, mod2=None, Kcat=None, MW=None, selected_rxns=None, verboseMode=2,
                    objVal=None, RHS=None, solver=None, medval=None, FUN=None):
    """Over expression of one gene with ccFBA: the gene gx must use at least the fraction y_pct of the budget C (RHS).

    gx:     the gene to be over expressed (e.g. 'b2323')
    y_pct:  fraction of C devoted to gene gx (at most 1), e.g. 0.05 for 5%
    FUN:    ccFBA function: cfba_moment (default), ccFBA_mr or getccFBA_sol (python versions of cfba_moment,
            cfba_moment_mr, getccFBA_sol). It must return a dictionary with obj, status, rxnFlx (net flux of each
            reaction of model) and geneConc (gene -> concentration).
    model, mod2, Kcat, MW, selected_rxns, verboseMode, RHS, solver, medval: passed to FUN (mod2 only when FUN
            has this parameter, i.e. getccFBA_sol). RHS (the budget C) is required. objVal is accepted for
            compatibility with R and not used.
    The fluxes are computed with the bounds of model (cfba_moment, ccFBA_mr) or mod2 (getccFBA_sol, default
    mod2irrev(model)); the minimum total flux step uses the bounds of model.

    Steps (as in R):
    1. WTSol: FUN without over expression. MC of gx = concentration * MW.
    2. If gx already uses at least y_pct of C, OExSol = WTSol (OverExprEff False). Otherwise FUN is called again
       with the lower bound of the concentration of gx = RHS*y_pct/MW (OverExprEff True).
       When gx has no concentration variable in the problem of FUN (getccFBA_sol: the gene is not in the chosen
       isoenzyme of any reaction; R stops with an error in this case) the over expressed protein cannot catalyze any
       reaction, so it only takes its part of the budget: FUN is called with RHS*(1-y_pct).
    3. Minimum total flux (MTF, sybil's algorithm "mtf"): the fluxes of the reactions that have at least one gene
       with a (non zero) molecular weight, and of the objective reaction, are fixed to the fluxes of OExSol; the sum
       of the absolute fluxes of all reactions is minimized.

    Returns a dictionary:
        Sol:         minimum total flux solution: status, obj (sum of absolute fluxes), rxnFlx (reaction -> flux)
        OExSol:      solution returned from FUN with over expression of gene gx
        WTSol:       wild type solution before over expression
        OverExprEff: True when the over expression changed the problem (gx used less than y_pct of C)
        objVal:      flux of the objective reaction of OExSol (dictionary when there are several)
        gxMC_WT, gxMC_OEx: MC (g/gDW) of gene gx in WTSol and OExSol (not returned in R)
    """
    if FUN is None:
        from .cfba_moment import cfba_moment
        FUN = cfba_moment
    allGenes = [g.id for g in model.genes]
    if gx not in allGenes:
        raise ValueError('simGeneOverExpr: Error: gene %s is not in the model' % gx)
    if y_pct > 1:
        raise ValueError('simGeneOverExpr: Error: y_pct = %.3f more than 1' % y_pct)
    if RHS is None:
        raise ValueError('simGeneOverExpr: Error: RHS (the budget C) is required')
    mw_gx = MW[MW.iloc[:, 0] == gx].iloc[:, 1]
    if len(mw_gx) == 0 or pd.isna(mw_gx.iloc[0]) or float(mw_gx.iloc[0]) == 0:
        raise ValueError('simGeneOverExpr: Error: gene %s has no molecular weight in MW' % gx)
    mw_gx = float(mw_gx.iloc[0])

    args = dict(model=model, Kcat=Kcat, MW=MW, verboseMode=verboseMode, RHS=RHS, medval=medval)
    params = inspect.signature(FUN).parameters
    if 'mod2' in params and mod2 is not None:
        args['mod2'] = mod2
    if solver is not None:
        args['solver'] = solver
    if selected_rxns is not None:
        args['selected_rxns'] = selected_rxns

    ####
    WTSol = FUN(**args)
    if WTSol['status'] != 'optimal':
        print('simGeneOverExpr: Warning: wild type solution status: %s' % WTSol['status'])
    has_var = gx in WTSol['geneConc']
    gxMC = (WTSol['geneConc'][gx] if has_var else 0.0) * mw_gx
    if verboseMode > 0:
        print('simGeneOverExpr: MC of gene %s =%.9f , %.3f %% of C by default' % (gx, gxMC, 100 * gxMC / RHS))
    if gxMC / RHS >= y_pct:
        OExSol = WTSol
        OverExprEff = False
        gxMC_OEx = gxMC
    else:
        ##optimal Solution contains less than y_pct of MC of gene X
        if has_var:
            OExSol = FUN(**args, geneConcLB={gx: (RHS * y_pct) / mw_gx})  # Only Lower bound Update 20/4/2021
            gxMC_OEx = OExSol['geneConc'][gx] * mw_gx
        else:
            if verboseMode > 0:
                print('simGeneOverExpr: gene %s has no concentration variable (it is not used in any enzyme), '
                      'its MC is taken from the budget: RHS = %g' % (gx, RHS * (1 - y_pct)))
            OExSol = FUN(**dict(args, RHS=RHS * (1 - y_pct)))
            gxMC_OEx = RHS * y_pct
        if verboseMode > 0:
            print('simGeneOverExpr: MC of gene %s =%.6f , %.3f %% of C after over expression'
                  % (gx, gxMC_OEx, 100 * gxMC_OEx / RHS))
        OverExprEff = True
    if OExSol['status'] != 'optimal':
        print('simGeneOverExpr: Warning: over expression solution status: %s' % OExSol['status'])

    ####### MinTotal flux : Hint: from ATMFBA  ######
    ###-----------------MTF step over nonmeasured expression ----------
    if verboseMode > 0:
        print('simGeneOverExpr: calculate min total flux')
    mwv = pd.to_numeric(MW.iloc[:, 1], errors='coerce')
    msrdGenes = set(MW.iloc[:, 0][(mwv != 0) & mwv.notna()].astype(str))
    hasGpr = [r for r in model.reactions if any(g.id in msrdGenes for g in r.genes)]  # GPR with a measured gene
    if verboseMode > 0:
        print('simGeneOverExpr: min total flux #Rxns:%d, to be affected=%d ' % (len(model.reactions),
                                                                               len(model.reactions) - len(hasGpr)))
    OExFlux = OExSol['rxnFlx']
    objRxns = [r for r in model.reactions if r.objective_coefficient != 0]
    with model:  # changes are reverted at the end
        for r in set(hasGpr) | set(objRxns):
            r.bounds = (OExFlux[r.id], OExFlux[r.id])
        model.objective = model.problem.Objective(0, direction='min')
        model.objective.set_linear_coefficients(
            {v: 1 for r in model.reactions for v in (r.forward_variable, r.reverse_variable)})
        mtf = model.optimize(raise_error=False)
        Sol = {'status': mtf.status, 'obj': mtf.objective_value,
               'rxnFlx': mtf.fluxes.to_dict() if mtf.status == 'optimal' else None}
    if Sol['status'] != 'optimal':
        print('simGeneOverExpr: Warning: min total flux status: %s' % Sol['status'])
    objVal = {r.id: OExFlux[r.id] for r in objRxns}
    if len(objVal) == 1:
        objVal = list(objVal.values())[0]
    return {'Sol': Sol, 'OExSol': OExSol, 'WTSol': WTSol, 'OverExprEff': OverExprEff, 'objVal': objVal,
            'gxMC_WT': gxMC, 'gxMC_OEx': gxMC_OEx}
