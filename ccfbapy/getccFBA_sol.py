# 08/05/2021
# by Abdelmoneim Amer Desouki
# part of ccFBApy
# Python version of getccFBA_sol (sybilccFBA, R/getccFBA_sol.R): gets fluxes from the ccFBA matrix
import numpy as np
import pandas as pd

from .getccFBA_mat import getccFBA_mat, solve_ccFBA_mat, _gene_bounds


def flx_MC_tables(mrMat, x):
    """Protein cost (molecular crowding, MC) of reaction-gene pairs from a solution x of the matrix mrMat.
    Returns (flx_MC, rg_MC): flx_MC has the columns of mrMat['rxnGeneCol'] and flx, grConc (concentration of the
    gene part used by the reaction), enzConc (enzyme concentration of the reaction, NaN for genes of isoenzymes that
    were not chosen), geneConc, rgMC = gMW*grConc, rxnMC = MW of the enzyme * enzConc; rg_MC is flx_MC sorted in
    decreasing order of rxnMC, with rgMCorder (order by rgMC) and rxnMCorder.
    The enzyme concentration is x[exprCol] / r1gstoich as in the kcat constraint (R divides by the stoichiometry of
    each gene, which gives Inf for genes of isoenzymes that were not chosen and a wrong cost for complexes with
    stoichiometry)."""
    rgm = mrMat['rxnGeneCol']
    g1 = mrMat['gpr1iso'].set_index('rsn')
    tng = g1['termNG'].reindex(rgm['revRxn']).values
    st1 = g1['stoich'].reindex(rgm['revRxn']).values
    r1gstoich = np.array([np.nan if pd.isna(t) else (1.0 if t > 1 else float(s)) for t, s in zip(tng, st1)])
    r1gstoich[rgm['stoich'].astype(float).values == 0] = np.nan  # gene not in the chosen enzyme

    def take(cols, na=np.nan):
        cols = pd.array(cols, dtype='Int64')
        out = np.full(len(cols), na, dtype=float)
        ok = ~np.asarray(cols.isna())
        out[ok] = x[np.asarray(cols[ok], dtype=int)]
        return out

    flx_MC = rgm.copy()
    flx_MC['flx'] = x[rgm['rxnCol'].values]
    flx_MC['grConc'] = take(rgm['grCol'], na=0.0)
    flx_MC['enzConc'] = take(rgm['gpr.exprCol']) / r1gstoich
    flx_MC['geneConc'] = take(rgm['gCol'], na=0.0)
    flx_MC['rgMC'] = rgm['gMW'].values * flx_MC['grConc'].values
    flx_MC['rxnMC'] = rgm['gpr.MW'].values * flx_MC['enzConc'].values
    rg_MC = flx_MC.sort_values('rgMC', ascending=False, kind='stable', na_position='last')
    rg_MC.insert(0, 'rgMCorder', np.arange(1, len(rg_MC) + 1))
    rg_MC = rg_MC.sort_values('rxnMC', ascending=False, kind='stable', na_position='last')
    rg_MC.insert(0, 'rxnMCorder', np.arange(1, len(rg_MC) + 1))
    return flx_MC, rg_MC


def getccFBA_sol(model, mod2=None, Kcat=None, MW=None, selected_rxns=None, verboseMode=2, objVal=None,
                 RHS=None, solver=None, medval=None, geneConcLB=None, geneConcUB=None,
                 cpx_stoich=None, C_mu_coef=0, mrMat=None, CSList=None, cs=1):
    """Builds the ccFBA matrix (getccFBA_mat) and solves it; can take a matrix mrMat and solve it.

    Parameters as in getccFBA_mat. When mrMat is given, only the gene concentration bounds (geneConcLB/geneConcUB)
    are applied to it and the bounds of the fluxes are taken from mod2 (default mrMat['model']).
    CSList and cs are accepted for compatibility with R (they are only used in a table that R does not return).

    Returns a dictionary:
        sol:       result of solve_ccFBA_mat (obj, status, x: values of all columns, prob, ...)
        metflx:    net flux of the reactions of model (Series), metflx_ir: fluxes of mod2 (Series)
        geneConc:  gene -> concentration (genes with a concentration variable), as the other python functions
                   (R returns a data frame gene, geneconc)
        flx_MC, rg_MC: protein cost of reaction-gene pairs (see flx_MC_tables)
        mrMat:     the matrix, prob: the optlang problem
        obj, status, rxnFlx (dict of metflx): same keys as cfba_moment and ccFBA_mr
    """
    if mrMat is None:
        # as in R, objVal and C_mu_coef are not passed to getccFBA_mat
        mrMat = getccFBA_mat(model=model, mod2=mod2, Kcat=Kcat, MW=MW, verboseMode=verboseMode, RHS=RHS,
                             cpx_stoich=cpx_stoich, geneConcLB=geneConcLB, geneConcUB=geneConcUB, medval=medval,
                             selected_rxns=selected_rxns)
    else:  # apply geneConc bounds
        allGenes = [g.id for g in model.genes]
        gix = mrMat['geneCol']['g'].values
        mrMat = dict(mrMat)
        mrMat['clb'] = mrMat['clb'].copy()
        mrMat['cub'] = mrMat['cub'].copy()
        mrMat['clb'][mrMat['geneCol']['Col'].values] = _gene_bounds(geneConcLB, allGenes, 0)[gix]
        mrMat['cub'][mrMat['geneCol']['Col'].values] = _gene_bounds(geneConcUB, allGenes, 1000)[gix]
    mod2 = mrMat['model'] if mod2 is None else mod2

    ##Solve #########################
    sol = solve_ccFBA_mat(mrMat, mod2=mod2, verboseMode=verboseMode)
    if verboseMode > 0:
        print(sol['status'])
    x = sol['x']
    flx_MC, rg_MC = flx_MC_tables(mrMat, x)
    geneConc = dict(zip(mrMat['geneCol']['gene'], x[mrMat['geneCol']['Col'].values]))
    metflx_ir = pd.Series(x[:mod2.react_num], index=mod2.react_id)  # fluxes of reactions in irreversible model
    metflx = pd.Series(mod2.net_flux(x), index=mod2.rxn_id_orig)  # fluxes of model
    return {'sol': sol, 'metflx': metflx, 'metflx_ir': metflx_ir, 'geneConc': geneConc,
            'rg_MC': rg_MC, 'flx_MC': flx_MC, 'mrMat': mrMat, 'prob': sol['prob'],
            'obj': sol['obj'], 'status': sol['status'], 'rxnFlx': metflx.to_dict()}
