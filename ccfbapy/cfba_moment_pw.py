# MOMENT method with pairwise OR, like the MATLAB implementation of MOMENT.
# by Abdelmoneim Amer Desouki
# part of ccFBApy, python version of cfba_moment_pw (sybilccFBA, R/cfba_moment_pw.R)
#
# vi : mmol/gDW/hr, rxn catalyzed by gi with Kcati
# Kcat: 1/s in the input, converted to 1/hr
# gi: mmol/gDW
# MWi: g/mmol
# C : g protein/gDW (metabolic)
import numpy as np
import pandas as pd
import scipy.sparse as sp

from .getccFBA_mat import mod2irrev, solve_ccFBA_mat, _build_problem
from .gpr1iso import _gpr_terms

SPONT = 's0001'  # spontaneous reactions (E. coli models)


def _kcat_lookup(Kcat):
    """(reaction id, dirxn) -> kcat (1/h); the first value is used when a reaction is repeated."""
    kmap = {}
    for rid, d, v in zip(Kcat['rxn_id'], Kcat['dirxn'], Kcat['val']):
        kmap.setdefault((str(rid), int(d)), float(v))
    return kmap


def cfba_moment_pw(model, mod2=None, Kcat=None, MW=None, selected_rxns=None, verboseMode=2, objVal=None,
                   RHS=None, solver=None, medval=None, pairwise=True, interface=None):
    """MOMENT with pairwise OR (like the MATLAB implementation of MOMENT).

    model:   cobra Model
    mod2:    irreversible model, mod2irrev(model) (default). The bounds of its reactions are used.
    Kcat:    DataFrame rxn_id, val (kcat in 1/s), dirxn (1 forward, -1 reverse); rxn_id with or without 'R_'
    MW:      DataFrame, first column gene, second column molecular weight (g/mmol). When None, no budget.
    selected_rxns: ids of reactions of mod2 to constrain (default: all reactions of mod2). An id of a reaction of
             model that was split selects both directions.
    verboseMode: 0 silent, 1 status, >2 each reaction
    objVal:  when given, the problem is to find the minimum budget that gives this objective value
    RHS:     the budget C (g/gDW), e.g. 0.27 for E. coli
    solver:  accepted for compatibility with R, not used (see interface)
    medval:  kcat (1/h) for reactions without kcat, default median of Kcat (in 1/h)
    pairwise: True (default): rules with more than two isoenzymes are written as a chain of pairwise sums, as in R.
             False: one row v <= kcat*(sum of the isoenzymes), as in cfba_moment (same solution).
    interface: optlang interface used to solve, default cobra's default solver

    The kcat of reaction '<id>_b' of mod2 is the value of '<id>' with dirxn -1, the others use dirxn 1.
    Reactions without GPR or with the GPR 's0001' are not constrained; the isoenzyme 's0001' is removed from OR rules.

    Returns a dict: obj, status, rxnFlx (net flux of the reactions of model), flx_ir (fluxes of mod2),
    geneConc (gene -> concentration), x (values of all columns), geneCol (DataFrame gene, Col),
    prob (optlang problem), mat (the constraint matrix, can be solved again with solve_ccFBA_mat after
    changing the bounds of mat['model']), sol (result of solve_ccFBA_mat).
    """
    if Kcat is None:
        raise ValueError('Kcat is required')
    if mod2 is None:
        mod2 = mod2irrev(model)
    n, m = mod2.react_num, mod2.met_num
    Kcat = Kcat.copy()
    Kcat['val'] = Kcat['val'].astype(float) * 60 * 60  # convert to /hr
    if MW is not None:
        MW = MW.copy()
        MW[MW.columns[1]] = MW.iloc[:, 1].astype(float)
    if medval is None:
        medval = float(np.median(Kcat['val']))
    kmap = _kcat_lookup(Kcat)

    ix = {r: i for i, r in enumerate(mod2.react_id)}
    if selected_rxns is None or len(selected_rxns) == 0:
        sel = list(range(n))
    else:
        sel = []
        orig = {r: i for i, r in enumerate(mod2.rxn_id_orig)}
        for r in selected_rxns:
            r = getattr(r, 'id', r)
            if r in ix:
                sel.append(ix[r])
            elif r in orig:  # reaction of model: both directions
                sel += sorted(set(mod2.rev2irrev[orig[r]]))
            else:
                raise ValueError("selected_rxns: '%s' is not a reaction of mod2" % r)
        sel = list(dict.fromkeys(sel))

    # columns: fluxes of mod2, genes, auxiliary variables
    cnames = ['R_' + x for x in mod2.react_id]
    genes = [g for g in mod2.allGenes if g != SPONT]
    geneCol = pd.DataFrame({'gene': genes, 'Col': n + np.arange(len(genes))})
    gcol = dict(zip(geneCol['gene'], geneCol['Col']))
    cnames += ['g_' + g for g in genes]
    colid = n + len(genes)
    rnames = ['S_' + x for x in mod2.met_id]
    rowind = m
    entries = []  # (row, col, value)
    rtype, rlb, rub = ['E'] * m, [0.0] * m, [0.0] * m

    def add_col(name):
        nonlocal colid
        cnames.append(name)
        colid += 1
        return colid - 1

    def add_row(name, cols, vals, typ='U', lb=0.0, ub=0.0):
        nonlocal rowind
        entries.extend((rowind, c, v) for c, v in zip(cols, vals))
        rnames.append(name)
        rtype.append(typ)
        rlb.append(lb)
        rub.append(ub)
        rowind += 1
        return rowind - 1

    rxnkcatRow = []
    for vind in sel:
        r = mod2.react_id[vind]
        r0 = mod2.rxn_id_orig[mod2.irrev2rev[vind]]
        dirxn = -1 if mod2.is_bwd[vind] else 1
        key = (r0, dirxn) if (r0, dirxn) in kmap else ('R_' + r0, dirxn)
        kval, ksrc = (kmap[key], 'Kcat') if key in kmap else (medval, 'medianVal')

        terms = [t for t in _gpr_terms(model.reactions[int(mod2.irrev2rev[vind])].gpr)]
        # remove the spontaneous 'gene' from the rules (R does this for the rules of iAF1260 and iJO1366)
        if len(terms) > 1:
            terms = [t for t in terms if t != [SPONT]]
        terms = [[g for g in t if g != SPONT] for t in terms]
        terms = [t for t in terms if t]
        if not terms:  # no GPR or only s0001
            continue
        if verboseMode > 2:
            print(r, terms)
        if len(terms) == 1:  # no OR (only one term)
            if len(terms[0]) == 1:  # 1-one gene
                ri = add_row('rkcat_' + r, [vind, gcol[terms[0][0]]], [1.0, -kval])
            else:  # 2-straight ANDs: aux <= each gene, v <= kcat*aux
                aux = add_col('cpx_' + r)
                ri = add_row('rkcat_' + r, [vind, aux], [1.0, -kval])
                for g in terms[0]:
                    add_row('min_%s_%s' % (g, r), [aux, gcol[g]], [1.0, -1.0])
            rxnkcatRow.append((vind, r, dirxn, kval, ksrc, ri, 1))
            continue
        # 3-OR / 4-mixed OR and AND: one column per term (the gene, or an auxiliary variable for an AND term)
        tcols = []
        for p, t in enumerate(terms):
            if len(t) == 1:
                tcols.append(gcol[t[0]])
            else:
                aux = add_col('cpx_%s_%d' % (r, p))
                for g in t:
                    add_row('min_%s_%s_%d' % (g, r, p), [aux, gcol[g]], [1.0, -1.0])
                tcols.append(aux)
        tcols = sorted(set(tcols))  # R: which(row_vals != 0), i.e. in the order of the columns
        if len(tcols) <= 2 or not pairwise:
            ri = add_row('rkcat_' + r, [vind] + tcols, [1.0] + [-kval] * len(tcols))
        else:  # pairwise: aux1 <= t1 + t2, aux_k <= t_(k+1) + aux_(k-1), v <= kcat*(t_last + aux)
            o_aux = add_col('auxOR_%s_1' % r)
            add_row('AuxOR_%s_1' % r, [tcols[0], tcols[1], o_aux], [-1.0, -1.0, 1.0])
            for k, c in enumerate(tcols[2:-1], start=2):
                n_aux = add_col('auxOR_%s_%d' % (r, k))
                add_row('AuxOR_%s_%d' % (r, k), [c, o_aux, n_aux], [-1.0, -1.0, 1.0])
                o_aux = n_aux
            ri = add_row('rkcat_' + r, [tcols[-1], o_aux, vind], [-kval, -kval, 1.0])
        rxnkcatRow.append((vind, r, dirxn, kval, ksrc, ri, len(tcols)))

    ncol = colid
    clb = np.zeros(ncol)
    cub = np.full(ncol, 1000.0)
    clb[:n] = mod2.lowbnd
    cub[:n] = mod2.uppbnd
    obj_cf = np.zeros(ncol)
    obj_cf[:n] = mod2.obj_coef
    obj_dir = mod2.obj_dir
    cc_ind = mod_objc_ind = None
    if MW is not None:  # Sum(gi * MWi) <= C [g protein/gDW]
        mw_map = {}
        for g, v in zip(MW.iloc[:, 0], MW.iloc[:, 1]):
            mw_map.setdefault(str(g).strip(), float(v))
        has = [(gcol[g], mw_map[g]) for g in genes if g in mw_map and not np.isnan(mw_map[g])]  # genes with MW
        lcind = [c for c, _ in has]
        lnz = [v for _, v in has]
        if objVal is not None:  # minimum budget that gives objVal
            oc = np.nonzero(mod2.obj_coef)[0]
            mod_objc_ind = add_row('ObjC_', list(oc), list(mod2.obj_coef[oc]), typ='L', lb=float(objVal), ub=np.inf)
            obj_cf[:] = 0
            obj_cf[lcind] = lnz
            obj_dir = 'min'
        elif RHS is not None:  # Add solvent constraint
            cc_ind = add_row('MW', lcind, lnz, typ='U', lb=0.0, ub=float(RHS))
        elif verboseMode > 0:
            print('cfba_moment_pw: Warning: RHS and objVal are None, no budget constraint')

    S = mod2.S.tocoo()
    er = np.array([e[0] for e in entries], dtype=int)
    ec = np.array([e[1] for e in entries], dtype=int)
    ev = np.array([e[2] for e in entries], dtype=float)
    LHS = sp.csr_matrix((np.concatenate([S.data, ev]), (np.concatenate([S.row, er]), np.concatenate([S.col, ec]))),
                        shape=(rowind, ncol))
    rxnkcatRow = pd.DataFrame(rxnkcatRow, columns=['vind', 'rxn_id', 'dirxn', 'kval', 'kcat_src', 'rowind', 'nterms'])
    mat = {'LHS': LHS, 'rlb': np.array(rlb), 'rub': np.array(rub), 'rtype': rtype, 'rnames': rnames,
           'clb': clb, 'cub': cub, 'cnames': cnames, 'obj_cf': obj_cf, 'obj_dir': obj_dir,
           'geneCol': geneCol, 'rxnkcatRow': rxnkcatRow, 'model': mod2, 'cc_ind': cc_ind,
           'mod_objc_ind': mod_objc_ind}
    if verboseMode > 0:
        print('cfba_moment_pw: LP %d rows, %d columns (%d genes), %d constrained reactions'
              % (rowind, ncol, len(genes), len(rxnkcatRow)))
    sol = solve_ccFBA_mat(mat, interface=interface, verboseMode=verboseMode)
    if verboseMode > 0:
        print('cfba_moment_pw: %s, objective=%s' % (sol['status'], sol['obj']))
    return {'obj': sol['obj'], 'status': sol['status'], 'rxnFlx': sol['rxnFlx'], 'flx_ir': sol['flx_ir'],
            'geneConc': sol['geneConc'], 'x': sol['x'], 'geneCol': geneCol, 'prob': sol['prob'], 'mat': mat,
            'sol': sol}
