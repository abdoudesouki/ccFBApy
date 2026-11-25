# ccFBA with gene allocation: isoenzymes and protein complexes without choosing one isoenzyme (getGpr1iso).
# part of ccFBApy
#
# Problem of MOMENT with rules like (A and B) or (A and C): each complex (term) t gets a variable with
# cpx_t <= g for each of its genes, separately, so the concentration of the shared subunit A is used by both
# complexes at the same time and the cost of the enzyme is under estimated. getGpr1iso avoids this by keeping
# only the cheapest isoenzyme of each reaction (a preprocessing step that removes the other isoenzymes).
#
# Here the concentration of each gene is allocated to the enzymes (terms) that use it:
#   e_rt : concentration of enzyme (term) t of reaction r                 [mmol/gDW]
#   v_r  <= kcat_r * sum_t e_rt          (forward and backward direction, with their own kcat)
#   sum_{(r,t) containing g} st_gt * e_rt = g_g   for each gene g (st_gt: number of subunits of g in t)
#   sum_g MW_g * g_g <= C
# A subunit shared by several complexes, or by several reactions (multi-functional enzyme), is counted once for
# each use, so the cost is never under estimated, and the LP chooses the isoenzymes (a mixture when it is
# better, e.g. when a gene is knocked out or limited).
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .getccFBA_mat import mod2irrev, solve_ccFBA_mat, _build_problem, _gene_bounds
from .gpr1iso import _gpr_terms, _parse_stoich


def _term_stoich(genes, st_map):
    """Subunits of each gene of a term (sorted genes) from cpx_stoich; 1 when not given."""
    st = st_map.get(','.join(genes))
    sv = None
    if st is not None:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            sv = _parse_stoich(st, len(genes), ','.join(genes))
    return [1.0] * len(genes) if sv is None else sv


def getccFBA_iso_mat(model, mod2=None, Kcat=None, MW=None, selected_rxns=None, verboseMode=1, objVal=None,
                     RHS=None, medval=None, geneConcLB=None, geneConcUB=None, cpx_stoich=None, C_mu_coef=0,
                     gene_alloc='E', rm_genes=None):
    """ccFBA problem with gene allocation, as a constraint matrix (solve it with solve_ccFBA_mat).

    Isoenzymes are not removed (no getGpr1iso): every term of a GPR rule (gene or complex) is an enzyme with its
    own concentration variable, and the concentration of each gene is shared between all the enzymes that contain
    it (of the same or of other reactions). This removes the under estimation of the cost of rules like
    (A and B) or (A and C), and also includes multi-functional enzymes and complex stoichiometry.

    model, mod2, Kcat, MW, selected_rxns, objVal, RHS, medval, geneConcLB, geneConcUB, cpx_stoich, C_mu_coef:
        as in getccFBA_mat (Kcat in 1/s, MW in g/mmol, medval in 1/h, RHS = budget C in g/gDW; RHS None: minimum
        budget that gives objVal, default the FBA optimum).
    gene_alloc: 'E' (default): each gene concentration equals the sum of its uses; 'U': the sum of its uses is at
        most the gene concentration (an over expressed protein may stay unused).
    rm_genes: genes removed from the OR rules (the terms that contain them are dropped when other terms remain),
        e.g. ['s0001'] as in R's cfba_moment. Default None: keep them (a gene without MW has no cost).

    Returns a dict in the format of getccFBA_mat (0-based indices): LHS, rlb, rub, rtype, rnames, clb, cub,
    cnames, obj_cf, obj_dir, model (mod2), fluxind, geneCol (gene, Col, MW, nenz), enzCol (one row per enzyme:
    rsn, rxn, term, genes, stoich, nterms, MW (of the enzyme, with stoichiometry), Col), rxnkcatRow (one row per
    constrained direction: vind, rxn_id, dirxn, kval (1/h), kcat_src, rowind), allocRow, kcat_ir, kcat_src_ir,
    geneMW, cc_ind, mod_objc_ind.
    Column names: R_<rxn> (fluxes of mod2), g_<gene>, enz_<rxn>_<t> (term t of the rule, 1-based).
    Row names: S_<met>, rkcat_<rxn of mod2>, galloc_<gene>, cc or ObjC_.
    """
    if Kcat is None:
        raise ValueError('Kcat is required')
    if gene_alloc not in ('E', 'U'):
        raise ValueError("gene_alloc must be 'E' or 'U'")
    Kcat = Kcat.copy()
    Kcat['val'] = Kcat['val'].astype(float) * 60 * 60  # 1/s -> 1/h
    if medval is None:
        medval = float(np.median(Kcat['val']))
    if mod2 is None:
        mod2 = mod2irrev(model)
    if mod2.rxn_id_orig != [r.id for r in model.reactions] or mod2.allGenes != [g.id for g in model.genes]:
        raise ValueError('mod2 does not correspond to model')
    n, m = mod2.react_num, mod2.met_num
    react_id = [r.id for r in model.reactions]
    rxn_ix = {r: i for i, r in enumerate(react_id)}
    allGenes = [g.id for g in model.genes]
    rm_genes = set(rm_genes or [])

    if selected_rxns is None or len(selected_rxns) == 0:
        selected_rxns = [r.id for r in model.reactions if str(r.gpr).strip() != '']
    selected_rxns = list(dict.fromkeys(getattr(r, 'id', r) for r in selected_rxns))
    bad = [r for r in selected_rxns if r not in rxn_ix]
    if bad:
        raise ValueError('selected_rxns not in model: %s' % ', '.join(bad[:10]))

    mw_map = {}
    if MW is not None:
        for g, v in zip(MW.iloc[:, 0], MW.iloc[:, 1]):
            mw_map.setdefault(str(g).strip(), float(v))
    st_map = {}
    if cpx_stoich is not None:
        for g, s in zip(cpx_stoich['genes'], cpx_stoich['stoich']):
            st_map.setdefault(str(g).strip().replace(' ', ''), s)

    # kcat of each direction (as getccFBA_mat: with or without 'R_', the last value of a repeated reaction)
    has_src = 'src' in Kcat.columns
    kmap = {}
    for _, row in Kcat.iterrows():
        kmap[(str(row['rxn_id']), int(row['dirxn']))] = (float(row['val']), row['src'] if has_src else 'Kcat')
    kcat_ir = np.full(n, medval)
    kcat_src_ir = np.array(['medianVal'] * n, dtype=object)
    for j in range(n):
        r0 = react_id[mod2.irrev2rev[j]]
        d = -1 if mod2.is_bwd[j] else 1
        key = (r0, d) if (r0, d) in kmap else ('R_' + r0, d)
        if key in kmap:
            kcat_ir[j], kcat_src_ir[j] = kmap[key]

    # enzymes: one per term of the rule of each selected reaction
    enz = []
    for r in selected_rxns:
        rsn = rxn_ix[r]
        terms = _gpr_terms(model.reactions[rsn].gpr)
        if rm_genes and len(terms) > 1:
            kept = [t for t in terms if not rm_genes.intersection(t)]
            terms = kept if kept else terms
        if not terms:
            raise ValueError('selected reaction without GPR: %s' % r)
        for k, t in enumerate(terms):
            sv = _term_stoich(t, st_map)
            mw_t = sum(mw_map.get(g, 0.0) * s for g, s in zip(t, sv) if not np.isnan(mw_map.get(g, 0.0)))
            enz.append({'rsn': rsn, 'rxn': r, 'term': k + 1, 'genes': ','.join(t), 'stoich': ','.join(
                '%g' % s for s in sv), 'nterms': len(terms), 'MW': mw_t, '_g': t, '_s': sv})
    enzCol = pd.DataFrame(enz)
    used = sorted({g for e in enz for g in e['_g']}, key=allGenes.index)
    geneCol = pd.DataFrame({'gene': used, 'Col': n + np.arange(len(used)),
                            'MW': [mw_map.get(g, np.nan) for g in used],
                            'g': [allGenes.index(g) for g in used]})
    gcol = dict(zip(geneCol['gene'], geneCol['Col']))
    ne = len(enzCol)
    enzCol['Col'] = n + len(used) + np.arange(ne)
    geneCol['nenz'] = [sum(g in e['_g'] for e in enz) for g in used]
    ncol = n + len(used) + ne
    cnames = ['R_' + x for x in mod2.react_id] + ['g_' + g for g in used] + \
             ['enz_%s_%d' % (r, t) for r, t in zip(enzCol['rxn'], enzCol['term'])]

    rows, cols, vals = [], [], []
    rnames = ['S_' + x for x in mod2.met_id]
    rtype, rlb, rub = ['E'] * m, [0.0] * m, [0.0] * m
    nrow = m

    def add_row(name, cidx, cval, typ, lb, ub):
        nonlocal nrow
        rows.extend([nrow] * len(cidx))
        cols.extend(cidx)
        vals.extend(cval)
        rnames.append(name)
        rtype.append(typ)
        rlb.append(lb)
        rub.append(ub)
        nrow += 1
        return nrow - 1

    # 1- kcat: v_j - kcat_j * sum_t e_rt <= 0, for each direction of the selected reactions
    enz_of = enzCol.groupby('rsn')['Col'].apply(list).to_dict()
    kr = []
    for r in selected_rxns:
        rsn = rxn_ix[r]
        for j in sorted(set(mod2.rev2irrev[rsn])):
            ec = enz_of[rsn]
            ri = add_row('rkcat_' + mod2.react_id[j], [j] + ec, [1.0] + [-kcat_ir[j]] * len(ec), 'U', 0.0, 0.0)
            kr.append((j, mod2.react_id[j], -1 if mod2.is_bwd[j] else 1, kcat_ir[j], kcat_src_ir[j], ri))
    rxnkcatRow = pd.DataFrame(kr, columns=['vind', 'rxn_id', 'dirxn', 'kval', 'kcat_src', 'rowind'])

    # 2- gene allocation: sum_(r,t) st_gt * e_rt - g = 0 (or <= 0)
    uses = {g: ([], []) for g in used}
    for e, c in zip(enz, enzCol['Col']):
        for g, s in zip(e['_g'], e['_s']):
            uses[g][0].append(c)
            uses[g][1].append(s)
    ar = []
    for g in used:
        c, s = uses[g]
        ri = add_row('galloc_' + g, c + [gcol[g]], s + [-1.0], gene_alloc, 0.0 if gene_alloc == 'E' else -np.inf, 0.0)
        ar.append((g, ri, len(c)))
    allocRow = pd.DataFrame(ar, columns=['gene', 'rowind', 'nuses'])

    # bounds
    clb = np.zeros(ncol)
    cub = np.full(ncol, 1000.0)
    clb[:n], cub[:n] = mod2.lowbnd, mod2.uppbnd
    LB = _gene_bounds(geneConcLB, allGenes, 0)
    UB = _gene_bounds(geneConcUB, allGenes, 1000)
    clb[geneCol['Col'].values] = LB[geneCol['g'].values]
    cub[geneCol['Col'].values] = UB[geneCol['g'].values]

    # 3- budget
    gMW = geneCol['MW'].values
    has = ~np.isnan(gMW)
    lcind = list(geneCol['Col'].values[has])
    lnz = list(gMW[has])
    if verboseMode > 0 and (~has).sum():
        print('Warning : some genes have no MW:%d (%s), they have no cost'
              % ((~has).sum(), ', '.join(geneCol['gene'].values[~has][:10])))
    obj_cf = np.zeros(ncol)
    oc = np.nonzero(mod2.obj_coef)[0]
    cc_ind = mod_objc_ind = None
    if RHS is not None:
        cidx, cval = lcind, lnz
        if C_mu_coef != 0:  # C = RHS + C_mu_coef * objective
            cidx = lcind + list(oc)
            cval = lnz + list(-C_mu_coef * mod2.obj_coef[oc])
        cc_ind = add_row('cc', cidx, cval, 'U', 0.0, float(RHS))
        obj_cf[:n] = mod2.obj_coef
        obj_dir = mod2.obj_dir
    else:  # minimum budget that gives objVal
        if objVal is None:
            fba = _build_problem(mod2.S, [0.0] * m, [0.0] * m, ['E'] * m, mod2.lowbnd, mod2.uppbnd,
                                 cnames[:n], rnames[:m], mod2.obj_coef, mod2.obj_dir)
            if fba.optimize() != 'optimal':
                raise RuntimeError('Failed to calculate the objective value (FBA)')
            objVal = fba.objective.value
            if verboseMode > 0:
                print('Objective value to be used: %f' % objVal)
        if mod2.obj_dir == 'max':
            mod_objc_ind = add_row('ObjC_', list(oc), list(mod2.obj_coef[oc]), 'L', float(objVal), np.inf)
        else:
            mod_objc_ind = add_row('ObjC_', list(oc), list(mod2.obj_coef[oc]), 'U', -np.inf, float(objVal))
        obj_cf[lcind] = lnz
        obj_dir = 'min'

    S = mod2.S.tocoo()
    LHS = sp.csr_matrix((np.concatenate([S.data, vals]),
                         (np.concatenate([S.row, rows]), np.concatenate([S.col, cols]))), shape=(nrow, ncol))
    if verboseMode > 0:
        print('ccFBA_iso: LP %d rows, %d columns: %d genes, %d enzymes (%d reactions, %d with isoenzymes)'
              % (nrow, ncol, len(used), ne, len(selected_rxns), int((enzCol.groupby('rsn').size() > 1).sum())))
    gm = np.full(len(allGenes), np.nan)
    for i, g in enumerate(allGenes):
        gm[i] = mw_map.get(g, np.nan)
    return {'LHS': LHS, 'rlb': np.array(rlb), 'rub': np.array(rub), 'rtype': rtype, 'rnames': rnames,
            'clb': clb, 'cub': cub, 'cnames': cnames, 'obj_cf': obj_cf, 'obj_dir': obj_dir, 'model': mod2,
            'fluxind': np.arange(n), 'geneCol': geneCol.drop(columns='g'),
            'enzCol': enzCol.drop(columns=['_g', '_s']), 'rxnkcatRow': rxnkcatRow, 'allocRow': allocRow,
            'kcat_ir': kcat_ir, 'kcat_src_ir': kcat_src_ir, 'geneMW': gm, 'cc_ind': cc_ind,
            'mod_objc_ind': mod_objc_ind}


def ccFBA_iso(model, mod2=None, Kcat=None, MW=None, selected_rxns=None, verboseMode=1, objVal=None, RHS=None,
              medval=None, geneConcLB=None, geneConcUB=None, cpx_stoich=None, C_mu_coef=0, gene_alloc='E',
              rm_genes=None, mat=None, interface=None):
    """Build (getccFBA_iso_mat) and solve the ccFBA problem with gene allocation.

    Arguments as getccFBA_iso_mat. mat: a matrix from getccFBA_iso_mat to be solved again (its flux bounds are
    taken from mod2, default mat['model']; geneConcLB/geneConcUB are applied to it). interface: optlang interface.

    Returns a dict: obj, status, rxnFlx (net flux of the reactions of model), flx_ir (fluxes of mod2), geneConc
    (gene -> concentration, mmol/gDW), enzConc (DataFrame enzCol with conc: concentration of each enzyme (term),
    MC = MW * conc (g/gDW) and flx: net flux of its reaction), rxnMC (Series: protein cost of each reaction,
    g/gDW; with gene_alloc 'E' the costs add up to the budget used), x, prob, mat.
    """
    if mat is None:
        mat = getccFBA_iso_mat(model, mod2=mod2, Kcat=Kcat, MW=MW, selected_rxns=selected_rxns,
                               verboseMode=verboseMode, objVal=objVal, RHS=RHS, medval=medval,
                               geneConcLB=geneConcLB, geneConcUB=geneConcUB, cpx_stoich=cpx_stoich,
                               C_mu_coef=C_mu_coef, gene_alloc=gene_alloc, rm_genes=rm_genes)
    elif geneConcLB is not None or geneConcUB is not None:
        mat = dict(mat)
        allGenes = [g.id for g in model.genes]
        gi = [allGenes.index(g) for g in mat['geneCol']['gene']]
        mat['clb'], mat['cub'] = mat['clb'].copy(), mat['cub'].copy()
        mat['clb'][mat['geneCol']['Col'].values] = _gene_bounds(geneConcLB, allGenes, 0)[gi]
        mat['cub'][mat['geneCol']['Col'].values] = _gene_bounds(geneConcUB, allGenes, 1000)[gi]
    mod2 = mat['model'] if mod2 is None else mod2
    sol = solve_ccFBA_mat(mat, mod2=mod2, interface=interface, verboseMode=verboseMode)
    x = sol['x']
    enzConc = mat['enzCol'].copy()
    enzConc['conc'] = x[enzConc['Col'].values]
    enzConc['MC'] = enzConc['MW'] * enzConc['conc']
    enzConc['flx'] = [sol['rxnFlx'][r] for r in enzConc['rxn']]
    rxnMC = enzConc.groupby('rxn', sort=False)['MC'].sum().sort_values(ascending=False)
    if verboseMode > 0:
        print('ccFBA_iso: %s, objective=%s' % (sol['status'], sol['obj']))
    return {'obj': sol['obj'], 'status': sol['status'], 'rxnFlx': sol['rxnFlx'], 'flx_ir': sol['flx_ir'],
            'geneConc': sol['geneConc'], 'enzConc': enzConc, 'rxnMC': rxnMC, 'x': x, 'prob': sol['prob'],
            'mat': mat}
