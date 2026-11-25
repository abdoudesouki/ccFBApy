
# 15/2/2024
# by Abdelmoneim Amer Desouki
# part of ccFBApy
# Python version of getccFBA_mat (sybilccFBA, R/getccFBA_mat.R)

#given  Kcat vector (in 1/Sec), and molecular weights  (Maximize biomass)
#1-Add gi vars (gene concentrations)
#2-add gpr constraints: add f/b for reversible rxns react_rev
#		vi/kcati <= gi
#       vi cat by g1 & g2 :  add aux variable then Aux <= all(gi), vi<= Kcat*Aux
#       vi cat by g1 | g2 :   OR is eliminated by gpr1iso, choose the smallest enzyme
#      complex gpr recursively
##... ----------
#3-	mr:multiple reactions (multiple function enzymes
# 		for genes catalyzing more n rxn(n>1): add gij: j 1:n  such that Sum(gij) <= gi
#4-add Enzimatic constraint Sum(MW*gi)<=C
#
# All indices (rows, columns, reactions, genes) are 0-based, while they are 1-based in R.

import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .gpr1iso import getGpr1iso


class IrrevModel:
    """Irreversible version of a cobra model (python equivalent of sybil's mod2irrev / modelorg_irrev).

    A reaction that can carry negative flux (lower_bound < 0) is split in a forward reaction
    '<id>_f' and a backward reaction '<id>_b' (placed next to each other), the others keep their id.
    rev: optional list of reaction ids that are split too although their lower bound is >= 0 (as
    setting react_rev(model)=TRUE before sybil's mod2irrev); their backward reaction gets the bounds
    (0, -lower_bound), e.g. (0, 0), and can be opened later with uppbnd.
    Attributes (n: number of reactions of the irreversible model):
        react_id, react_name, lowbnd, uppbnd, obj_coef (length n), S (sparse, metabolites x n), met_id
        irrev2rev: index of the original reaction of each reaction
        rev2irrev: array (reactions of model x 2): forward and backward index (both equal the
                   index of the reaction when it is not split)
        matchrev:  index of the opposite direction, -1 when the reaction is not split
        is_bwd:    True for the backward ('_b') reactions
        react_rev: per reaction of the original model: True when it was split
        allGenes, rxnGeneMat (sparse, n x genes), gpr (rules of the original reactions)
        obj_dir:   'max' or 'min' (objective direction of the model)
    Bounds and objective coefficients can be changed before solving (as with sybil's
    uppbnd(mod2)/lowbnd(mod2)), but reactions that were not split cannot become reversible.
    """

    def __init__(self, model, rev=None):
        rxns = model.reactions
        nr = len(rxns)
        self.mod_id = model.id
        self.met_id = [m.id for m in model.metabolites]
        met_ix = {m.id: i for i, m in enumerate(model.metabolites)}
        self.allGenes = [g.id for g in model.genes]
        gene_ix = {g: i for i, g in enumerate(self.allGenes)}
        self.react_rev = np.array([r.lower_bound < 0 for r in rxns])
        if rev is not None:
            rev = [getattr(r, 'id', r) for r in rev]
            ix = {r.id: i for i, r in enumerate(rxns)}
            bad = [r for r in rev if r not in ix]
            if bad:
                raise ValueError("rev: reactions not in model: %s" % ', '.join(bad[:10]))
            self.react_rev[[ix[r] for r in rev]] = True
        self.obj_dir = model.objective_direction

        react_id, react_name, irrev2rev, lb, ub, obj, is_bwd = [], [], [], [], [], [], []
        rev2irrev = np.zeros((nr, 2), dtype=int)
        rows, cols, vals = [], [], []
        grow, gcol = [], []
        for i, r in enumerate(rxns):
            c = r.objective_coefficient
            dirs = [(r.id + '_f', 1, max(0.0, r.lower_bound), max(0.0, r.upper_bound), False),
                    (r.id + '_b', -1, max(0.0, -r.upper_bound), max(0.0, -r.lower_bound), True)] \
                if self.react_rev[i] else [(r.id, 1, r.lower_bound, r.upper_bound, False)]
            for rid, sgn, l, u, b in dirs:
                k = len(react_id)
                react_id.append(rid)
                react_name.append(r.name)
                irrev2rev.append(i)
                lb.append(l)
                ub.append(u)
                obj.append(sgn * c)
                is_bwd.append(b)
                for met, s in r.metabolites.items():
                    rows.append(met_ix[met.id])
                    cols.append(k)
                    vals.append(sgn * s)
                for g in r.genes:
                    grow.append(k)
                    gcol.append(gene_ix[g.id])
            k = len(react_id)
            rev2irrev[i] = (k - 2, k - 1) if self.react_rev[i] else (k - 1, k - 1)

        n = len(react_id)
        self.react_id = react_id
        self.react_name = react_name
        self.irrev2rev = np.array(irrev2rev, dtype=int)
        self.rev2irrev = rev2irrev
        self.is_bwd = np.array(is_bwd)
        self.matchrev = np.full(n, -1, dtype=int)
        self.matchrev[rev2irrev[self.react_rev, 0]] = rev2irrev[self.react_rev, 1]
        self.matchrev[rev2irrev[self.react_rev, 1]] = rev2irrev[self.react_rev, 0]
        self.lowbnd = np.array(lb, dtype=float)
        self.uppbnd = np.array(ub, dtype=float)
        self.obj_coef = np.array(obj, dtype=float)
        self.S = sp.csr_matrix((vals, (rows, cols)), shape=(len(self.met_id), n))
        self.rxnGeneMat = sp.csr_matrix((np.ones(len(grow)), (grow, gcol)), shape=(n, len(self.allGenes)))
        self.gpr = [str(r.gpr) for r in rxns]
        self.rxn_id_orig = [r.id for r in rxns]

    @property
    def react_num(self):
        return len(self.react_id)

    @property
    def met_num(self):
        return len(self.met_id)

    def net_flux(self, x):
        """Fluxes of the irreversible reactions -> net flux (forward - backward) of the original reactions."""
        x = np.asarray(x)[:self.react_num]
        net = x[self.rev2irrev[:, 0]].copy()
        net[self.react_rev] -= x[self.rev2irrev[self.react_rev, 1]]
        return net

    def __repr__(self):
        return "<IrrevModel %s: %d reactions (%d split), %d metabolites>" % (
            self.mod_id, self.react_num, int(self.react_rev.sum()), self.met_num)


def mod2irrev(model, rev=None):
    """Irreversible model (like mod2irrev of sybil), see IrrevModel.
    rev: reaction ids to be split even when their lower bound is >= 0."""
    return IrrevModel(model, rev=rev)


def _gene_bounds(val, genes, default):
    """geneConcLB/geneConcUB: None, dict gene -> value, or a sequence in the order of model.genes."""
    out = np.full(len(genes), float(default))
    if val is None:
        return out
    if isinstance(val, dict):
        ix = {g: i for i, g in enumerate(genes)}
        for g, v in val.items():
            g = getattr(g, 'id', g)
            if g not in ix:
                raise ValueError("gene '%s' in gene concentration bounds is not in the model" % g)
            out[ix[g]] = v
        return out
    val = np.asarray(val, dtype=float)
    if len(val) != len(genes):
        raise ValueError("gene concentration bounds: expected %d values (one per gene), got %d"
                         % (len(genes), len(val)))
    return val.copy()


def _fmt(x):
    """number -> text as R's paste() (1 -> '1', 2.5 -> '2.5')."""
    return format(float(x), ".15g")


def _build_problem(LHS, rlb, rub, rtype, clb, cub, cnames, rnames, obj, obj_dir, interface=None):
    """optlang problem from a constraint matrix, row types as in sybil:
    'E': row == rlb, 'U': row <= rub, 'L': row >= rlb, 'D': rlb <= row <= rub, 'F': free."""
    if interface is None:
        import cobra
        interface = cobra.Configuration().solver
    prob = interface.Model()
    LHS = sp.csr_matrix(LHS)
    vars_ = [interface.Variable(cnames[j], lb=clb[j], ub=cub[j]) for j in range(LHS.shape[1])]
    prob.add(vars_)
    cons = []
    for i in range(LHS.shape[0]):
        t = rtype[i]
        lo = rlb[i] if t in ('E', 'L', 'D') else None
        up = rlb[i] if t == 'E' else (rub[i] if t in ('U', 'D') else None)
        if lo is not None and not np.isfinite(lo):
            lo = None
        if up is not None and not np.isfinite(up):
            up = None
        cons.append(interface.Constraint(0, lb=lo, ub=up, name=rnames[i]))
    prob.add(cons)
    prob.update()
    for i, c in enumerate(cons):
        s, e = LHS.indptr[i], LHS.indptr[i + 1]
        if e > s:
            c.set_linear_coefficients({vars_[j]: v for j, v in zip(LHS.indices[s:e], LHS.data[s:e])})
    obj = np.asarray(obj, dtype=float)
    prob.objective = interface.Objective(0, direction=obj_dir)
    prob.objective.set_linear_coefficients({vars_[j]: obj[j] for j in np.nonzero(obj)[0]})
    prob.update()
    return prob


def getccFBA_mat(model, mod2=None, Kcat=None, MW=None, selected_rxns=None, verboseMode=2, objVal=None,
                 RHS=None, solver=None, medval=None,
                 geneConcLB=None, geneConcUB=None, cpx_stoich=None, C_mu_coef=0):
    """ccFBA model as a constraint matrix (python version of getccFBA_mat of sybilccFBA).

    Improvements over MOMENT: 1- multi-functional enzymes: a gene catalyzing n > 1 reactions gets one
    variable per reaction and their sum equals the gene concentration. 2- only one enzyme (the
    cheapest isoenzyme, see getGpr1iso) is used for each reaction. 3- optionally protein complex
    stoichiometry (cpx_stoich).

    model:      cobra Model
    mod2:       irreversible model (mod2irrev(model)); can be given to save the time of building it.
                Reactions with lower_bound < 0 are split in '<id>_f' and '<id>_b'.
    Kcat:       DataFrame with columns 'rxn_id', 'val' (kcat in 1/s) and 'dirxn' (1 forward,
                -1 reverse), optionally 'src' (source of the value). rxn_id is the reaction id, with
                or without the SBML prefix 'R_'.
    MW:         DataFrame, first column gene id, second column molecular weight (g/mmol)
    selected_rxns: reaction ids to be constrained; default: all reactions with a GPR.
    RHS:        the budget C (g protein/gDW), e.g. 0.27 for E. coli. When RHS is None the problem is
                to find the minimum budget that gives the objective value objVal (when objVal is None,
                it is the optimum of FBA).
    medval:     kcat (1/h) for reactions without kcat, default: median of Kcat.
    geneConcLB, geneConcUB: gene concentration bounds (mmol/gDW), dict gene -> value or a sequence in
                the order of model.genes. Default 0 and 1000.
    cpx_stoich: DataFrame with columns 'genes' and 'stoich' (see getGpr1iso).
    C_mu_coef:  C as a linear function of the objective (growth): C = RHS + C_mu_coef * mu.
    solver:     accepted for compatibility with R, not used (the model's solver interface is used).

    Returns a dictionary (row/column indices are 0-based):
        LHS (scipy sparse, rows x columns), rlb, rub, rnames, rtype (rows), clb, cub, cnames (columns),
        obj_cf, obj_dir: objective; fluxind, geneind: columns of the fluxes and of the genes;
        rxnkcatRow, minGiRow, mfeRow: kcat, complex and multi-functional enzyme constraints;
        rgst, mod_cpx_mw, gpr1iso: output of getGpr1iso; rxnExprCol: flux column, enzyme column and
        row of each kcat constraint; rxnGeneCol: reaction-gene map used to get the enzyme costs;
        geneCol, mfeCol: gene variables; kcat_ir, kcat_src_ir: kcat of each irreversible reaction;
        geneMW; model (the irreversible model mod2); mod_objc_ind: row of the constraint on the
        objective (RHS None); cc_ind: row of the capacity constraint (RHS given).
    Column names: R_<rxn> (fluxes of mod2), g1_<gene>_1 (gene catalyzing one reaction),
        gmr_<gene>_<n> (gene catalyzing n reactions), gr_<gene>_<rxn> (part of a multi-functional
        gene used by rxn), cpx_<rxn> (complex of rxn).
    Row types (sybil): 'E' row == rlb, 'U' row <= rub, 'L' row >= rlb.
    Use solve_ccFBA_mat() to solve the problem.
    """
    #Multiple OR are not done pairwise but evaluated in one line in contrast to MOMENT original implementation.
    #   objVal: when set the minimum C that gives this growth, will be returned as objective
    #			units: mmol/gDW/h
    # vi : mmol/gDW/hr  , rxn catal by gi with Kcati
    # Kcat: hr-1
    # gi: mmol/gDW
    # MWi: g/mmol
    # C : g/gDW protein (metabolic)
    if Kcat is None:
        raise ValueError("Kcat is required")
    if MW is None:
        raise ValueError("MW is required (it is used to choose the isoenzymes)")
    Kcat = Kcat.copy()
    Kcat['val'] = Kcat['val'].astype(float) * 60 * 60  # convert to 1/hr
    MW = MW.copy()
    MW[MW.columns[1]] = MW.iloc[:, 1].astype(float)
    if mod2 is None:
        mod2 = mod2irrev(model)
    if mod2.rxn_id_orig != [r.id for r in model.reactions] or mod2.allGenes != [g.id for g in model.genes]:
        raise ValueError("mod2 does not correspond to model")
    if medval is None:
        medval = float(np.median(Kcat['val']))

    react_id = [r.id for r in model.reactions]
    rxn_ix = {r: i for i, r in enumerate(react_id)}
    allGenes = [g.id for g in model.genes]
    nR, nG = len(react_id), len(allGenes)
    if selected_rxns is None or len(selected_rxns) == 0:  # as in R: all reactions with a GPR
        selected_rxns = [r.id for r in model.reactions if str(r.gpr).strip() != ""]
    selected_rxns = list(dict.fromkeys(getattr(r, 'id', r) for r in selected_rxns))
    bad = [r for r in selected_rxns if r not in rxn_ix]
    if bad:
        raise ValueError("selected_rxns not in model: %s" % ', '.join(bad[:10]))
    bad = [r for r in selected_rxns if str(model.reactions[rxn_ix[r]].gpr).strip() == ""]
    if bad:
        raise ValueError("selected_rxns without GPR: %s" % ', '.join(bad[:10]))
    selflg = np.zeros(nR, dtype=bool)
    selflg[[rxn_ix[r] for r in selected_rxns]] = True

    if verboseMode > 0:
        print("Preprocessing isoenzymes ... ")
    with warnings.catch_warnings():  # genes without MW are reported below
        warnings.filterwarnings('ignore', message='.*have no molecular weight.*')
        r1i = getGpr1iso(model, MW=MW, cpx_stoich=cpx_stoich)
    gpr1iso = r1i['gpr1iso'].copy()
    mod_cpx_mw = r1i['mod_cpx_mw']
    rgst = r1i['rgst']
    if gpr1iso['stoich'].isna().any():
        raise ValueError("There are NA values in gpr1iso")

    ### start making constraint matrix ---------------------------------------------------
    if verboseMode > 0:
        print("Preparing LP ... ")
    n = mod2.react_num
    m = mod2.met_num
    rnames = ['S_' + x for x in mod2.met_id]
    rlb = [0.0] * m
    rub = [0.0] * m
    rtype = ['E'] * m
    cnames = ['R_' + x for x in mod2.react_id]
    entries = {}  # (row, col) -> value, added to S(mod2)
    nrow = m

    def setval(rws, cls, vals):
        for i, j, v in zip(np.broadcast_to(rws, np.shape(cls)), cls, np.broadcast_to(vals, np.shape(cls))):
            entries[(int(i), int(j))] = float(v)

    # add all variables once
    #Add variables for all genes : geneCol
    rgst = rgst.copy()
    rgst.loc[~selflg, :] = 0
    nz = rgst.values != 0
    genecnt = nz.sum(axis=0)
    selectedGenes = genecnt != 0
    ng = int(selectedGenes.sum())
    geneind = np.full(nG, -1, dtype=int)  # -1 (NA in R) for unused genes
    geneind[selectedGenes] = n + np.arange(ng)
    rg_ind = np.argwhere(nz.T)[:, ::-1]  # (rxn, gene) pairs in R's order (by gene, then reaction)
    gsel = np.nonzero(selectedGenes)[0]
    geneCol = pd.DataFrame({'gene': [allGenes[g] for g in gsel], 'rxn': None, 'Col': geneind[gsel],
                            'cnt': genecnt[gsel],
                            'vname': ['%s_%s_%d' % ('g1' if genecnt[g] == 1 else 'gmr', allGenes[g], genecnt[g])
                                      for g in gsel],
                            'r': pd.array([pd.NA] * len(gsel), dtype='Int64'), 'g': gsel})
    first_r = {}
    for r, g in rg_ind:
        first_r.setdefault(g, r)
    one = geneCol['cnt'] == 1
    geneCol.loc[one, 'r'] = [first_r[g] for g in geneCol.loc[one, 'g']]
    ncol = n + ng
    cnames += list(geneCol['vname'])

    rg_mfe = rg_ind[genecnt[rg_ind[:, 1]] > 1] if len(rg_ind) else rg_ind
    mfeCol = pd.DataFrame({'gene': [allGenes[g] for g in rg_mfe[:, 1]],
                           'rxn': [react_id[r] for r in rg_mfe[:, 0]],
                           'Col': ncol + np.arange(len(rg_mfe)), 'cnt': genecnt[rg_mfe[:, 1]],
                           'rsn': rg_mfe[:, 0], 'g': rg_mfe[:, 1]})
    mfeCol['vname'] = ['gr_%s_%s' % (g, r) for g, r in zip(mfeCol['gene'], mfeCol['rxn'])]
    ncol += len(mfeCol)
    cnames += list(mfeCol['vname'])

    rNG = nz.sum(axis=1)
    mgrflg = np.nonzero((rNG > 1) & selflg)[0]  # Add auxiliary variables to multiple gene rules
    gprCol = pd.DataFrame({'rsn': mgrflg, 'rxn': [react_id[r] for r in mgrflg],
                           'Col': ncol + np.arange(len(mgrflg)),
                           'vname': ['cpx_' + react_id[r] for r in mgrflg]})
    ncol += len(gprCol)
    cnames += list(gprCol['vname'])

    # column of the enzyme (expression) of each reaction: complex (cpx_), part of a
    # multi-functional gene (gr_) or gene (g1_). Only for selected reactions (NA for the others).
    gpr_col = dict(zip(gprCol['rsn'], gprCol['Col']))
    mfe_genes = set(mfeCol['gene'])
    mfe_col_rsn = dict(zip(mfeCol['rsn'][::-1], mfeCol['Col'][::-1]))  # first match by rsn
    gene_col = dict(zip(geneCol['gene'], geneCol['Col']))
    exprCol = []
    for rsn, genes in zip(gpr1iso['rsn'], gpr1iso['genes']):
        if not selflg[rsn]:
            exprCol.append(pd.NA)
        elif rsn in gpr_col:
            exprCol.append(gpr_col[rsn])
        elif genes in mfe_genes:
            exprCol.append(mfe_col_rsn.get(rsn, pd.NA))
        else:
            exprCol.append(gene_col.get(genes, pd.NA))
    gpr1iso['exprCol'] = pd.array(exprCol, dtype='Int64')
    if verboseMode > 1:
        print('Enzyme variables:', pd.Series([cnames[c][:2] for c in gpr1iso['exprCol'].dropna()]).value_counts().to_dict())

    ### ----------------------------------------
    # kcat of forward (pv) and reverse (nv) direction of each reaction of the model
    has_src = 'src' in Kcat.columns
    kcat_pv = np.full(nR, medval, dtype=float)
    kcat_nv = np.full(nR, medval, dtype=float)
    kcat_src_pv = np.array(['medianVal'] * nR, dtype=object)
    kcat_src_nv = kcat_src_pv.copy()
    kmap = {}
    for _, row in Kcat.iterrows():  # last value wins when a reaction is repeated (as in R)
        kmap[(str(row['rxn_id']), int(row['dirxn']))] = (row['val'], row['src'] if has_src else 'Kcat')
    for i, r in enumerate(react_id):
        for d, kv, ks in ((1, kcat_pv, kcat_src_pv), (-1, kcat_nv, kcat_src_nv)):
            key = (r, d) if (r, d) in kmap else ('R_' + r, d)
            if key in kmap:
                kv[i], ks[i] = kmap[key]
    kcat_ir = np.where(mod2.is_bwd, kcat_nv[mod2.irrev2rev], kcat_pv[mod2.irrev2rev])
    kcat_src_ir = np.where(mod2.is_bwd, kcat_src_nv[mod2.irrev2rev], kcat_src_pv[mod2.irrev2rev])

    ### Add constraint
    #AddRow Vind - Kcat*colind <= 0 {all selected rxns}
    ### ----------------------------------------
    g1 = gpr1iso.set_index('rsn')
    rxns_ind = np.array([rxn_ix[r] for r in selected_rxns], dtype=int)
    rev = mod2.react_rev[rxns_ind]
    rxnkcatRow = pd.DataFrame({'rsn': rxns_ind, 'rxn_id': selected_rxns,
                               'vind': mod2.rev2irrev[rxns_ind, 0],
                               'vind_b': pd.array(np.where(rev, mod2.rev2irrev[rxns_ind, 1], -1), dtype='Int64'),
                               'stoich': g1.loc[rxns_ind, 'stoich'].values,
                               'termNG': g1.loc[rxns_ind, 'termNG'].values,
                               'exprCol': g1.loc[rxns_ind, 'exprCol'].values})
    rxnkcatRow.loc[~rev, 'vind_b'] = pd.NA
    rxnkcatRow['kval'] = kcat_ir[rxnkcatRow['vind']]
    rxnkcatRow['r1gstoich'] = [1.0 if t > 1 else float(s) for t, s in zip(rxnkcatRow['termNG'], rxnkcatRow['stoich'])]
    rxnkcatRow['rkval'] = np.where(rev, kcat_ir[mod2.rev2irrev[rxns_ind, 1]], np.nan)
    rxnkcatRow['rowind'] = nrow + np.arange(len(rxnkcatRow))
    ## 1-rkcat_ constr
    consNames = ['rkcat_' + mod2.react_id[v] for v in rxnkcatRow['vind']]
    rxnkcatRow['consNames'] = consNames
    if rxnkcatRow['exprCol'].isna().any():
        raise RuntimeError("No enzyme column for: %s" % ', '.join(rxnkcatRow.loc[rxnkcatRow['exprCol'].isna(), 'rxn_id']))
    exc = rxnkcatRow['exprCol'].values.astype(int)
    setval(rxnkcatRow['rowind'].values, rxnkcatRow['vind'].values, 1.0)
    setval(rxnkcatRow['rowind'].values, exc, -1 * rxnkcatRow['kval'].values / rxnkcatRow['r1gstoich'].values)
    nrow += len(rxnkcatRow)
    rnames += consNames
    rlb += [0.0] * len(rxnkcatRow)
    rub += [0.0] * len(rxnkcatRow)
    rtype += ['U'] * len(rxnkcatRow)
    ###
    flg = rev
    nb = int(flg.sum())
    rxnkcatRow['consNames_b'] = None
    rxnkcatRow['rowind_b'] = pd.array([pd.NA] * len(rxnkcatRow), dtype='Int64')
    if nb > 0:
        vb = rxnkcatRow.loc[flg, 'vind_b'].values.astype(int)
        rxnkcatRow.loc[flg, 'consNames_b'] = ['rbkcat_' + mod2.react_id[v] for v in vb]
        rxnkcatRow.loc[flg, 'rowind_b'] = nrow + np.arange(nb)
        rb = rxnkcatRow.loc[flg, 'rowind_b'].values.astype(int)
        setval(rb, vb, 1.0)
        setval(rb, exc[flg], -1 * rxnkcatRow.loc[flg, 'rkval'].values / rxnkcatRow.loc[flg, 'r1gstoich'].values)
        nrow += nb
        rnames += list(rxnkcatRow.loc[flg, 'consNames_b'])
        rlb += [0.0] * nb
        rub += [0.0] * nb
        rtype += ['U'] * nb

    rxnExprCol = pd.concat([
        pd.DataFrame({'flxCol': rxnkcatRow['vind'].values, 'exprCol': exc, 'rowind': rxnkcatRow['rowind'].values}),
        pd.DataFrame({'flxCol': rxnkcatRow.loc[flg, 'vind_b'].values.astype(int), 'exprCol': exc[flg],
                      'rowind': rxnkcatRow.loc[flg, 'rowind_b'].values.astype(int)})], ignore_index=True)
    ###------------------------------------------------------------------------
    ## 2-min_ constr for complexes termNG > 1, for each gene: mapping and to minimum
    #     gene >= stoich * cpx   (i.e. cpx <= gene/stoich)
    if verboseMode >= 2:
        print('Adding constraints for mapping AND to min')
    cpx_rxns = set(np.nonzero((rNG > 1) & selflg)[0])
    rg_r1cpx = np.array([x for x in rg_ind if x[0] in cpx_rxns], dtype=int).reshape(-1, 2)
    mfe_col_rg = {(r, g): c for r, g, c in zip(mfeCol['rsn'], mfeCol['gene'], mfeCol['Col'])}
    minGiRow = pd.DataFrame({'rsn': rg_r1cpx[:, 0], 'gsn': rg_r1cpx[:, 1],
                             'rxn': [react_id[r] for r in rg_r1cpx[:, 0]],
                             'gene': [allGenes[g] for g in rg_r1cpx[:, 1]],
                             'stoich': rgst.values[rg_r1cpx[:, 0], rg_r1cpx[:, 1]],
                             'exprCol': g1.loc[rg_r1cpx[:, 0], 'exprCol'].values})
    minGiRow['geneCol'] = [mfe_col_rg[(r, g)] if g in mfe_genes else gene_col[g]
                           for r, g in zip(minGiRow['rsn'], minGiRow['gene'])]
    if len(minGiRow) > 0:
        minGiRow['rowind'] = nrow + np.arange(len(minGiRow))
        minGiRow['consNames'] = ['min_%s_%s_%s' % (g, r, _fmt(s)) for g, r, s in
                                 zip(minGiRow['gene'], minGiRow['rxn'], minGiRow['stoich'])]
        setval(minGiRow['rowind'].values, minGiRow['exprCol'].values.astype(int), -1 * minGiRow['stoich'].values)
        setval(minGiRow['rowind'].values, minGiRow['geneCol'].values, 1.0)
        nrow += len(minGiRow)
        rnames += list(minGiRow['consNames'])
        rlb += [0.0] * len(minGiRow)
        rub += [0.0] * len(minGiRow)
        rtype += ['L'] * len(minGiRow)
    else:
        minGiRow['rowind'] = pd.Series(dtype=int)
        minGiRow['consNames'] = pd.Series(dtype=object)
    ## 3-mfe_ constr---------------------------------------------------------------------
    #Add multiple reactions constraints:  sum(gr_gene_r) - gmr_gene = 0
    if verboseMode >= 2:
        print('Adding multiple reactions constraints')
    mfe_g = [allGenes[g] for g in np.nonzero(genecnt > 1)[0]]
    mfeRow = pd.DataFrame({'gene': mfe_g, 'rowind': nrow + np.arange(len(mfe_g)),
                           'Col': [gene_col[g] for g in mfe_g]})
    mfe_row_map = mfeRow.set_index('gene')
    mfeCol['grow.gene'] = mfeCol['gene']
    mfeCol['grow.rowind'] = mfe_row_map.loc[mfeCol['gene'], 'rowind'].values if len(mfeCol) else []
    mfeCol['grow.Col'] = mfe_row_map.loc[mfeCol['gene'], 'Col'].values if len(mfeCol) else []
    if len(mfeRow) > 0:
        setval(mfeCol['grow.rowind'].values, mfeCol['Col'].values, 1.0)  ## geneRxn
        setval(mfeRow['rowind'].values, mfeRow['Col'].values, -1.0)  ## total gene
        nrow += len(mfeRow)
        rnames += ['mfe_' + g for g in mfeRow['gene']]
        rlb += [0.0] * len(mfeRow)
        rub += [0.0] * len(mfeRow)
        rtype += ['E'] * len(mfeRow)

    ##############II. Add MW constraint------------------
    ## Sum(gi. MWi)<= C[gdwpr/gDW]
    fluxind = np.arange(n)
    clb = np.zeros(ncol)
    cub = np.zeros(ncol)
    clb[fluxind] = mod2.lowbnd
    cub[fluxind] = mod2.uppbnd
    cub[n:] = 1000  # C is less than 1, can be 2017
    LB = _gene_bounds(geneConcLB, allGenes, 0)
    UB = _gene_bounds(geneConcUB, allGenes, 1000)
    clb[geneCol['Col'].values] = LB[geneCol['g'].values]
    cub[geneCol['Col'].values] = UB[geneCol['g'].values]

    mw_map = {}
    for g, v in zip(MW.iloc[:, 0], MW.iloc[:, 1]):
        mw_map.setdefault(str(g).strip(), float(v))
    gMW = np.array([mw_map.get(g, np.nan) for g in allGenes])
    nomw = (geneind >= 0) & np.isnan(gMW)
    if nomw.sum() > 0 and verboseMode > 0:
        print('Warning : some genes have no MW:%d (%s), they are not in the capacity constraint'
              % (nomw.sum(), ', '.join(np.array(allGenes)[nomw][:10])))
    cflg = (geneind >= 0) & ~np.isnan(gMW)
    lcind = geneind[cflg]
    lnz = gMW[cflg]
    obj_cf = np.zeros(ncol)
    if RHS is not None:  #test value of upper bound (C)
        if verboseMode >= 2:
            print('Adding capacity constraint.')
        cc_ind = nrow
        setval(cc_ind, lcind, lnz)
        if C_mu_coef != 0:  ## C = RHS + C_mu_coef*objective
            oc = np.nonzero(mod2.obj_coef)[0]
            setval(cc_ind, oc, -C_mu_coef * mod2.obj_coef[oc])
        nrow += 1
        rnames.append('cc')
        rlb.append(0.0)
        rub.append(float(RHS))
        rtype.append('U')
        obj_cf[fluxind] = mod2.obj_coef
        mod_objc_ind = None
        obj_dir = mod2.obj_dir  # maximize the objective of the model (e.g. growth)
    else:  # Add constraint on FBA objective: objective >= objVal, minimize the used budget
        if objVal is None:  #calculate minimum cost and gene concentration
            fba = _build_problem(mod2.S, [0.0] * m, [0.0] * m, ['E'] * m, mod2.lowbnd, mod2.uppbnd,
                                 cnames[:n], rnames[:m], mod2.obj_coef, mod2.obj_dir)
            st = fba.optimize()
            if st != 'optimal':
                raise RuntimeError("Failed to calculate objective value, sol status: %s" % st)
            objVal = fba.objective.value
            if verboseMode > 0:
                print('Objective value to be used: %f' % objVal)
        if objVal is None or np.isnan(objVal):
            raise ValueError('objVal parameter is NA.')
        mod_objc_ind = nrow
        oc = np.nonzero(mod2.obj_coef)[0]
        setval(mod_objc_ind, oc, mod2.obj_coef[oc])
        nrow += 1
        rnames.append('ObjC_')
        # R adds this row with type 'U' (and bounds objVal, Inf) in a wrong row; the intent is objective >= objVal
        rlb.append(float(objVal) if mod2.obj_dir == 'max' else -np.inf)
        rub.append(np.inf if mod2.obj_dir == 'max' else float(objVal))
        rtype.append('L' if mod2.obj_dir == 'max' else 'U')
        cc_ind = None
        obj_dir = 'min'  # minimize the budget
        obj_cf[lcind] = lnz

    # LHS = [S(mod2) 0; ccFBA rows]
    S = mod2.S.tocoo()
    er = np.fromiter((k[0] for k in entries), dtype=int, count=len(entries))
    ec = np.fromiter((k[1] for k in entries), dtype=int, count=len(entries))
    ev = np.fromiter(entries.values(), dtype=float, count=len(entries))
    LHS = sp.csr_matrix((np.concatenate([S.data, ev]), (np.concatenate([S.row, er]), np.concatenate([S.col, ec]))),
                        shape=(nrow, ncol))

    #######################################################
    rg_ind2 = np.argwhere(mod2.rxnGeneMat.toarray().T != 0)[:, ::-1]  # (irrev rxn, gene), in R's order
    r2, g2 = rg_ind2[:, 0], rg_ind2[:, 1]
    rev_r = mod2.irrev2rev[r2]
    g1r = g1.reindex(rev_r)
    rxnGeneCol = pd.DataFrame({'rsn': r2, 'gsn': g2, 'rxn_id': [mod2.react_id[i] for i in r2],
                               'gene': [allGenes[g] for g in g2], 'revRxn': rev_r, 'rNG': rNG[rev_r],
                               'genercnt': genecnt[g2], 'rxnCol': fluxind[r2],
                               'gCol': pd.array([g if g >= 0 else pd.NA for g in geneind[g2]], dtype='Int64'),
                               'gpr.exprCol': g1r['exprCol'].values,
                               'gpr.MW': g1r['MW'].values, 'gpr.mod_cpx_MW1s': g1r['mod_cpx_MW1s'].values,
                               'gpr.genes': g1r['genes'].values,
                               'kcat': kcat_ir[r2], 'kcat_src': kcat_src_ir[r2], 'gMW': gMW[g2],
                               'stoich': rgst.values[rev_r, g2]})
    mfe_col_rg_ix = {(r, g): c for r, g, c in zip(mfeCol['rsn'], mfeCol['g'], mfeCol['Col'])}
    grCol = [mfe_col_rg_ix.get((r, g), pd.NA) if genecnt[g] > 1 else (geneind[g] if geneind[g] >= 0 else pd.NA)
             for r, g in zip(rev_r, g2)]
    rxnGeneCol['grCol'] = pd.array(grCol, dtype='Int64')
    cn = np.array(cnames, dtype=object)

    def colname(c):
        return [None if pd.isna(x) else cn[int(x)] for x in c]
    rxnGeneCol['exprCname'] = colname(rxnGeneCol['gpr.exprCol'])
    rxnGeneCol['gColname'] = colname(rxnGeneCol['gCol'])
    rxnGeneCol['grColname'] = colname(rxnGeneCol['grCol'])

    ###
    ## cnames: rxn_id,g_gene,cmplx_rxn,g_gene_rxn,iso_rxn
    assert LHS.shape == (len(rlb), len(cnames)) == (len(rub), len(clb))
    out = {'LHS': LHS,
           'rlb': np.array(rlb), 'rub': np.array(rub), 'rnames': rnames, 'rtype': rtype,  # rows
           'clb': clb, 'cub': cub, 'cnames': cnames,  # columns
           'fluxind': fluxind,
           'geneind': pd.Series(pd.array([g if g >= 0 else pd.NA for g in geneind], dtype='Int64'), index=allGenes),
           'rxnkcatRow': rxnkcatRow, 'minGiRow': minGiRow, 'mfeRow': mfeRow,
           'rgst': rgst, 'mod_cpx_mw': mod_cpx_mw,  # complex stoichiometry
           'rxnExprCol': rxnExprCol,  # flxCol,exprCol,rowind,reaction total expression
           'rxnGeneCol': rxnGeneCol,  # rxn_id,gene,Col : used to get MCrowding
           'geneCol': geneCol, 'mfeCol': mfeCol, 'gpr1iso': gpr1iso,  # low level structure, used for debugging
           'kcat_ir': kcat_ir, 'kcat_src_ir': kcat_src_ir, 'geneMW': gMW,
           'model': mod2,
           'obj_cf': obj_cf, 'obj_dir': obj_dir,
           'mod_objc_ind': mod_objc_ind,  # row index of constraint on FBA model objective
           'cc_ind': cc_ind,  # capacity constraint index
           }
    if verboseMode > 0:
        print("LP: %d rows, %d columns (%d genes, %d multi-functional parts, %d complexes)"
              % (nrow, ncol, len(geneCol), len(mfeCol), len(gprCol)))
    return out


def solve_ccFBA_mat(mat, mod2=None, interface=None, verboseMode=1, prob=None):
    """Solve the problem returned by getccFBA_mat.

    The bounds of the flux columns are taken from the irreversible model (mod2, default mat['model']),
    so they can be changed there before solving (e.g. uppbnd of an uptake reaction), as in R.
    interface: optlang interface (default: cobra's default solver).
    prob: problem of an earlier call with the same mat (Res['prob']) to be reused; only the bounds of
          the flux columns are updated (faster when only the bounds of mod2 change).
    Returns a dictionary: obj, status, x (values of all columns), rxnFlx (net flux of each reaction of
    the original model), flx_ir (fluxes of mod2), geneConc (gene -> concentration of the gene
    variables), prob (the optlang problem)."""
    mod2 = mat['model'] if mod2 is None else mod2
    n = mod2.react_num
    clb, cub = mat['clb'].copy(), mat['cub'].copy()
    clb[:n] = mod2.lowbnd
    cub[:n] = mod2.uppbnd
    if verboseMode > 0:
        print("Solving LP ...")
    if prob is None:
        prob = _build_problem(mat['LHS'], mat['rlb'], mat['rub'], mat['rtype'], clb, cub, mat['cnames'],
                              mat['rnames'], mat['obj_cf'], mat['obj_dir'], interface)
    else:
        if len(prob.variables) != len(mat['cnames']):
            raise ValueError("prob does not correspond to mat")
        for j in range(n):
            v = prob.variables[j]
            if v.lb != clb[j] or v.ub != cub[j]:
                v.set_bounds(clb[j], cub[j])
        prob.update()
    status = prob.optimize()
    x = np.array([v.primal for v in prob.variables], dtype=float) if status == 'optimal' else \
        np.full(len(mat['cnames']), np.nan)
    obj = prob.objective.value if status == 'optimal' else np.nan
    rxnFlx = dict(zip(mod2.rxn_id_orig, mod2.net_flux(x)))
    geneConc = dict(zip(mat['geneCol']['gene'], x[mat['geneCol']['Col'].values]))
    return {'obj': obj, 'status': status, 'x': x, 'rxnFlx': rxnFlx, 'flx_ir': x[:n],
            'geneConc': geneConc, 'prob': prob}
