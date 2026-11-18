#get molecular weight of enzyme complexes and the isoenzyme with minimum MW.
# By Dr. Abdelmoneim Amer Desouki
# 31/3/2024
# Python version of getGpr1iso (sybilccFBA, R/getGpr1iso.R)
import ast
import warnings

import numpy as np
import pandas as pd


def _gpr_terms(gpr):
    """Terms of a cobra GPR in disjunctive normal form: list of sorted lists of unique gene ids.
    'A and B or C' -> [['A','B'], ['C']]. Duplicate terms are removed (first occurrence is kept)."""
    if gpr is None or str(gpr).strip() == "" or gpr.body is None:
        return []

    def dnf(node):
        if isinstance(node, ast.Expression):
            return dnf(node.body)
        if isinstance(node, ast.Name):
            return [[node.id]]
        if isinstance(node, ast.BoolOp) and isinstance(node.op, ast.Or):
            return [t for v in node.values for t in dnf(v)]
        if isinstance(node, ast.BoolOp) and isinstance(node.op, ast.And):
            terms = [[]]
            for v in node.values:
                terms = [t + u for t in terms for u in dnf(v)]
            return terms
        raise ValueError("Unsupported element in GPR rule: %s" % ast.dump(node))

    out = []
    for t in dnf(gpr.body):
        t = sorted(set(t))
        if t not in out:
            out.append(t)
    return out


def _fmt(x):
    """number -> text as R's paste() would show it (NA for missing)."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "NA"
    return format(float(x), ".15g")


def _parse_stoich(st, ngenes, genes):
    """Stoichiometry text '2,1,1' -> list of floats. Returns None when not available."""
    if st is None or (isinstance(st, float) and np.isnan(st)):
        return None
    parts = [p.strip() for p in str(st).split(",")]
    if any(p in ("", "NA", "nan", "NaN") for p in parts):
        warnings.warn("Stoichiometry '%s' of complex %s contains NA, all subunits are set to 1" % (st, genes))
        return None
    if len(parts) != ngenes:
        raise ValueError("cpx_stoich: stoichiometry '%s' has %d values but complex '%s' has %d genes"
                         % (st, len(parts), genes, ngenes))
    return [float(p) for p in parts]


def getGpr1iso(model, MW, cpx_stoich=None):
    """Preprocessing of the model: molecular weight of every enzyme (complex) of each GPR, and
    choice of the isoenzyme with minimum molecular weight, so that the OR's are removed from the GPR's.

    model:      cobra Model
    MW:         DataFrame, first column gene id, second column molecular weight (g/mmol)
    cpx_stoich: optional DataFrame with (at least) the columns 'genes' and 'stoich':
                'genes':  comma separated gene ids of the complex, sorted (e.g. 'b2221,b2222')
                'stoich': comma separated number of subunits of each gene, same order (e.g. '2,2')
                Complexes that are not found get stoichiometry one for all subunits.

    Returns a dictionary:
        gpr1iso:    one row per reaction with a GPR: the chosen term (complex or isoenzyme)
        mod_cpx_mw: all terms of all GPR's with their molecular weight
        rgst:       DataFrame reactions x genes, stoichiometry of the genes of the chosen term
    Indices: 'rsn' is the 0-based index of the reaction in model.reactions.
    gprtype: r1g (one gene), r1cpx (one complex), rmg (isoenzymes: OR of genes),
             rmcpx (OR of complexes / mixed).
    Genes without a molecular weight get MW NaN (NA in R); like in R, a term with unknown
    MW is preferred (e.g. the spontaneous 's0001' in 'b0875 or s0001').
    """
    react_id = [r.id for r in model.reactions]
    allGenes = [g.id for g in model.genes]
    gene_ix = {g: i for i, g in enumerate(allGenes)}

    mw_map = {}
    for g, v in zip(MW.iloc[:, 0], MW.iloc[:, 1]):
        g = str(g).strip()
        if g not in mw_map:  # as R's MW[genes,2]: first row of the gene is used
            mw_map[g] = float(v)

    st_map = {}
    if cpx_stoich is not None:
        for g, s in zip(cpx_stoich['genes'], cpx_stoich['stoich']):
            g = str(g).strip()
            if g not in st_map:  # as R's match(): first row is used
                st_map[g] = s

    ### 0. Preprocessing evaluate complex MW's ---------------------------------------
    rows = []
    for rsn, r in enumerate(model.reactions):
        terms = _gpr_terms(r.gpr)
        if len(terms) == 0:
            continue
        nOR = len(terms) - 1
        nAND = sum(len(t) - 1 for t in terms)
        rl = str(r.gpr)
        if nOR == 0:
            gprtype = 'r1g' if len(terms[0]) == 1 else 'r1cpx'
        else:
            gprtype = 'rmg' if nAND == 0 else 'rmcpx'
        if gprtype == 'rmg':
            # R numbers the isoenzymes in the order of the genes in the model (column order of
            # rxnGeneMat); termNo is used to choose between isoenzymes with equal MW.
            terms = sorted(terms, key=lambda t: gene_ix.get(t[0], len(gene_ix)))
        for j, t in enumerate(terms):
            rows.append({'rsn': rsn, 'rxn': r.id,
                         'bool_expr': ' and '.join(t) if gprtype == 'rmcpx' else rl,
                         'genes': ','.join(t), 'termNo': j + 1, 'nOR': nOR, 'nAND': nAND,
                         'gprtype': gprtype})
    cols = ['rsn', 'rxn', 'bool_expr', 'genes', 'termNo', 'nOR', 'nAND', 'gprtype']
    mod_cpx = pd.DataFrame(rows, columns=cols)
    order = {'r1g': 0, 'r1cpx': 1, 'rmg': 2, 'rmcpx': 3}  # same blocks as in R
    mod_cpx = mod_cpx.assign(_o=mod_cpx['gprtype'].map(order)).sort_values(
        ['_o', 'rsn', 'termNo'], kind='stable').drop(columns='_o').reset_index(drop=True)

    mod_cpx['stoich'] = [st_map.get(g, np.nan) for g in mod_cpx['genes']]
    mod_cpx['stoich'] = [np.nan if (isinstance(s, float) and np.isnan(s)) else str(s).replace(' ', '')
                         for s in mod_cpx['stoich']]

    termNG = []
    MW1s = []    # all ones
    MWtxt = []
    stMW = []    # using stoich.
    missing = set()
    for genes, st in zip(mod_cpx['genes'], mod_cpx['stoich']):
        gl = genes.split(',')
        mws = [mw_map.get(g, np.nan) for g in gl]
        missing.update(g for g in gl if g not in mw_map)
        termNG.append(len(gl))
        MWtxt.append(','.join(_fmt(x) for x in mws))
        MW1s.append(float(np.sum(mws)))  # NaN when a gene has no MW (as sum() in R)
        sv = _parse_stoich(st, len(gl), genes)
        stMW.append(np.nan if sv is None else float(np.sum(np.array(mws) * np.array(sv))))
    if missing:
        warnings.warn("%d gene(s) have no molecular weight: %s" % (len(missing), ', '.join(sorted(missing))))

    mod_cpx_mw = mod_cpx.copy()
    mod_cpx_mw['termNG'] = termNG
    mod_cpx_mw['mod_cpx_MW1s'] = MW1s
    mod_cpx_mw['mod_st_MW'] = stMW
    mod_cpx_mw['mod_cpx_MWtxt'] = MWtxt
    mod_cpx_mw['MW'] = np.where(np.isnan(stMW), MW1s, stMW)

    ### Preprocessing 2. choose best iso enzyme for each rxn  --------------------------------
    # each iso gpr: choose the term with min(MW). As in R (min() of a group with NA is NA), when a
    # term has unknown MW, the term(s) with unknown MW are chosen.
    keep = np.zeros(len(mod_cpx_mw), dtype=bool)
    for rsn, grp in mod_cpx_mw.groupby('rsn', sort=False):
        ix = grp.index.values
        if grp['nOR'].iloc[0] == 0:
            keep[ix] = True
            continue
        mwv = grp['MW'].values
        if np.isnan(mwv).any():
            keep[ix[np.isnan(mwv)]] = True
        else:
            keep[ix[np.isclose(mwv, mwv.min(), rtol=1e-12, atol=0)]] = True
    gpr1iso = mod_cpx_mw[keep].copy()
    # equal iso enzymes: take the first term
    if gpr1iso['rsn'].duplicated().any():
        first = gpr1iso.groupby('rsn')['termNo'].transform('min')
        gpr1iso = gpr1iso[(gpr1iso['nOR'] == 0) | (gpr1iso['termNo'] == first)]

    st = gpr1iso['stoich'].astype(object)
    st[st.isna() & (gpr1iso['termNG'] == 1)] = '1'
    gpr1iso['stoich'] = st

    rgst = pd.DataFrame(np.zeros((len(react_id), len(allGenes))), index=react_id, columns=allGenes)
    rgst_v = rgst.values
    for i in gpr1iso.index:
        genes = gpr1iso.at[i, 'genes'].split(',')
        sv = _parse_stoich(gpr1iso.at[i, 'stoich'], len(genes), gpr1iso.at[i, 'genes'])
        if sv is None:
            sv = [1.0] * len(genes)
            gpr1iso.at[i, 'stoich'] = ','.join(['1'] * len(genes))
        rsn = gpr1iso.at[i, 'rsn']
        for g, s in zip(genes, sv):
            rgst_v[rsn, gene_ix[g]] = s
    rgst = pd.DataFrame(rgst_v, index=react_id, columns=allGenes)

    return {'gpr1iso': gpr1iso, 'mod_cpx_mw': mod_cpx_mw, 'rgst': rgst}
