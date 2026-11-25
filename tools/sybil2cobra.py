# Read the data sets of the R package (R/sybilccFBA/data/*.RData) in Python: development tool of ccFBApy,
# used by tools/make_package_data.py to export them to SBML and CSV (the package ccfbapy does not need rdata).
# sybil models (class modelorg) are converted to cobra models, data frames to pandas DataFrames.
# Needs the package rdata (pip install rdata).
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp


def _str(x):
    """first element of an R vector as a python string ('' when missing)."""
    if x is None:
        return ''
    x = np.asarray(x).ravel()
    return '' if len(x) == 0 or (isinstance(x[0], float) and np.isnan(x[0])) else str(x[0])


def _dgc(m):
    """R Matrix dgCMatrix/lgCMatrix (as read by rdata) -> scipy csc matrix."""
    nr, nc = (int(d) for d in m.Dim)
    x = getattr(m, 'x', None)
    x = np.ones(len(m.i)) if x is None else np.asarray(x, dtype=float)
    return sp.csc_matrix((x, np.asarray(m.i), np.asarray(m.p)), shape=(nr, nc))


def sybil2cobra(mo):
    """sybil modelorg (as read by rdata) -> cobra Model.
    Ids are kept as in sybil (e.g. reaction 'EX_glc__D(e)', metabolite 'glc__D[e]'), so that the reaction
    ids of the R examples can be used. Bounds, objective, GPR rules, names, compartments, formulas and
    charges are copied. react_rev of sybil is not kept: in cobra a reaction is reversible when its
    lower bound is < 0 (in the models of the package react_rev is the same as lowbnd < 0)."""
    import cobra
    model = cobra.Model(_str(mo.mod_id))
    model.name = _str(mo.mod_name)
    comps = [str(c) for c in np.asarray(mo.mod_compart).ravel()]
    model.compartments = {c: c for c in comps}
    met_attr = getattr(mo, 'met_attr', None)
    met_attr = met_attr if isinstance(met_attr, pd.DataFrame) and len(met_attr) == int(mo.met_num[0]) else None
    mets = []
    for i, mid in enumerate(np.asarray(mo.met_id)):
        ci = int(np.asarray(mo.met_comp)[i]) - 1
        m = cobra.Metabolite(str(mid), name=str(np.asarray(mo.met_name)[i]),
                             compartment=comps[ci] if 0 <= ci < len(comps) else None)
        if met_attr is not None:
            if 'chemicalFormula' in met_attr and isinstance(met_attr['chemicalFormula'].iloc[i], str):
                m.formula = met_attr['chemicalFormula'].iloc[i]
            if 'charge' in met_attr and pd.notna(met_attr['charge'].iloc[i]):
                m.charge = int(met_attr['charge'].iloc[i])
        mets.append(m)
    model.add_metabolites(mets)

    S = _dgc(mo.S)
    lb, ub = np.asarray(mo.lowbnd, dtype=float), np.asarray(mo.uppbnd, dtype=float)
    rids, rnames = np.asarray(mo.react_id), np.asarray(mo.react_name)
    gpr = np.asarray(mo.gpr)
    rxns = []
    for j, rid in enumerate(rids):
        r = cobra.Reaction(str(rid), name=str(rnames[j]), lower_bound=float(lb[j]), upper_bound=float(ub[j]))
        s, e = S.indptr[j], S.indptr[j + 1]
        r.add_metabolites({mets[i]: float(v) for i, v in zip(S.indices[s:e], S.data[s:e])})
        rxns.append(r)
    model.add_reactions(rxns)
    # genes in the order of allGenes (as in R); otherwise cobra adds them in the order of a set, which changes
    # between python sessions
    model.genes.extend([cobra.Gene(str(g)) for g in np.asarray(mo.allGenes).ravel()])
    for j, r in enumerate(model.reactions):
        g = str(gpr[j]) if isinstance(gpr[j], str) else ''
        if g.strip():
            r.gene_reaction_rule = g
    oc = np.asarray(mo.obj_coef, dtype=float)
    model.objective = {model.reactions[j]: float(oc[j]) for j in np.nonzero(oc)[0]}
    return model


def read_RData(fname):
    """Objects of an .RData file: sybil models -> cobra Model, data frames -> DataFrame, other objects as
    returned by rdata. Returns a dictionary name -> object."""
    try:
        import rdata
    except ImportError as e:
        raise ImportError("reading .RData files needs the package rdata: pip install rdata"
                          ) from e
    with warnings.catch_warnings():  # rdata warns for the S4 classes it has no constructor for
        warnings.filterwarnings('ignore', message='Missing constructor for R class')
        objs = rdata.conversion.convert(rdata.parser.parse_file(fname))
    out = {}
    for k, v in objs.items():
        cls = getattr(v, 'class', None)
        if cls is not None and 'modelorg' in [str(c) for c in np.asarray(cls).ravel()]:
            v = sybil2cobra(v)
        elif isinstance(v, pd.DataFrame):  # R factors -> strings (as stringsAsFactors=FALSE)
            v = v.copy()
            for c in v.columns:
                if isinstance(v[c].dtype, pd.CategoricalDtype):
                    v[c] = v[c].astype(object).where(v[c].notna(), None)
            v.columns = [str(c) for c in v.columns]
        out[k] = v
    return out
