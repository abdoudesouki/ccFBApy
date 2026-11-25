# Read the MOMENT model exported from MATLAB and solve it.
# by Abdelmoneim Amer Desouki
# part of ccFBApy, Python version of readmodel (sybilccFBA, R/readmodel.R)
import os

import numpy as np
import pandas as pd
import scipy.sparse as sp

# example files of the R package (R: system.file("extdata", package="sybilccFBA")), shipped with ccfbapy
EXTDATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'extdata')


def read_moment_files(path=EXTDATA):
    """Read the files of the MOMENT model exported from MATLAB (separator ';'):
    mat.csv (row;column;value, 1-based, no header), mets.csv, rxns.csv, rbnds.csv, cbnds.csv.
    Default folder: the extdata folder of the package (ccfbapy.datasets.extdata_path()).
    Returns a dict with mat, mets, rxns, rbnds, cbnds (DataFrames),
    which can be passed to readmodel(**files).

    mat.csv has no header line. The R source reads it with read.table(..., header=T), which takes the first
    entry of the matrix as the header (and so drops it); here it is read as data."""
    def rd(name, **kw):
        # no quoting, as quote="" in R (names contain quotes, e.g. "2',3'-cyclic")
        return pd.read_csv(os.path.join(path, name), sep=';', quoting=3, **kw)
    mat = rd('mat.csv', header=None, names=['i', 'j', 'val'])
    # 6 lines of mets.csv (the alcohols added at the end, e.g. 'isobutanol') have one field more: it is dropped
    mets = rd('mets.csv', keep_default_na=False, engine='python', on_bad_lines=lambda fields: fields[:4])
    return {'mat': mat, 'mets': mets, 'rxns': rd('rxns.csv', keep_default_na=False),
            'rbnds': rd('rbnds.csv'), 'cbnds': rd('cbnds.csv')}


def readmodel(mat, mets, rxns, rbnds, cbnds, solver="glpk", C=0.27, lp_fname=None, verboseMode=1):
    """Create the LP of the MOMENT model exported from MATLAB, solve it and return the fluxes.

    mat:    constraint matrix as triplets: row, column, value (1-based indices, as in MATLAB/R)
    mets:   metabolites (one per equality row: the first len(mets) rows are S v = b)
    rxns:   reactions (one per flux column: the first len(rxns) columns); column 'ocf' is the objective
    rbnds:  bounds of the rows (lb, ub): metabolite rows equal lb, the other rows are <= ub
    cbnds:  bounds of the columns (lb, ub)
    solver: optlang interface: 'glpk' (R: "glpkAPI"), 'cplex' (R: "cplexAPI"), 'gurobi', ... or an interface module
    C:      budget, the upper bound of the last row (the enzyme capacity constraint), 0.27 as in R (48% of the
            protein is metabolic). None: use the value of rbnds (0.56 in the example files).
    lp_fname: when given, the problem is written to this file in LP format (R always writes an LP file).
    The objective is maximized.

    The sizes are taken from the inputs (R has them fixed in the code: 6705 rows, 1674 metabolites,
    4991 columns, 3234 fluxes, which are the sizes of the example files).

    Returns a dict: x (values of all columns, as the returned vector of R), obj, status, rxnFlx (Series: values of
    the flux columns with the reaction ids of rxns), prob (the optlang problem). Columns are named with the
    reaction ids (fluxes) and 'x<j>' (other columns, j 1-based); rows with the metabolite ids and 'r<i>'; the last
    row is named 'MW'.
    """
    import optlang
    interface = solver
    if isinstance(solver, str):
        name = {'glpkapi': 'glpk', 'cplexapi': 'cplex'}.get(solver.lower(), solver.lower())
        interface = getattr(optlang, name + '_interface', None)
        if interface is None:
            import importlib
            interface = importlib.import_module('optlang.%s_interface' % name)

    nr, nc, nf, nm = len(rbnds), len(cbnds), len(rxns), len(mets)
    mat = pd.DataFrame(np.asarray(mat)[:, :3], columns=['i', 'j', 'val'])
    ii = mat['i'].astype(int).values - 1
    jj = mat['j'].astype(int).values - 1
    if ii.min() < 0 or ii.max() >= nr or jj.min() < 0 or jj.max() >= nc:
        raise ValueError('mat has indices out of the size of rbnds (%d rows) / cbnds (%d columns)' % (nr, nc))
    # as LHS[i,j] = val in R: the last value of a repeated (i, j) is kept
    keep = ~pd.DataFrame({'i': ii, 'j': jj}).duplicated(keep='last').values
    LHS = sp.csr_matrix((mat['val'].astype(float).values[keep], (ii[keep], jj[keep])), shape=(nr, nc))

    cobj = np.concatenate([pd.to_numeric(rxns['ocf']).values.astype(float), np.zeros(nc - nf)])
    rb = rbnds.iloc[:, :2].astype(float).values
    cb = cbnds.iloc[:, :2].astype(float).values
    rupper = rb[:, 1].copy()
    if C is not None:
        rupper[nr - 1] = C  # set to limit in paper 48% of protein is metabolic

    rids = [str(x).strip() for x in rxns.iloc[:, 0]]
    cnames = rids + ['x%d' % (j + 1) for j in range(nf, nc)]
    rnames = [str(x).strip() for x in mets.iloc[:, 0]] + ['r%d' % (i + 1) for i in range(nm, nr - 1)] + ['MW']
    if len(set(cnames)) < nc or len(set(rnames)) < nr:  # names must be unique in optlang
        cnames = ['x%d' % (j + 1) for j in range(nc)]
        rnames = ['r%d' % (i + 1) for i in range(nr)]

    if verboseMode > 0:
        print('readmodel: %d rows (%d metabolites), %d columns (%d fluxes), %d non zeros'
              % (nr, nm, nc, nf, LHS.nnz))
    prob = interface.Model(name='Moment')
    vars_ = []
    for j in range(nc):
        lb, ub = cb[j]
        vars_.append(interface.Variable(cnames[j], lb=None if np.isinf(lb) else lb, ub=None if np.isinf(ub) else ub))
    prob.add(vars_)
    # rows 1..nm: S v = b (fixed to lb), the other rows: <= ub
    cons = [interface.Constraint(0, lb=rb[i, 0], ub=rb[i, 0], name=rnames[i]) if i < nm else
            interface.Constraint(0, ub=rupper[i], name=rnames[i]) for i in range(nr)]
    prob.add(cons)
    prob.update()
    for i, c in enumerate(cons):
        s, e = LHS.indptr[i], LHS.indptr[i + 1]
        if e > s:
            c.set_linear_coefficients({vars_[j]: v for j, v in zip(LHS.indices[s:e], LHS.data[s:e])})
    prob.objective = interface.Objective(0, direction='max')
    prob.objective.set_linear_coefficients({vars_[j]: cobj[j] for j in np.nonzero(cobj)[0]})
    prob.update()

    if lp_fname is not None:
        if verboseMode > 0:
            print('Writing problem to file: %s ...' % lp_fname)
        with open(lp_fname, 'w') as f:
            f.write(prob.to_lp())

    status = prob.optimize()
    if status == 'optimal':
        x = np.array([v.primal for v in prob.variables], dtype=float)
        obj = prob.objective.value
    else:
        x = np.full(nc, np.nan)
        obj = np.nan
    rxnFlx = pd.Series(x[:nf], index=rids)
    if verboseMode > 0:
        msg = ['%s=%f' % (lbl, rxnFlx[r]) for lbl, r in
               (('GLC upt', 'REV_EX_glc(e)'), ('AC', 'REV_EX_ac(e)'), ('Pyr', 'REV_EX_pyr(e)'),
                ('fruc', 'REV_EX_fru(e)'), ('Lac', 'REV_EX_lac_L(e)'), ('galt', 'REV_EX_galt(e)')) if r in rxnFlx.index]
        print('readmodel: %s, objective=%s%s' % (status, obj, (', ' + ', '.join(msg)) if msg else ''))
    return {'x': x, 'obj': obj, 'status': status, 'rxnFlx': rxnFlx, 'prob': prob}
