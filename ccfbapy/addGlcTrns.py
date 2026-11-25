# Add a glucose transport constraint (hexose transporters of S. cerevisiae) to the problem.
# by Abdelmoneim Amer Desouki
# part of ccFBApy, python version of addGlcTrns (sybilccFBA, R/addGlcTrns.R)
#
# The R function is not complete: the columns of the transporters are not added (that code is commented out),
# so trnsCol is NULL and the new row cannot be built. This version adds the columns as intended in the R code.

# name: (gene, kcat, Km); capacity = kcat*Si/(Si+Km). Values as in R (Gal2 only with galactose).
HXT = {'Hxt1': ('YHR094C', 41.0, 107.0),
       'Hxt2': ('YMR011W', 16.1, 2.9),
       'Hxt3': ('YDR345C', 18.5, 29.0),
       'Hxt4': ('YHR092C', 12.0, 6.2),
       'Hxt5': ('YHR096C', 14.0, 10.0),
       'Hxt6': ('YDR343C', 11.4, 1.5),
       'Hxt7': ('YDR342C', 11.7, 1.3),
       'Gal2': ('YLR081W', 17.5, 1.5)}
DEFAULT_TRANSPORTERS = ['Hxt1', 'Hxt2', 'Hxt3', 'Hxt4', 'Hxt5', 'Hxt6', 'Hxt7']


def glc_trns_capacity(Si=0.2, transporters=None):
    """Capacity of each transporter for the glucose concentration Si: kcat*Si/(Si+Km). Returns a dict name -> value."""
    transporters = DEFAULT_TRANSPORTERS if transporters is None else transporters
    return {t: HXT[t][1] * Si / (Si + HXT[t][2]) for t in transporters}


def _flux_var(prob, glc_rxn):
    """Variable of the (forward) flux of glc_rxn in prob: the reaction id (cobra model: cfba_moment, ccFBA_mr),
    'R_' + id (matrix of getccFBA_mat) or the forward direction '<id>_f' of a split reaction."""
    for name in (glc_rxn, 'R_' + glc_rxn, glc_rxn + '_f', 'R_' + glc_rxn + '_f'):
        if name in prob.variables:
            return prob.variables[name]
    raise ValueError("glucose transport reaction '%s' not found in the problem" % glc_rxn)


def addGlcTrns(prob, mod2=None, glc_rxn='R_GLCt1', Si=0.2, transporters=None, trns_ub=1000, trnsCol=None):
    """Add a glucose transport constraint to the problem: an upper bound on glucose consumption.

    prob:     optlang problem, e.g. res['prob'] of cfba_moment / ccFBA_mr, or solve_ccFBA_mat(...)['prob'].
              It is changed in place.
    mod2:     accepted for compatibility with R, not needed (the flux variable is found by name).
    glc_rxn:  id of the glucose transport reaction (R: "R_GLCt1" of iMM904). Its forward flux is constrained.
    Si:       glucose concentration (as the Km values), default 0.2 as in R.
    transporters: names of the transporters (keys of HXT), default Hxt1..Hxt7 (as in R; add 'Gal2' when
              galactose is present).
    trns_ub:  upper bound of the new transporter variables (default 1000, as in the R code).
    trnsCol:  optional dict transporter name -> name of an existing variable of prob to use instead of a new
              variable (e.g. {'Hxt1': 'g_YHR094C'} to use the gene concentration of cfba_moment).

    Adds, for each transporter t without an entry in trnsCol, the variable 'trns_<t>' (0 <= x_t <= trns_ub), and
    the constraint 'glcTrns':  v(glc_rxn) - sum_t capacity_t * x_t <= 0, with capacity_t = kcat_t*Si/(Si+Km_t).
    Returns prob (as R).
    """
    transporters = DEFAULT_TRANSPORTERS if transporters is None else list(transporters)
    bad = [t for t in transporters if t not in HXT]
    if bad:
        raise ValueError('unknown transporters: %s (known: %s)' % (', '.join(bad), ', '.join(HXT)))
    if 'glcTrns' in prob.constraints:
        raise ValueError("the problem already has a constraint 'glcTrns'")
    interface = prob.interface
    cap = glc_trns_capacity(Si, transporters)
    trnsCol = {} if trnsCol is None else dict(trnsCol)
    new_vars = []
    cols = {}
    for t in transporters:
        if t in trnsCol:
            if trnsCol[t] not in prob.variables:
                raise ValueError("variable '%s' (trnsCol of %s) not in the problem" % (trnsCol[t], t))
            cols[t] = prob.variables[trnsCol[t]]
        else:
            v = interface.Variable('trns_' + t, lb=0, ub=trns_ub)
            new_vars.append(v)
            cols[t] = v
    if new_vars:
        prob.add(new_vars)
    glcRxn = _flux_var(prob, glc_rxn)
    c = interface.Constraint(0, ub=0, name='glcTrns')
    prob.add(c)
    prob.update()
    coef = {cols[t]: -cap[t] for t in transporters}
    coef[glcRxn] = 1
    c.set_linear_coefficients(coef)
    prob.update()
    # Add to crowding constraint (same budget), Molecular weights required? (open question in the R code)
    return prob
