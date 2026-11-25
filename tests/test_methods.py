import pandas as pd
import pytest

import ccfbapy
from ccfbapy import (cfba_moment, ccFBA_mr, cfba_moment_pw, getccFBA_mat, solve_ccFBA_mat, getccFBA_sol,
                     ccFBA_iso, getGpr1iso, mod2irrev, simulate_EColi, simGeneOverExpr, addGlcTrns,
                     glc_trns_capacity)


def test_version():
    assert ccfbapy.__version__
    assert ccfbapy.cfba_moment_mr is ccfbapy.ccFBA_mr


def test_moment_overestimates_shared_subunit(toy, toy_kw):
    # A is used by both complexes: 2 complexes for the price of 3 genes
    assert cfba_moment(toy, **toy_kw)['obj'] == pytest.approx(240)
    assert ccFBA_mr(toy, **toy_kw)['obj'] == pytest.approx(240)
    assert cfba_moment_pw(toy, **toy_kw)['obj'] == pytest.approx(240)


def test_moment_does_not_change_model(toy, toy_kw):
    n = len(toy.solver.variables)
    cfba_moment(toy, **toy_kw)
    assert len(toy.solver.variables) == n


def test_gpr1iso_and_iso(toy, toy_kw):
    assert solve_ccFBA_mat(getccFBA_mat(toy, **toy_kw), verboseMode=0)['obj'] == pytest.approx(180)
    res = ccFBA_iso(toy, **toy_kw)
    assert res['status'] == 'optimal'
    assert res['obj'] == pytest.approx(180)
    assert res['rxnMC'].sum() == pytest.approx(0.1)
    assert getccFBA_sol(toy, **toy_kw)['obj'] == pytest.approx(180)


def test_knock_out_isoenzyme(toy, toy_kw):
    assert solve_ccFBA_mat(getccFBA_mat(toy, geneConcUB={'B': 0}, **toy_kw), verboseMode=0)['obj'] == \
        pytest.approx(0, abs=1e-9)
    res = ccFBA_iso(toy, geneConcUB={'B': 0}, **toy_kw)
    assert res['obj'] == pytest.approx(180)
    assert res['geneConc']['B'] == pytest.approx(0, abs=1e-12)
    # C heavier and B limited: A-B up to the limit, then A-C (3600 * (0.02 + 0.06 / 2.5) = 158.4)
    kw = dict(toy_kw, MW=pd.DataFrame({'gene': ['A', 'B', 'C'], 'mw': [1.0, 1.0, 1.5]}))
    res = ccFBA_iso(toy, geneConcUB={'B': 0.02}, **kw)
    assert res['obj'] == pytest.approx(158.4)
    assert list(res['enzConc']['conc'].round(6)) == [0.02, 0.024]
    assert solve_ccFBA_mat(getccFBA_mat(toy, geneConcUB={'B': 0.02}, **kw), verboseMode=0)['obj'] == \
        pytest.approx(72)


def test_iso_reuse_matrix(toy, toy_kw):
    mat = ccfbapy.getccFBA_iso_mat(toy, **toy_kw)
    assert ccFBA_iso(toy, mat=mat, geneConcUB={'A': 0.01}, verboseMode=0)['obj'] == pytest.approx(36)


def test_min_budget(toy, toy_kw):
    kw = dict(toy_kw, RHS=None)
    res = ccFBA_iso(toy, objVal=180, **kw)
    assert res['obj'] == pytest.approx(0.1)
    assert solve_ccFBA_mat(getccFBA_mat(toy, objVal=180, **kw), verboseMode=0)['obj'] == pytest.approx(0.1)


def test_getGpr1iso(toy, toy_kw):
    r = getGpr1iso(toy, MW=toy_kw['MW'])
    g = r['gpr1iso']
    assert list(g['rxn']) == ['R1']
    assert g['gprtype'].iloc[0] == 'rmcpx'
    assert g['MW'].iloc[0] == pytest.approx(2)
    assert r['rgst'].loc['R1'].sum() == 2


def test_mod2irrev(toy):
    mod2 = mod2irrev(toy)
    assert mod2.react_id == ['EX_s_f', 'EX_s_b', 'R1', 'EX_p']
    assert list(mod2.net_flux([1, 3, 2, 2])) == [-2, 2, 2]


def test_simulate_EColi(toy, toy_kw):
    mod2 = mod2irrev(toy)
    res = simulate_EColi(toy, mod2, kcat=toy_kw['Kcat'], mw=toy_kw['MW'], budget_C=0.1, CSList=['EX_s_b'],
                         atpmz=False, verboseMode=0)
    assert res['gr_data']['flx'].iloc[0] == pytest.approx(180)
    assert res['uptake']['flx'].iloc[0] == pytest.approx(180)


def test_simGeneOverExpr(toy, toy_kw):
    kw = dict(toy_kw)
    RR = simGeneOverExpr('C', 0.5, toy, RHS=kw.pop('RHS'), medval=3600, **kw)
    assert RR['OverExprEff']
    assert RR['gxMC_OEx'] == pytest.approx(0.05)


def test_addGlcTrns(toy):
    cap = glc_trns_capacity()
    assert sum(cap.values()) == pytest.approx(4.7926, abs=1e-4)
    m = toy.copy()
    addGlcTrns(m.solver, glc_rxn='R1', trns_ub=1)
    assert m.slim_optimize() == pytest.approx(sum(cap.values()))
    assert 'glcTrns' in m.solver.constraints
