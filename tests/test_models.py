"""Tests with the models and tables of sybilccFBA shipped with the package (SBML and CSV files)."""
import sys

import pytest

import ccfbapy
from ccfbapy import datasets


@pytest.fixture(scope='module')
def iaf():
    m = datasets.load_model('iAF1260')
    m.reactions.get_by_id('R_EX_glc_e_').lower_bound = -1000
    m.reactions.get_by_id('R_EX_o2_e_').lower_bound = -1000
    m.reactions.get_by_id('R_ATPM').lower_bound = 0
    return m


def test_no_rdata():
    assert 'rdata' not in sys.modules
    assert not hasattr(ccfbapy, 'read_RData')


def test_datasets():
    assert set(datasets.list_datasets()) == set(datasets.MODELS) | set(datasets.TABLES)
    assert len(datasets.load_dataset('kcat')) == 508
    assert len(datasets.load_dataset('mw')) == 1260
    assert len(datasets.load_dataset('mw_iML1515')) == 1515
    assert list(datasets.load_dataset('yst_kcat').columns) == ['rxn_id', 'val', 'dirxn']
    assert len(datasets.load_dataset('yst_mw')) == 5899
    with pytest.raises(ValueError):
        datasets.load_model('kcat')
    with pytest.raises(ValueError):
        datasets.model_path('iAF1260', source='R')


def test_models():
    m = datasets.load_model('iMM904')
    assert (len(m.reactions), len(m.metabolites), len(m.genes)) == (1577, 1228, 904)
    assert m.slim_optimize() == pytest.approx(0.287866, abs=1e-6)
    b = datasets.load_model('iMM904', source='bigg')
    assert 'EX_glc__D_e' in b.reactions


def test_kcat_list():
    m = datasets.load_model('iML1515')
    assert len(datasets.load_kcat_list(m)) == 651


def test_iAF1260_methods(iaf):
    kw = dict(Kcat=datasets.load_dataset('kcat'), MW=datasets.load_dataset('mw'), RHS=0.27, medval=3600 * 22.6,
              verboseMode=0)
    mod2 = ccfbapy.mod2irrev(iaf)
    pw = ccfbapy.cfba_moment_pw(iaf, mod2, **kw)['obj']
    assert pw == pytest.approx(0.508335, abs=1e-5)       # MATLAB MOMENT model: 0.508309
    mat = ccfbapy.solve_ccFBA_mat(ccfbapy.getccFBA_mat(iaf, mod2=mod2, **kw), verboseMode=0)['obj']
    iso = ccfbapy.ccFBA_iso(iaf, mod2=mod2, **kw)
    mr = ccfbapy.ccFBA_mr(iaf, **kw)['obj']
    assert mat <= iso['obj'] + 1e-7 <= mr + 2e-7
    ko = ccfbapy.ccFBA_iso(iaf, mat=iso['mat'], geneConcUB={'b0118': 0}, verboseMode=0)['obj']
    assert ko > 0.25                                     # acnA takes over
