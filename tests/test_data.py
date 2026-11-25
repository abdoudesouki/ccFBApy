import os

import pytest

import ccfbapy
from ccfbapy import datasets


def test_tables():
    kcat = datasets.load_kcat()
    assert list(kcat.columns) == ['rxn_id', 'val', 'dirxn'] and len(kcat) == 508
    assert len(datasets.load_mw()) == 1548
    assert {'genes', 'stoich'} <= set(datasets.load_cpx_stoich().columns)
    assert len(datasets.load_yeast_kcat()) == 938
    assert len(datasets.load_yeast_mw()) == 5899
    assert len(datasets.load_kcat_list()) > 600
    assert os.path.exists(datasets.extdata_path('NC_000913.faa'))
    with pytest.raises(FileNotFoundError):
        datasets.extdata_path('no_such_file')


def test_calc_MW():
    g = ccfbapy.calc_MW()
    assert len(g) == 4145
    assert g.loc[g['Synonym'] == 'b0002', 'mw'].iloc[0] == pytest.approx(89119.9356)


def test_readmodel():
    res = ccfbapy.readmodel(**ccfbapy.read_moment_files(), verboseMode=0)
    assert res['status'] == 'optimal'
    assert res['obj'] == pytest.approx(0.50831, abs=1e-4)
