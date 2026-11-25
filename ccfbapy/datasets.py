"""Data sets shipped with ccFBApy (SBML and CSV files; no R or rdata needed).

- ``load_dataset(name)``: the data sets of the R package sybilccFBA (R: ``data(name)``): the models ``iAF1260``,
  ``iJO1366``, ``iML1515``, ``iMM904`` (``cobra.Model``, ids as in sybil) and the tables ``kcat``, ``mw``,
  ``mw_iML1515``, ``yst_kcat``, ``yst_mw`` (``pandas.DataFrame``). They were exported from the ``.RData`` files of
  sybilccFBA to SBML and CSV with ``tools/make_package_data.py`` of the repository.
- ``load_model(name, source='sybilccFBA')``: one of the four models; ``source='bigg'``: the SBML models of the
  repository folder ``models/`` (BiGG ids, e.g. ``EX_glc__D_e``).
- ``load_kcat()``, ``load_mw()``: the E. coli kcat values (1/s) and molecular weights (g/mmol) in the format of the
  ``Kcat`` and ``MW`` arguments (``data/Ecoli_kcat.csv``, ``data/Ecoli_MW.csv`` of the repository).
- ``load_yeast_kcat()``, ``load_yeast_mw()``: the same for S. cerevisiae.
- ``load_cpx_stoich()``: protein complex stoichiometry (``cpx_stoich_me.csv``).
- ``load_kcat_list(model=None)``: kcat_max of ``allKcats_upd34_dd_h.csv`` as a Kcat table (as in the example of
  ``simulate_EColi``).
- ``model_path(name, source)``, ``extdata_path(fname=None)``, ``data_path(fname=None)``: paths of the files.
"""
import os

import pandas as pd

_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

MODELS = ('iAF1260', 'iJO1366', 'iML1515', 'iMM904')
TABLES = ('kcat', 'mw', 'mw_iML1515', 'yst_kcat', 'yst_mw')
DATASETS = MODELS + TABLES
SOURCES = ('sybilccFBA', 'bigg')


def _read_csv(path, **kw):
    # round_trip: the values are read exactly as written (the default parser can change the last digit)
    return pd.read_csv(path, float_precision='round_trip', **kw)


def data_path(fname=None):
    """Folder of the data files of the package, or the path of one of its files."""
    return _DATA if fname is None else os.path.join(_DATA, fname)


def extdata_path(fname=None):
    """Folder of the example files of sybilccFBA (inst/extdata), or the path of one of them, e.g.
    extdata_path('cpx_stoich_me.csv'). The list of files: os.listdir(extdata_path())."""
    d = os.path.join(_DATA, 'extdata')
    if fname is None:
        return d
    p = os.path.join(d, fname)
    if not os.path.exists(p):
        raise FileNotFoundError("'%s' is not in %s (files: %s)" % (fname, d, ', '.join(sorted(os.listdir(d)))))
    return p


def list_datasets():
    """Names of the data sets of sybilccFBA that load_dataset can read."""
    return list(DATASETS)


def model_path(name, source='sybilccFBA'):
    """Path of the SBML file (gzip) of a model: source 'sybilccFBA' (the R data set, ids as in sybil) or 'bigg'
    (the model of the folder models/ of the repository, BiGG ids)."""
    if name not in MODELS:
        raise ValueError("unknown model '%s', available: %s" % (name, ', '.join(MODELS)))
    if source not in SOURCES:
        raise ValueError("unknown source '%s', available: %s" % (source, ', '.join(SOURCES)))
    return os.path.join(_DATA, source, name + '.xml.gz')


def load_model(name, source='sybilccFBA'):
    """Metabolic model as a cobra.Model: iAF1260, iJO1366, iML1515 (E. coli) or iMM904 (S. cerevisiae).

    source='sybilccFBA' (default): the data set of the R package (R: data(name)), with its ids (e.g.
    'R_EX_glc_e_' in iAF1260, 'EX_glc__D(e)' in iML1515), bounds and GPR rules, so that the results are the
    same as in R. source='bigg': the SBML model of the folder models/ of the repository (BiGG ids, e.g.
    'EX_glc__D_e'); it is another version of the model (other bounds, some other ids)."""
    import cobra
    return cobra.io.read_sbml_model(model_path(name, source))


def load_dataset(name):
    """Data set of sybilccFBA (as data(name) in R): the models (cobra.Model, see load_model) and the tables
    kcat, mw, mw_iML1515, yst_kcat (columns rxn_id, val, dirxn), yst_mw (DataFrame)."""
    if name not in DATASETS:
        raise ValueError("unknown data set '%s', available: %s" % (name, ', '.join(DATASETS)))
    if name in MODELS:
        return load_model(name)
    df = _read_csv(os.path.join(_DATA, 'sybilccFBA', name + '.csv'))
    if name == 'yst_mw':
        df['LENGTH'] = df['LENGTH'].astype('Int64')
    return df


def load_kcat():
    """E. coli kcat values used in MOMENT (1/s): DataFrame rxn_id (with the prefix 'R_'), val, dirxn."""
    return _read_csv(os.path.join(_DATA, 'Ecoli_kcat.csv'))


def load_mw():
    """E. coli molecular weights (g/mmol) of 1548 genes: DataFrame Synonym (gene), mw."""
    return _read_csv(os.path.join(_DATA, 'Ecoli_MW.csv'))


def load_yeast_kcat():
    """S. cerevisiae kcat values (1/s): DataFrame rxn_id (a number, see the documentation of yst_kcat), val,
    dirxn."""
    return _read_csv(os.path.join(_DATA, 'yeast_imm_kcat.csv'))


def load_yeast_mw():
    """S. cerevisiae molecular weights (g/mmol): DataFrame gene, mw, LENGTH."""
    return _read_csv(os.path.join(_DATA, 'yeast_mw.csv'))


def load_cpx_stoich():
    """Protein complex stoichiometry (cpx_stoich_me.csv): DataFrame with the columns genes and stoich (and
    Enzyme.ID, Gene.composition, ...), for the cpx_stoich argument."""
    return pd.read_csv(extdata_path('cpx_stoich_me.csv'))


def load_kcat_list(model=None, column='kcat_max'):
    """kcat values of allKcats_upd34_dd_h.csv as a Kcat table (rxn_id, val, dirxn, src), as in the example of
    simulate_EColi: rows with an ijo_id, val from `column`, missing src set to 'Max'. When model is given, only
    its reactions are kept."""
    kl = _read_csv(extdata_path('allKcats_upd34_dd_h.csv'))
    kl = kl[kl['ijo_id'].notna()]
    kcat = pd.DataFrame({'rxn_id': kl['ijo_id'], 'val': kl[column], 'dirxn': kl['dirxn'], 'src': kl['src']})
    if model is not None:
        kcat = kcat[kcat['rxn_id'].isin([r.id for r in model.reactions])]
    kcat = kcat.copy()
    kcat.loc[kcat['src'].isna(), 'src'] = 'Max'
    return kcat.reset_index(drop=True)
