"""Prepare the data files of the package ccfbapy (ccfbapy/data). Development tool, run from the root folder of the
repository before building the package:

    python tools/make_package_data.py            # all steps
    python tools/make_package_data.py --no-rdata # without the export of the R data sets (no rdata needed)

Steps:
1. copy the CSV tables of data/ to ccfbapy/data/ (Ecoli_kcat.csv, Ecoli_MW.csv, yeast_imm_kcat.csv, yeast_mw.csv);
2. copy the example files of the R package (R/sybilccFBA/inst/extdata) to ccfbapy/data/extdata/;
3. compress the SBML models of models/ (BiGG ids) to ccfbapy/data/bigg/<model>.xml.gz;
4. export the data sets of the R package (R/sybilccFBA/data/*.RData) to ccfbapy/data/sybilccFBA/:
   the models iAF1260, iJO1366, iML1515, iMM904 as SBML (<model>.xml.gz, with the ids of sybil) and the tables
   kcat, mw, mw_iML1515, yst_kcat, yst_mw as CSV. Each export is read again and compared with the R data set.
   This is the only step that needs the package rdata (pip install rdata); ccfbapy itself reads the SBML and CSV
   files and does not need rdata.

In the R data sets of iAF1260 and iMM904 the metabolites have no compartment (met_comp is NA, ids such as
'M_glc_D_c[NA]'). SBML needs one: it is taken from the id (the letter before '[', e.g. 'c'); the ids are kept.
SBML writes GPR rules without redundant parentheses ('(A and B) and C' -> 'A and B and C': the same rule), and a
missing metabolite charge is read back as 0.
"""
import argparse
import gzip
import os
import re
import shutil
import sys
import warnings

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKG_DATA = os.path.join(ROOT, 'ccfbapy', 'data')
MODELS = ['iAF1260', 'iJO1366', 'iML1515', 'iMM904']
TABLES = ['kcat', 'mw', 'mw_iML1515', 'yst_kcat', 'yst_mw']
COMPARTMENTS = {'c': 'cytosol', 'e': 'extracellular', 'p': 'periplasm', 'm': 'mitochondrion', 'x': 'peroxisome',
                'r': 'endoplasmic reticulum', 'n': 'nucleus', 'v': 'vacuole', 'g': 'golgi apparatus'}


def copy_tables():
    for f in ['Ecoli_kcat.csv', 'Ecoli_MW.csv', 'yeast_imm_kcat.csv', 'yeast_mw.csv']:
        shutil.copy2(os.path.join(ROOT, 'data', f), os.path.join(PKG_DATA, f))
    src = os.path.join(ROOT, 'R', 'sybilccFBA', 'inst', 'extdata')
    dst = os.path.join(PKG_DATA, 'extdata')
    os.makedirs(dst, exist_ok=True)
    for f in sorted(os.listdir(src)):
        shutil.copy2(os.path.join(src, f), os.path.join(dst, f))
    print('copied the tables of data/ and %d files of inst/extdata' % len(os.listdir(src)))


def compress_bigg():
    dst = os.path.join(PKG_DATA, 'bigg')
    os.makedirs(dst, exist_ok=True)
    for m in MODELS:
        src = os.path.join(ROOT, 'models', m + '.xml')
        out = os.path.join(dst, m + '.xml.gz')
        with open(src, 'rb') as fi, gzip.GzipFile(out, 'wb', compresslevel=9, mtime=0) as fo:
            shutil.copyfileobj(fi, fo)
        print('%s -> %s (%d kB)' % (os.path.relpath(src, ROOT), os.path.relpath(out, ROOT), os.path.getsize(out) // 1024))


def _gpr_terms(gpr):
    """GPR rule -> sorted list of its AND terms (sorted gene lists), i.e. the rule without redundant parentheses."""
    sys.path.insert(0, ROOT)
    from ccfbapy.gpr1iso import _gpr_terms as terms
    return sorted(terms(gpr))


def _compare(a, b):
    """Differences between two cobra models (ids, order, bounds, stoichiometry, GPR, objective, names)."""
    bad = []
    if [r.id for r in a.reactions] != [r.id for r in b.reactions]:
        bad.append('reaction ids')
    if [m.id for m in a.metabolites] != [m.id for m in b.metabolites]:
        bad.append('metabolite ids')
    if [g.id for g in a.genes] != [g.id for g in b.genes]:
        bad.append('genes')
    if bad:
        return bad
    for ra, rb in zip(a.reactions, b.reactions):
        if ra.bounds != rb.bounds:
            bad.append('bounds of ' + ra.id)
        if {m.id: v for m, v in ra.metabolites.items()} != {m.id: v for m, v in rb.metabolites.items()}:
            bad.append('stoichiometry of ' + ra.id)
        if _gpr_terms(ra.gpr) != _gpr_terms(rb.gpr):  # SBML removes redundant parentheses: compare the terms
            bad.append('GPR of ' + ra.id)
        if ra.name != rb.name or ra.objective_coefficient != rb.objective_coefficient:
            bad.append('name/objective of ' + ra.id)
    for ma, mb in zip(a.metabolites, b.metabolites):
        # a missing charge (None in the R data set of iAF1260) is read back from SBML as 0
        if (ma.name, ma.formula, ma.charge or 0) != (mb.name, mb.formula, mb.charge or 0):
            bad.append('metabolite ' + ma.id)
    return bad[:10]


def export_rdata():
    import cobra
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sybil2cobra import read_RData
    src = os.path.join(ROOT, 'R', 'sybilccFBA', 'data')
    dst = os.path.join(PKG_DATA, 'sybilccFBA')
    os.makedirs(dst, exist_ok=True)
    for name in MODELS:
        model = read_RData(os.path.join(src, name + '.RData'))[name]
        fixed = 0
        for m in model.metabolites:
            if m.compartment is None or m.compartment not in COMPARTMENTS:
                mm = re.search(r'_([a-z])\[[^\]]*\]$', m.id) or re.search(r'\[([a-z])\]$', m.id)
                m.compartment = mm.group(1) if mm else 'c'
                fixed += 1
        model.compartments = {c: COMPARTMENTS.get(c, c) for c in sorted({m.compartment for m in model.metabolites})}
        out = os.path.join(dst, name + '.xml')
        cobra.io.write_sbml_model(model, out)
        with open(out, 'rb') as fi, gzip.GzipFile(out + '.gz', 'wb', compresslevel=9, mtime=0) as fo:
            shutil.copyfileobj(fi, fo)
        os.remove(out)
        back = cobra.io.read_sbml_model(out + '.gz')
        diff = _compare(model, back)
        if diff:
            raise RuntimeError('%s: the exported SBML differs from the R data set: %s' % (name, diff))
        fa, fb = model.slim_optimize(), back.slim_optimize()
        if not np.isclose(fa, fb):
            raise RuntimeError('%s: FBA %s (R data set) != %s (SBML)' % (name, fa, fb))
        print('%s.RData -> %s (%d kB): %d reactions, %d metabolites, %d genes, FBA %.6f%s' % (
            name, os.path.relpath(out + '.gz', ROOT), os.path.getsize(out + '.gz') // 1024, len(back.reactions),
            len(back.metabolites), len(back.genes), fb,
            ', compartment of %d metabolites taken from the id' % fixed if fixed else ''))
    for name in TABLES:
        obj = read_RData(os.path.join(src, name + '.RData'))[name]
        if name == 'yst_kcat':  # numeric matrix in R
            obj = pd.DataFrame(np.asarray(obj.values, dtype=float), columns=['rxn_id', 'val', 'dirxn'])
            obj['rxn_id'] = obj['rxn_id'].astype(int)
            obj['dirxn'] = obj['dirxn'].astype(int)
        obj = obj.reset_index(drop=True)
        for c in obj.columns:
            if str(obj[c].dtype) in ('Int32', 'Int64') and not obj[c].isna().any():
                obj[c] = obj[c].astype(int)
            elif str(obj[c].dtype).startswith('string'):
                obj[c] = obj[c].astype(object)
        out = os.path.join(dst, name + '.csv')
        obj.to_csv(out, index=False)  # floats are written with full precision
        back = pd.read_csv(out, float_precision='round_trip')  # exact (the default parser can change the last bit)
        for c in obj.columns:
            a, b = obj[c], back[c]
            ok = np.allclose(a.astype(float), b.astype(float), rtol=0, atol=0, equal_nan=True) \
                if np.issubdtype(b.dtype, np.number) else (a.astype(str).values == b.astype(str).values).all()
            if not ok:
                raise RuntimeError('%s: column %s differs after the export' % (name, c))
        print('%s.RData -> %s: %d rows, columns %s' % (name, os.path.relpath(out, ROOT), len(back), list(back.columns)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--no-rdata', action='store_true', help='do not export the R data sets (needs rdata)')
    args = ap.parse_args()
    warnings.filterwarnings('ignore', message='Missing constructor for R class')
    os.makedirs(PKG_DATA, exist_ok=True)
    copy_tables()
    compress_bigg()
    if not args.no_rdata:
        export_rdata()


if __name__ == '__main__':
    main()
