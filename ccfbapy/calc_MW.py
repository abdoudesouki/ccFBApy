# Calculate molecular weights of proteins from the genome (.faa) file.
# by Abdelmoneim Amer Desouki
# part of ccFBApy, Python version of calc_MW (sybilccFBA, R/calc_MW.R)
import os

import numpy as np
import pandas as pd

# example files of the R package (R: system.file("extdata", package="sybilccFBA")), shipped with ccfbapy
EXTDATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'extdata')

# amino acids counted, in the order of the columns of the result (as in R)
AALST = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "W", "V", "Y"]
WATER_MW = 18.015  # g/mol


def _find(fname):
    """File name as given, or the file of the same name in the extdata folder of the package."""
    if os.path.exists(fname) or os.path.isabs(fname):
        return fname
    alt = os.path.join(EXTDATA, fname)
    return alt if os.path.exists(alt) else fname


def _read_faa(fname, aa_mw):
    """Proteins of one .faa file: list of (gsn, gi, mw, counts)."""
    rows = []
    gsn, gi, counts = 0, None, None

    def add():
        mw = sum(counts[ch] * (aa_mw[ch] - WATER_MW) for ch in AALST) + WATER_MW  # water of the first AA kept
        rows.append((gsn, gi, mw, [counts[ch] for ch in AALST]))

    with open(fname) as con:
        for ln in con:
            ln = ln.rstrip('\r\n')
            if ln.startswith('>'):
                if gsn > 0:
                    add()
                gi = ln.split('|')[1]
                counts = dict.fromkeys(AALST, 0)
                gsn += 1
            elif gsn > 0:
                for ch in AALST:  # counts of each AA
                    counts[ch] += ln.count(ch)
    if gsn > 0:
        add()
    return rows


def calc_MW(aa_fname="aa.txt", ptt_fname="test2.ptt", faa_fname="NC_000913.faa", nchrm=None):
    """Molecular weights of the proteins of a genome, from the protein sequences (.faa).

    aa_fname:  file with the amino acids (tab separated, columns AA_1 (one letter code) and MW (g/mol))
    ptt_fname: protein table (tab separated, with the column PID: protein id, and the gene names, e.g. Synonym)
    faa_fname: protein sequences (FASTA, header '>gi|<PID>|...'), one file per chromosome: a file name or a list
    nchrm:     number of chromosomes (number of .faa files used), default: all files in faa_fname
    File names that do not exist are looked up in the extdata folder of the package
    (ccfbapy.datasets.extdata_path(); the defaults are the E. coli
    files there).

    Molecular weight = sum of the weights of the residues - water (18.015 g/mol) * (number of residues - 1).
    Only the 20 standard amino acids are counted.

    Returns a DataFrame: the columns of the .ptt file merged (on PID) with gsn (serial number of the protein in
    its .faa file), mw (molecular weight in g/mol) and the number of each amino acid (A, C, ..., Y), sorted by PID.
    The ccFBA functions need the molecular weight in g/mmol: divide mw by 1000.
    """
    # Default units of Amino Acids and Water molecule are g/mol divide by 1000 to
    # convert to g/mmol the match default flux units mmol/gDW/hr
    aa = pd.read_csv(_find(aa_fname), sep='\t')
    aa_mw = dict(zip(aa['AA_1'], aa['MW'].astype(float)))
    missing = [ch for ch in AALST if ch not in aa_mw]
    if missing:
        raise ValueError('amino acids missing in %s: %s' % (aa_fname, ', '.join(missing)))
    ptt = pd.read_csv(_find(ptt_fname), sep='\t')

    faa = [faa_fname] if isinstance(faa_fname, str) else list(faa_fname)
    if nchrm is None:
        nchrm = len(faa)
    if nchrm > len(faa):
        raise ValueError('nchrm=%d but %d .faa file(s) given' % (nchrm, len(faa)))

    rows = []
    for chrm in range(nchrm):
        rows += _read_faa(_find(faa[chrm]), aa_mw)
    cnts = pd.DataFrame({'gi': [float(r[1]) for r in rows], 'gsn': [r[0] for r in rows],
                         'mw': [r[2] for r in rows]})
    cnts = pd.concat([cnts, pd.DataFrame([r[3] for r in rows], columns=AALST)], axis=1)

    ptt = ptt.copy()
    ptt['PID'] = pd.to_numeric(ptt['PID'], errors='coerce')
    geneCnt = ptt.merge(cnts, left_on='PID', right_on='gi').drop(columns='gi')
    geneCnt = geneCnt[['PID'] + [c for c in geneCnt.columns if c != 'PID']]  # as R's merge: key first
    geneCnt = geneCnt.sort_values('PID', kind='stable').reset_index(drop=True)
    if np.issubdtype(geneCnt['PID'].dtype, np.floating) and (geneCnt['PID'] % 1 == 0).all():
        geneCnt['PID'] = geneCnt['PID'].astype('int64')
    return geneCnt
