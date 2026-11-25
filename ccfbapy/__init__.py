"""ccFBApy: cost constrained flux balance analysis (MOMENT, ccFBA) with COBRApy.

Python port of the R package sybilccFBA (Abdelmoneim Amer Desouki). Main functions:

- cfba_moment: MOMENT (Adadi et al. 2012)
- cfba_moment_pw: MOMENT with pairwise OR, on the irreversible model (as the MATLAB implementation)
- ccFBA_mr: MOMENT with multi-functional enzymes (R: cfba_moment_mr)
- getccFBA_mat / solve_ccFBA_mat / getccFBA_sol: ccFBA as a constraint matrix (cheapest isoenzyme,
  multi-functional enzymes, complex stoichiometry)
- ccFBA_iso / getccFBA_iso_mat: ccFBA with gene allocation (isoenzymes without getGpr1iso)
- simulate_EColi, simGeneOverExpr, getGpr1iso, mod2irrev, calc_MW, readmodel, addGlcTrns
- datasets: data sets of sybilccFBA (models as SBML, kcat values and molecular weights as CSV, example files)

Documentation: the reference manual (docs/, reference_manual.html) of the repository.
"""
__version__ = '0.1.0'

from .cfba_moment import cfba_moment, remove_enzyme_layer
from .ccFBA_mr import ccFBA_mr
from .cfba_moment_pw import cfba_moment_pw
from .gpr1iso import getGpr1iso
from .getccFBA_mat import getccFBA_mat, solve_ccFBA_mat, mod2irrev, IrrevModel
from .getccFBA_sol import getccFBA_sol, flx_MC_tables
from .ccFBA_iso import ccFBA_iso, getccFBA_iso_mat
from .simulate_EColi import simulate_EColi
from .simGeneOverExpr import simGeneOverExpr
from .calc_MW import calc_MW
from .readmodel import readmodel, read_moment_files
from .addGlcTrns import addGlcTrns, glc_trns_capacity
from .genFuns import scatterPlot
from . import datasets

cfba_moment_mr = ccFBA_mr  # name of the R function

__all__ = ['cfba_moment', 'remove_enzyme_layer', 'ccFBA_mr', 'cfba_moment_mr', 'cfba_moment_pw', 'getGpr1iso',
           'getccFBA_mat', 'solve_ccFBA_mat', 'mod2irrev', 'IrrevModel', 'getccFBA_sol', 'flx_MC_tables',
           'ccFBA_iso', 'getccFBA_iso_mat', 'simulate_EColi', 'simGeneOverExpr', 'calc_MW', 'readmodel',
           'read_moment_files', 'addGlcTrns', 'glc_trns_capacity', 'scatterPlot',
           'datasets', '__version__']
