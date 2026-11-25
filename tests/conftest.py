import cobra
import pandas as pd
import pytest


@pytest.fixture
def toy():
    """s -> p by R1 with the rule (A and B) or (A and C); objective: export of p."""
    m = cobra.Model('toy')
    s, p = cobra.Metabolite('s', compartment='c'), cobra.Metabolite('p', compartment='c')
    ex_s = cobra.Reaction('EX_s', lower_bound=-1000, upper_bound=1000)
    ex_s.add_metabolites({s: -1})
    r1 = cobra.Reaction('R1', lower_bound=0, upper_bound=1000)
    r1.add_metabolites({s: -1, p: 1})
    r1.gene_reaction_rule = '(A and B) or (A and C)'
    ex_p = cobra.Reaction('EX_p', lower_bound=0, upper_bound=1000)
    ex_p.add_metabolites({p: -1})
    m.add_reactions([ex_s, r1, ex_p])
    m.objective = 'EX_p'
    return m


@pytest.fixture
def toy_kw():
    """kcat 1/s (3600/h), all MW 1, budget 0.1: one unit of enzyme costs 2, so the right flux is 180."""
    return dict(Kcat=pd.DataFrame({'rxn_id': ['R_R1'], 'val': [1.0], 'dirxn': [1]}),
                MW=pd.DataFrame({'gene': ['A', 'B', 'C'], 'mw': [1.0, 1.0, 1.0]}), RHS=0.1, verboseMode=0)
