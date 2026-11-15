

# 15/2/2024
# by Abdelmoneim Amer Desouki
# part of ccFBApy

#given  Kcat vector (in 1/Sec), and molecular weights  (Maximize biomass)
#1-Add gi vars (gene concentrations)
#2-add gpr constraints: add f/b for reversible rxns react_rev
#		vi/kcati <= gi
#       vi cat by g1 & g2 :  add aux variable then Aux <= all(gi), vi<= Kcat*Aux
#       vi cat by g1 | g2 :   OR is eliminated by gpr1iso, choose the smallest enzyme
#      complex gpr recursively
##... ----------
#3-	mr:multiple reactions (multiple function enzymes
# 		for genes catalyzing more n rxn(n>1): add gij: j 1:n  such that Sum(gij) <= gi
#4-add Enzimatic constraint Sum(MW*gi)<=C

def getccFBA_mat (model, Kcat=None, MW=None, selected_rxns=None, verboseMode=2, objVal=None,
                RHS=None, solver="GLPK", medval=None, #runFVA=False, fvaRxn=None,
                geneConcLB=None, geneConcUB=None,cpx_stoich = None,C_mu_coef=0):				
		
#Multiple OR are not done pairwise but evaluated in one line in contrast to MOMENT original implementation.
#    this can affect knockouts of isoenzymes
#addCols
#	Kcat measurement for set of rxns  (Unit: 1/s)
#	MW measurement for gene using readfaa.r, calc_MW or gene annotation
#   RHS : C  g protein/gDW  (metabolic)
#   objVal: when set the minimum C that gives this growth, will be returned as objective
#   cpx_stoich: giving the stoichiometry of complexes: data frame containing at least two columns 'genes','stoich' 1/10/2015
##               'gene'  : delimited string of genes(ordered by geneid(e.g bnumber)),
##               'stoich': number of subunits of each gene in the same order in 'genes'
##  C_mu_coef: used to have C as a linear function of mu (biomass) : C = RHS + C_mu_coef*Biomass

#			units: mmol/gDW/h
	# vi : mmol/gDW/hr  , rxn catal by gi with Kcati
	# Kcat: hr-1
	# gi: mmol/gDW
	# MWi: g/mmol (val of AA in g/mol, divide by 1000)
	# Calc_MW
	# C : g/gDW protein (metabolic)
    from  optlang import Variable, Constraint
    import numpy as np
    import pandas as pd
    
    geneCol = {}  # column in problem corresponding to given gene

    # Convert Kcat to 1/hr
    if Kcat is not None:
        Kcat=Kcat.copy()
        Kcat.iloc[:,1]= [float(k) * 60 * 60 for k in Kcat.iloc[:,1]]

    # Convert MW to g/mol instead of g/mmol
    if MW is not None:
        MW=MW.copy()
        MW.iloc[:, 1] = [float(m) for m in MW.iloc[:, 1]]

    if medval is None:
        medval = np.median(Kcat.iloc[:, 1])
    
    if selected_rxns is None:
        selected_rxns = [r.id for r in model.reactions if str(r.gpr) != "" and str(r.gpr) != "s0001" ]
    
    print("Preprocessing isoenzymes ... ");
	r1i=getGpr1iso(model,MW=MW,cpx_stoich=cpx_stoich)
	gpr1iso=r1i['gpr1iso']
	mod_cpx_mw=r1i['mod_cpx_mw']; rgst=r1i['rgst'];
    if(any(isnan(gpr1iso[,'stoich']))){
        stop("There are NA values in gpr1iso")
    }
### start making constraint matrix ---------------------------------------------------
	print("Preparing LP ... ");
