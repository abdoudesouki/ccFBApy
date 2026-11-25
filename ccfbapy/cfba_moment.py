
# 4/2/2024
# by Abdelmoneim Amer Desouki
# part of ccFBApy

def cfba_moment(model, Kcat=None, MW=None, selected_rxns=None, verboseMode=2, objVal=None,
                RHS=None, solver="GLPK", medval=None, #runFVA=False, fvaRxn=None,
                geneConcLB=None, geneConcUB=None):

    # Kcat measurement for set of rxns (unit: 1/s) will be converted to (1/hr)
    # MW measurement for gene using readfaa.r
    # selected_rxns: reactions in the enzymes constraint, when it is None all reactions with gpr are used
    # medval: value of kcat for reactions with kcat given, when it is not specified median of known values is used
    # RHS: the upper limit on molecular crawding used by enzymes (e.g: 0.27 for Ecoli from Adadi et al 2011)
    # solver is GLPK.
    # geneConcLB, geneConcUB: Add gene concentration lower/upper bounds, dictionary of geneid to max concentration
        # can be used to simulate over expression of a gene.
    # vi : mmol/gDW/hr  , rxn catal by gi with Kcati
	# Kcat: hr-1
	# gi: mmol/gDW
	# MWi: g/mmol (val of AA in g/mol, divide by 1000)
  
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

    print("Preparing LP ...")
    # Work on a copy so the caller's model is not modified: the enzyme variables (g_*, cpx_*)
    # and kcat/MW constraints are added directly to the solver problem. Modifying the input model
    # made later copies of it (e.g. modelR=modelorg.copy()) already contain these variables,
    # and the next call failed with ContainerAlreadyContains ('g_b2307').
    model = model.copy()
    prob = model.solver # .problem
    remove_enzyme_layer(model)
    n = len(prob.variables)
    m = len(model.metabolites)

    colid = n   # Start id of gene columns
    rowind = m 
    aux_id = 0
    Vars={}
    Cons={}
    # Add variables for all genes
    for g in model.genes:
        if g != "s0001":
            geneCol[g.id]=(f"g_{g.id}",colid)#.append((g, colid))
            Vars[colid]={'name':f"g_{g.id}", 'lb':0, 'ub':1000}
            colid += 1

    if geneConcLB is not None:
        for g in geneConcLB.keys():
            Vars[geneCol[g][1]]['lb']=geneConcLB[g]
            
    if geneConcUB is not None:
        for g in geneConcUB.keys():
            Vars[geneCol[g][1]]['ub']=geneConcUB[g]
            
    for r in selected_rxns:
        rl = str(model.reactions.get_by_id(r).gpr)# gpr(model)[react_id(model) == r]
        # Remove s0001 gene from rules in iAF model
        if rl == "( b0875  or  s0001 )":
            rl = "b0875"
        if rl == "( b1377  or  b0929  or  b2215  or  b0241  or  s0001  or  b3875  or  b1319  or  b0957 )":
            # Add your desired logic here
            pass
        # Append to the list r_xnkcat_row
        #rxnkcat_row.append(rl)

#ng=length(allGenes(model))
    def convert_gpr(rl):
        if rl == '( b0875  or  s0001 )':
            return 'b0875'
        elif rl == '( b1377  or  b0929  or  b2215  or  b0241  or  s0001  or  b3875  or  b1319  or  b0957 )':
            return '( b1377  or  b0929  or  b2215  or  b0241 or  b3875  or  b1319  or  b0957 )'
        elif rl == '( s0001  or  b0451 )':
            return '(b0451)'
        elif rl == '( s0001  or  b3927 )':
            return '( b3927 )'
        elif rl == '(b3927 or s0001)':
            return 'b3927'
        elif rl == '(s0001 or b0957 or b3875 or b2215 or b0241 or b1319 or b1377 or b0929)':
            return '(b0957 or b3875 or b2215 or b0241 or b1319 or b1377 or b0929)'
        elif rl == '(s0001 or b0875)':
            return 'b0875'
        elif rl == '(b0451 or s0001)':
            return 'b0451'
        elif rl == '(b1377 or b0241 or b3875 or s0001 or b0957 or b2215 or b1319 or b0929)':
            return '(b1377 or b0241 or b3875 or b0957 or b2215 or b1319 or b0929)'
        elif rl == '(b1377 or b2215 or b1319 or b3875 or s0001 or b0957 or b0241 or b0929)':
            return '(b1377 or b2215 or b1319 or b3875 or b0957 or b0241 or b0929)'
        elif rl == '(s0001 or b0451)':
            return 'b0451'
        elif rl == '(s0001 or b3927)':
            return 'b3927'
        elif rl == '(s0001 or b1779)':
            return 'b1779'
        elif rl == '(b1779 or s0001)':
            return 'b1779'
        elif rl == '(b0875 or s0001)':
            return 'b0875'
        else:
            return rl.replace(")", " ) ").replace("(", " ( ").replace(" OR ", " or ").replace(" AND ", " and ")

    # Main loop
    rxnkcatRow=[]
    
    for r in selected_rxns:
        #print(prob.optimize(),'b4 rxn:'+r)
        rl = str(model.reactions.get_by_id(r).gpr)
        rl = convert_gpr(rl)
        #======
        # Kcat ids are SBML ids with 'R_' prefix (as in the R version, where sybil keeps it), but cobra
        # strips the prefix when reading SBML; add it unless the reaction id already has it.
        kcat_id = r if r.startswith('R_') else 'R_'+r
        tmp_kcat=Kcat[(Kcat['rxn_id']==kcat_id) & (Kcat["dirxn"]==1)]['val']
        if len(tmp_kcat)>0:
            kval = float(tmp_kcat.iloc[0])
        else:
            kval = medval
        tmp_kcat=Kcat[(Kcat['rxn_id']==kcat_id) & (Kcat["dirxn"]==-1)]['val']
        if len(tmp_kcat)>0:
            rkval = float(tmp_kcat.iloc[0])#Warning
        else:
            rkval = medval
        vind=model.reactions.index(r)    
        lnz=1
        # vind_b: index of reversible rxn, -1 when irreversible. R used vind_b>0 since R indices
        # start at 1; in python index 0 is valid, so all checks must be vind_b>=0.
        if not model.reactions.get_by_id(r).reversibility:            
            vind_b=-1            
        else:# Add row with median value to rev direction
           vind_b=model.reactions.index(r)
        if verboseMode>2: print("Process:", r,rl,vind,vind_b)    
        #======
        pr = []
        for x in rl.split(" or "):
            yl=[]
            for y in x.split(" and "):
                yl.append(y.replace("(", "").replace(")", "").replace(" ", ""))
            pr.append(yl)
        if len(pr)==1:# no OR (only one term) 
            #add row vi<=gi*kcat
            if len(pr[0]) == 1:
               #########################################1-one Gene----------------
               gene=pr[0][0]
               #get the corresponding col of the gene
               colind=geneCol[gene][1]
               #AddRow Vind - Kcat*colind <=0
               Cons[rowind]={'name':r+'_kcat','ub':0,
                        'lcind':[model.reactions[vind].forward_variable._index-1,colind],
                        'lnz':[lnz,-1*kval]}
               
               rowind=rowind+1;               
               rxnkcatRow.append({'rsn':vind, 'rxn_id' : r, 'dirxn':1,'rowind':rowind-1,'colind':colind,'val':kval});
               
               if vind_b>=0:
                    
                    Cons[rowind]={'name':r+'_b_kcat','ub':0,
                        'lcind':[model.reactions[vind].reverse_variable._index-1,colind],
                        'lnz':[lnz,-1*rkval]}
 
                    rowind=rowind+1;                    
                    rxnkcatRow.append({'rsn':vind_b, 'rxn_id' : r, 'dirxn':-1,'rowind':rowind-1,'colind':colind,'val':rkval});
            else:#straight ANDs
#########################################2-straight ANDs----------------
                if verboseMode>2: 
                    print(f"Rxn:%s gpr: %s, straight ANDs" %(r,rl));
                #add aux variable
                Vars[colid]={'name':f"cpx_{r}", 'lb':0, 'ub':1000}
                aux_id=colid
                colid=colid+1;
                Cons[rowind]={'name':r+'_kcat','ub':0,
                        'lcind':[model.reactions[vind].forward_variable._index-1,aux_id],
                        'lnz':[lnz,-1*kval]}
              
                rowind=rowind+1;               
                rxnkcatRow.append({'rsn':vind, 'rxn_id' : r, 'dirxn':1,'rowind':rowind-1,'colind':aux_id,'val':kval});
                if vind_b>=0:
                    Cons[rowind]={'name':r+'_b_kcat','ub':0,
                        'lcind':[model.reactions[vind].reverse_variable._index-1,aux_id],
                        'lnz':[lnz,-1*rkval]}
                    rowind=rowind+1;                    
                    rxnkcatRow.append({'rsn':vind_b, 'rxn_id' : r, 'dirxn':-1,'rowind':rowind-1,'colind':aux_id,'val':rkval});

                for gene in pr[0]:
                    #colind=as.numeric(geneCol[geneCol[,"gene"]==gene,"Col"])
                    colind=geneCol[gene][1]
                    Cons[rowind]={'name':r+'_'+gene,'ub':0, 'lcind':[aux_id,colind],  'lnz':[1,-1]}
                    rowind=rowind+1;                    
                    #rxnkcatRow.append({'rsn':vind_b, 'rxn_id' : r, 'dirxn':-1,'rowind':rowind-1,'colind':aux_id,'val':rkval});

        else:#OR/AND
#########################################3-straight OR----------------
           #add row for the rxn
           row_vals={}#np.zeros(colid+len(pr))
           #row_vals[vind]=lnz
           row_vals_b={}#np.zeros(colid+len(pr))
           #if vind_b>=0: row_vals_b[vind_b]=lnz
           for  p in range(len(pr)):#ORed gene list
               if len(pr[p])==1 :
                   gene=pr[p][0]
                   colind=geneCol[gene][1]
          	    #changeMatrixRow(lp=problem(prob),i=rxn_row_id,j=colind,val=-1*kval)
                   row_vals[colind]=-1*kval
                   rxnkcatRow.append({'rsn':vind, 'rxn_id' : r, 'dirxn':1,'rowind':None,'colind':colind,'val':kval});
                   if vind_b>=0: 
                        row_vals_b[colind]=-1*rkval
                        rxnkcatRow.append({'rsn':vind_b, 'rxn_id' : r, 'dirxn':-1,'rowind':None,'colind':colind,'val':rkval});
               else:#AND
		  #add auxiliary var then add rows
#########################################4-mixed OR and AND (sum of products)----------------		  
                 if verboseMode > 2: print(f"Rxn:%s gpr: %s, auxiliary ANDs"%(r,rl))
                 #cpxv=Variable(f"cpx_{r}_{p}", lb=0, ub=1000)
                 #prob._add_variable(cpxv)
                 Vars[colid]={'name':f"cpx_{r}_{p}", 'lb':0, 'ub':1000}
                 aux_id=colid
                 colid=colid+1;
                 row_vals[aux_id]=-1*kval
                 rxnkcatRow.append({'rsn':vind, 'rxn_id' : r, 'dirxn':1,'rowind':None,'colind':aux_id,'val':kval});

                 if vind_b>=0: 
                     row_vals_b[aux_id]=-1*rkval   
                     rxnkcatRow.append({'rsn':vind_b, 'rxn_id' : r, 'dirxn':-1,'rowind':None,'colind':aux_id,'val':rkval});
                 for g in pr[p]:
                     colind=geneCol[g][1]
                     Cons[rowind]={'name':'Aux_'+str(rowind),'ub':0,'lcind':[aux_id,colind], 'lnz':[1,-1]}
                     rowind=rowind+1; 
           row_vals[model.reactions[vind].forward_variable._index-1]= lnz
           if vind_b>=0: # as in R: row_vals_b[vind_b]=lnz (reverse flux must be in the '_b' row)
                row_vals_b[model.reactions[vind].reverse_variable._index-1]= lnz
           
           Cons[rowind]={'name':r,'ub':0, 'lcind':row_vals.keys(), 'lnz':row_vals.values()}
           rowind=rowind+1; 
           
           if vind_b>=0: 
                        Cons[rowind]={'name':r+'_b','ub':0, 'lcind':row_vals_b.keys(), 'lnz':row_vals_b.values()}
                        rowind=rowind+1;  
			            #rxnkcatRow[is.na(rxnkcatRow[,'rowind']) & rxnkcatRow[,'rxn_id']==r & rxnkcatRow[,'dirxn']==-1,'rowind']=rowind-1;
    ###############II. Add MW constraint------------------
## Sum(gi * MWi) <= C[gdwpr/gDW]
    
    if MW is not None:
      lnz=[]
      lcind=[]
      for g in geneCol.keys():
            colind=geneCol[g][1]
            #changeColsBnds(lp = problem(prob),colind,lb=geneConcLB[gIX[g]],ub=geneConcUB[gIX[g]])#16/4/2021 Update
            if len(MW[MW.iloc[:,0]==g])>0 :# consider only genes with computed MW.
                lnz.append(MW[MW.iloc[:,0]==g].iloc[:,1].values[0])
                if verboseMode > 3: print(g,lnz[-1])
                lcind.append(colind)
      if RHS is not None: #test value of upper bound (C)        
        Cons[rowind]={'name':'MW','ub':RHS, 'lcind':lcind, 'lnz':lnz}
        rowind=rowind+1;  
   
    #Add Variables
    if verboseMode > 0: print('Adding variables...')
    LVars=[]
    for i in range(n,colid):
        LVars.append(Variable(Vars[i]['name'], lb=Vars[i]['lb'], ub=Vars[i]['ub']))
    prob._add_variables(LVars)
    #Add Constraints
    if verboseMode > 0: print('Adding constraints...')
    #LCons=[]
    for i in range(m,rowind):
        coefficients = dict()     
        lcind=list(Cons[i]['lcind'])        
        lnz=list(Cons[i]['lnz'])        
        for j in range(len(lcind)):
          coefficients[prob.variables[lcind[j]]] = lnz[j]

        constraint = Constraint(0,name=Cons[i]['name'],  ub=Cons[i]['ub'])
        prob._add_constraint(constraint)
        #print(i)3rd time 3397
        # prob.update()
        constraint.set_linear_coefficients(coefficients=coefficients)
    #print('update prob..')
    prob.update()
    
    #get fluxes
    print("Solving LP ...")
    sol=prob.optimize()
    flx=prob.primal_values
    geneConc={str(g):flx[geneCol[str(g)][0]] for g in model.genes}
    #print(geneConc)
    rxnFlx={r.id:flx[r.forward_variable._name]-flx[r.reverse_variable._name] for r in model.reactions}
    #print(rxnFlx)
    obj=prob.objective.value
    
    df = pd.DataFrame(rxnkcatRow)         
    return({'prob':prob,'model':model,'rxnFlx':rxnFlx,'obj':obj,'geneConc':geneConc,'status':sol,'rxnkcatRow':df})


def remove_enzyme_layer(model):
    """Remove variables/constraints added by a previous call of cfba_moment or ccFBA_mr on this model
    (or on the model it was copied from), so that they can be added again cleanly.
    Variables: g_<gene>, gr_<gene>_<rxn>, cpx_*; constraints: <rxn>, <rxn>_b, <rxn>_kcat, <rxn>_b_kcat,
    <rxn>_<gene>, Aux_*, mfe_<gene>, MW. Reaction variables and metabolite constraints are never removed."""
    rxn_vars = set()
    for r in model.reactions:
        rxn_vars.add(r.forward_variable.name)
        rxn_vars.add(r.reverse_variable.name)
    gene_vars = {f"g_{g.id}" for g in model.genes}
    old_vars = [v for v in model.solver.variables if v.name not in rxn_vars and
                (v.name in gene_vars or v.name.startswith("cpx_") or v.name.startswith("gr_"))]

    met_ids = {m.id for m in model.metabolites}
    cons_names = {"MW"}
    cons_names.update(f"mfe_{g.id}" for g in model.genes)
    for r in model.reactions:
        cons_names.update([r.id, r.id + "_b", r.id + "_kcat", r.id + "_b_kcat"])
        cons_names.update(r.id + "_" + g.id for g in r.genes)
    old_cons = [c for c in model.solver.constraints if c.name not in met_ids and
                (c.name in cons_names or c.name.startswith("Aux_"))]

    if old_vars or old_cons:
        print("Removing %d variables and %d constraints left from a previous enzyme-constrained run"
              % (len(old_vars), len(old_cons)))
        model.solver.remove(old_cons + old_vars)
        model.solver.update()
