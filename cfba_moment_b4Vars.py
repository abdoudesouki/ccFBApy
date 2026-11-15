

def cfba_moment(model, Kcat=None, MW=None, selected_rxns=None, verboseMode=2, objVal=None,
                RHS=None, solver="GLPK", medval=None, runFVA=False, fvaRxn=None,
                geneConcLB=None, geneConcUB=None):
#if True:
#    Kcat=kcat.copy();MW=geneMW.copy();verboseMode=1;RHS=0.27;medval=3600*22.6
#    model=modelorg.copy();selected_rxns = None
# Kcat measurement for set of rxns (unit: 1/s)
    # MW measurement for gene using readfaa.r
    from  optlang import Variable, Constraint
    import numpy as np
    import pandas as pd
    #model=model.copy()
    geneCol = {}  # column in problem corresponding to given gene

    # Convert Kcat to 1/hr
    if Kcat is not None:
        Kcat.iloc[:,1]= [float(k) * 60 * 60 for k in Kcat.iloc[:,1]]

    # Convert MW to g/mol instead of g/mmol
    if MW is not None:
        MW.iloc[:, 1] = [float(m) for m in MW.iloc[:, 1]]

    if medval is None:
        medval = np.median(Kcat.iloc[:, 1])

    if selected_rxns is None:
        selected_rxns = [r.id for r in model.reactions if str(r.gpr) != "" and str(r.gpr) != "s0001" ]

    print("Preparing LP ...")
    prob = model.solver # .problem
    n = len(prob.variables)
    m = len(model.metabolites)

    colid = n   # Start id of gene columns
    rowind = m 
    aux_id = 0
    Vars=[]
    # Add variables for all genes
    for g in model.genes:
        if g != "s0001":
            geneCol[g.id]=(f"g_{g.id}",colid)#.append((g, colid))
            gv=Variable(f"g_{g.id}", lb=0, ub=1000)
            #print(g,gv)
            prob._add_variable(gv)
            colid += 1

    #rxnkcat_row = []
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
        tmp_kcat=Kcat[(Kcat['rxn_id']=='R_'+r) & (Kcat["dirxn"]==1)]['val']
        if len(tmp_kcat)>0:
            kval = float(tmp_kcat)
        else:
            kval = medval
        tmp_kcat=Kcat[(Kcat['rxn_id']==r) & (Kcat["dirxn"]==-1)]['val']
        if len(tmp_kcat)>0:
            rkval = float(tmp_kcat)
        else:
            rkval = medval
        vind=model.reactions.index(r)    
        lnz=1
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
               coefficients = dict()
               coefficients[model.reactions[vind].forward_variable] = lnz
               coefficients[prob.variables[colind]] = -1*kval #prob.variables[colind]._name
               constraint = Constraint(0,name=r+'_kcat',  ub=0)
               prob._add_constraint(constraint)
               prob.update()
               #print(coefficients)
               constraint.set_linear_coefficients(coefficients=coefficients)
               
               rowind=rowind+1;               
               rxnkcatRow.append({'rsn':vind, 'rxn_id' : r, 'dirxn':1,'rowind':rowind-1,'colind':colind,'val':kval});
               
               if vind_b>0:
                    coefficients = dict()
                    coefficients[model.reactions[vind].reverse_variable] = lnz
                    coefficients[prob.variables[colind]] = -1*rkval #prob.variables[colind]._name
                    constraint = Constraint(0,name=r+'_b_kcat',  ub=0)
                    prob._add_constraint(constraint)
                    prob.update()
                    #print(coefficients)
                    constraint.set_linear_coefficients(coefficients=coefficients)

                    rowind=rowind+1;                    
                    rxnkcatRow.append({'rsn':vind_b, 'rxn_id' : r, 'dirxn':-1,'rowind':rowind-1,'colind':colind,'val':rkval});
            else:#straight ANDs
#########################################2-straight ANDs----------------
                if verboseMode>2: 
                    print(f"Rxn:%s gpr: %s, straight ANDs" %(r,rl));
                #add aux variable
                #addCols(lp = problem(prob),1)
                cpxv=Variable(f"cpx_{r}", lb=0, ub=1000)
                #print(g,gv)
                prob._add_variable(cpxv)
                aux_id=colid
                colid=colid+1;
                coefficients = dict()
                coefficients[model.reactions[vind].forward_variable] = lnz
                coefficients[prob.variables[aux_id]] = -1*kval #prob.variables[colind]._name
                constraint = Constraint(0,name=r+'_kcat',  ub=0)
                prob._add_constraint(constraint)
                prob.update()
                #print(coefficients)
                constraint.set_linear_coefficients(coefficients=coefficients)
               
                rowind=rowind+1;               
                rxnkcatRow.append({'rsn':vind, 'rxn_id' : r, 'dirxn':1,'rowind':rowind-1,'colind':aux_id,'val':kval});
                if vind_b>=0:
                    coefficients = dict()
                    coefficients[model.reactions[vind].reverse_variable] = lnz
                    coefficients[prob.variables[aux_id]] = -1*rkval #prob.variables[colind]._name
                    constraint = Constraint(0,name=r+'_b_kcat',  ub=0)
                    prob._add_constraint(constraint)
                    prob.update()
                    #print(coefficients)
                    constraint.set_linear_coefficients(coefficients=coefficients)

                    rowind=rowind+1;                    
                    rxnkcatRow.append({'rsn':vind_b, 'rxn_id' : r, 'dirxn':-1,'rowind':rowind-1,'colind':aux_id,'val':rkval});

                for gene in pr[0]:
                    #colind=as.numeric(geneCol[geneCol[,"gene"]==gene,"Col"])
                    colind=geneCol[gene][1]
                    coefficients = dict()
                    coefficients[prob.variables[aux_id]] = 1
                    coefficients[prob.variables[colind]] = -1 #
                    constraint = Constraint(0,name=r+'_'+gene,  ub=0)
                    prob._add_constraint(constraint)
                    prob.update()
                    #print(coefficients)
                    constraint.set_linear_coefficients(coefficients=coefficients)

                    rowind=rowind+1;                    
                    #rxnkcatRow.append({'rsn':vind_b, 'rxn_id' : r, 'dirxn':-1,'rowind':rowind-1,'colind':aux_id,'val':rkval});

        else:#OR/AND
#########################################3-straight OR----------------
           #add row for the rxn
           row_vals=np.zeros(colid+len(pr))
           #row_vals[vind]=lnz
           row_vals_b=np.zeros(colid+len(pr))
           #if vind_b>=0: row_vals_b[vind_b]=lnz
           for  p in range(len(pr)):#ORed gene list
               if len(pr[p])==1 :
                   gene=pr[p][0]
                   colind=geneCol[gene][1]
          	    #changeMatrixRow(lp=problem(prob),i=rxn_row_id,j=colind,val=-1*kval)
                   row_vals[colind]=-1*kval
                   rxnkcatRow.append({'rsn':vind, 'rxn_id' : r, 'dirxn':1,'rowind':None,'colind':colind,'val':kval});
                   if vind_b>0: 
                        row_vals_b[colind]=-1*rkval
                        rxnkcatRow.append({'rsn':vind_b, 'rxn_id' : r, 'dirxn':-1,'rowind':None,'colind':colind,'val':rkval});
               else:#AND
		  #add auxiliary var then add rows
#########################################4-mixed OR and AND (sum of products)----------------		  
                 if verboseMode > 2: print(f"Rxn:%s gpr: %s, auxiliary ANDs"%(r,rl))
                 cpxv=Variable(f"cpx_{r}_{p}", lb=0, ub=1000)
                 prob._add_variable(cpxv)
                 aux_id=colid
                 colid=colid+1;
                 row_vals[aux_id]=-1*kval
                 rxnkcatRow.append({'rsn':vind, 'rxn_id' : r, 'dirxn':1,'rowind':None,'colind':aux_id,'val':kval});

                 if vind_b>=0: 
                     row_vals_b[aux_id]=-1*rkval   
                     rxnkcatRow.append({'rsn':vind_b, 'rxn_id' : r, 'dirxn':-1,'rowind':None,'colind':aux_id,'val':rkval});
                 for g in pr[p]:
                     colind=geneCol[g][1]
                     coefficients = dict()
                     coefficients[prob.variables[aux_id]] = 1
                     coefficients[prob.variables[colind]] = -1 #
                     constraint = Constraint(0,name='Aux_'+str(rowind), ub=0)
                     prob._add_constraint(constraint)
                     prob.update()
                     #print(coefficients)
                     constraint.set_linear_coefficients(coefficients=coefficients)
                     rowind=rowind+1; 
           #end for p
           coefficients = dict()
           L=[i for i,x in enumerate(row_vals) if x!=0]
           for i in L:
                   coefficients[prob.variables[i]] = row_vals[i]
                    
           coefficients[model.reactions[vind].forward_variable] = lnz
           constraint = Constraint(0,name=r, ub=0)
           prob._add_constraint(constraint)
           prob.update()
           #print(coefficients)
           constraint.set_linear_coefficients(coefficients=coefficients)
           rowind=rowind+1; 
           
           if vind_b>=0: 
                        coefficients = dict()
                        L=[i for i,x in enumerate(row_vals_b) if x!=0]
                        for i in L:
                           coefficients[prob.variables[i]] = row_vals_b[i]
                        
                        coefficients[model.reactions[vind].reverse_variable] = lnz
                        constraint = Constraint(0,name=r+'_b', ub=0)
                        prob._add_constraint(constraint)
                        prob.update()
                        #print(coefficients)
                        constraint.set_linear_coefficients(coefficients=coefficients)
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
        coefficients = dict()        
        for i in range(len(lcind)):
           coefficients[prob.variables[lcind[i]]] = lnz[i]

        constraint = Constraint(0,name='MW', lb=0, ub=RHS)
        prob._add_constraint(constraint)
        prob.update()
        #print(coefficients)
        constraint.set_linear_coefficients(coefficients=coefficients)
        rowind=rowind+1;  
   
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
    return({'prob':prob,'rxnFlx':rxnFlx,'obj':obj,'geneConc':geneConc,'status':sol,'rxnkcatRow':df})