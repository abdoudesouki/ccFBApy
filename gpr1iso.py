#get molecular weight of enzyme complexes and the isoenzyme with minimum MW.
# By Dr. Abdelmoneim Amer Desouki
# 31/3/2024 
import pandas as pd
import numpy as np
def getGpr1iso (model,MW,cpx_stoich=None):
#if True:
    cpx_stoich=None
    nOR=[len(str(r.gpr).split(' or '))-1 for r in model.reactions]
    nAND = [len(str(r.gpr).split(' and '))-1 for r in model.reactions]
    
    ### 0. Preprocessing evaluate complex MW's
    #single gene rule
    r1gflg =[len(str(model.reactions[i].gpr).strip()) > 0 and nOR[i]==0 and nAND[i]==0 
                        for i in range(len(model.reactions))]
    mod_cpx=[]
    r1gix=[i for i,x in enumerate(r1gflg) if x]
    for i,x in enumerate(r1gix):
         d={'rsn':x, 'rxn':model.reactions[x].id,'bool_expr':str(model.reactions[x].gpr),
         'genes':str(model.reactions[x].gpr).strip(),'termNo':1,'nOR':nOR[x],'nAND':nAND[x],'gprtype':'r1g'}
         mod_cpx.append(d)
    #one complex     
    r1cpxflg=[nOR[i]==0 and nAND[i]>0 for i in range(len(model.reactions))]
    r1cpxix=[i for i,x in enumerate(r1cpxflg) if x]
    for i,x in enumerate(r1cpxix):
         rl=str(model.reactions[x].gpr).replace("(", "").replace(")", "")
         r1cpx_genes=','.join(sorted(rl.strip().split(' and ')))
         d={'rsn':x, 'rxn':model.reactions[x].id,'bool_expr':str(model.reactions[x].gpr),
         'genes':r1cpx_genes,'termNo':1,'nOR':nOR[x],'nAND':nAND[x],'gprtype':'r1cpx'}
         mod_cpx.append(d)
    #iso enzymes     
    rmgflg=[nOR[i]>0 and nAND[i]==0 for i in range(len(model.reactions))]
    rmgix=[i for i,x in enumerate(rmgflg) if x]
    for i,x in enumerate(rmgix):
         #termNo:used when equal MW
         #rmg_genes=','.join(sorted(str(model.reactions[x].gpr).strip().split(' or ')))
        rmg_genes=str(model.reactions[x].gpr).strip().split(' or ')
        #print(rmg_genes)
        for j,g in enumerate(rmg_genes):
             d={'rsn':x, 'rxn':model.reactions[x].id,'bool_expr':str(model.reactions[x].gpr),
             'genes':g,'termNo':i,'nOR':nOR[x],'nAND':nAND[x],'gprtype':'rmg'}
             mod_cpx.append(d)
    #multiple complex     
    rmcpxflg = [nOR[i]>0 and nAND[i]>0 for i in range(len(model.reactions))]
    rmcpxix=[i for i,x in enumerate(rmcpxflg) if x]
    for i,x in enumerate(rmcpxix):
         rl=str(model.reactions[x].gpr)
         # Replace ')' with ' ) ' and '(' with ' ( '
         pr = []
         for j in rl.split(" or "):
            yl=[]
            for y in j.split(" and "):
                yl.append(y.replace("(", "").replace(")", "").replace(" ", ""))
            pr.append(yl)
         pr_genes=[','.join(sorted([g for g in p])) for p in pr]
         for j,z in enumerate(pr_genes):
             d={'rsn':x, 'rxn':model.reactions[x].id,'bool_expr':str(model.reactions[x].gpr),
             'genes':pr_genes[j],'termNo':j,'nOR':nOR[x],'nAND':nAND[x],'gprtype':'rmcpx'}
             mod_cpx.append(d)
    
    mod_cpx=pd.DataFrame(mod_cpx)
    #print(mod_cpx['gprtype'].value_counts)
    #test
    if cpx_stoich is not None:
        mod_cpx_st = mod_cpx.copy()
        for cpx in mod_cpx_st['genes']:
            if cpx in list(cpx_stoich.loc[:,'genes']):
                mod_cpx_st.loc[mod_cpx_st['genes']==cpx,'stoich'] = cpx_stoich.loc[cpx_stoich['genes']==cpx, 'stoich']
    else:
        mod_cpx_st = mod_cpx.copy()
        mod_cpx_st['stoich'] = ''#','.join(['1' for x in mod_cpx_st['genes']])
    
    
    # Set row names of MW dataframe
    MW.index = MW.iloc[:, 0]

    # Initialize lists for storing results
    mod_cpx_MW1s = []  # All ones
    mod_cpx_MWtxt = []
    mod_st_MW = []  # Using stoich

    # Iterate over rows in mod_cpx_st
    for i in range(mod_cpx_st.shape[0]):
        genes = mod_cpx_st.loc[i, 'genes']
        stoich = str(mod_cpx_st.loc[i, 'stoich'])
        
        # Split genes by comma and get corresponding molecular weights
        mws = [] 
        for gene in genes.split(','):
            try:
                mws.append(MW.loc[gene, 'mw'])
            except KeyError:
                print(i,gene, 'gene has no molecular weight, will be set to 0')
                mws.append(0)
        
        # Calculate total molecular weight for this complex
        mod_cpx_MWtxt.append(','.join(map(str, mws)))
        mod_cpx_MW1s.append(sum(mws))
        if stoich!='': 
            mod_st_MW.append(sum(mw * float(s) for mw, s in zip(mws, stoich.split(','))))
        else:
            mod_st_MW.append(sum(mw  for mw in mws))

    # Calculate termNG
    termNG = [len(gene) - len(gene.replace(',', '')) + 1 for gene in mod_cpx_st['genes']]
    
    # Create a dataframe for mod_cpx_mw
    mod_cpx_mw = pd.DataFrame({
        'genes': mod_cpx_st['genes'],
        # 'stoich': mod_cpx_st['stoich'],
        'nOR':mod_cpx_st['nOR'],
        'rsn':mod_cpx_st['rsn'],
        'rxn':mod_cpx_st['rxn'],
        'termNo':mod_cpx_st['termNo'],
        'nAND':mod_cpx_st['nAND'],
        'gprtype':mod_cpx_st['gprtype'],
        'termNG': termNG,
        'mod_cpx_MW1s': mod_cpx_MW1s,
        'mod_st_MW': mod_st_MW,
        'mod_cpx_MWtxt': mod_cpx_MWtxt
    }, columns=['genes', 'nOR', 'rsn','rxn','termNo','nAND','gprtype','termNG', 'mod_cpx_MW1s', 'mod_st_MW', 'mod_cpx_MWtxt'])

    # Calculate MW based on stoich (if available)
    mod_cpx_mw['MW'] = mod_cpx_mw.apply(lambda row: mod_cpx_MW1s if pd.isna(row['mod_st_MW']) else row['mod_st_MW'], axis=1)

    # Merge with cpx_stoich (if available)
    if cpx_stoich is not None:
        mod_cpx_mw = pd.merge(mod_cpx_mw, cpx_stoich, on='genes', how='left')
    else:
        mod_cpx_mw['stoich']= mod_cpx_st['stoich']
### Preprocessing 2. choose best iso enzyme for each rxn  --------------------------------	
    react_id=[r.id for r in model.reactions];
    allGenes=[str(g) for g in model.genes]			
# Filter rows where nOR > 0
    flg = mod_cpx_mw['nOR'] > 0
    #rxn_iso_cpx = mod_cpx_mw[flg].groupby('rsn')['MW'].min()
    mod_cpx_mw['MW_min'] = mod_cpx_mw.groupby('rsn')['MW'].transform("min")

    # Choose iso gpr with min(MW)
    gpr1iso = mod_cpx_mw[(mod_cpx_mw['nOR'] == 0) | 
                         (mod_cpx_mw['MW_min']==mod_cpx_mw['MW'])]
    # Handle equal iso enzymes
    if gpr1iso['rsn'].duplicated().any():
        #tmpxx=[i for i,x in enumerate(gpr1iso['rsn'].duplicated()) if x ]
        #for x in tmpxx: print(gpr1iso[gpr1iso['rsn'] ==x])
        gpr1iso['equal_iso'] = gpr1iso.groupby('rsn')['termNo'].transform("min")#[gpr1iso['nOR'] > 0]
        gpr1iso = gpr1iso[(gpr1iso['nOR'] == 0) | (gpr1iso['termNo']==gpr1iso['equal_iso'])]
                          #(gpr1iso.apply(lambda row: f"{row['rsn']}{row['termNo']}" in equal_iso, axis=1))]
    # Calculate termNG (if needed)
    # gpr1iso['termNG'] = gpr1iso['genes'].apply(lambda x: len(x) - len(x.replace(',', '')) + 1)

    # Set stoich to 1 where it's NA and termNG is 1
    gpr1iso.loc[(gpr1iso['stoich']=='') & (gpr1iso['termNG'] == 1), 'stoich'] = 1

    # Initialize rgst matrix
    rgst = np.zeros((len(react_id), len(allGenes)))
    rgst = pd.DataFrame(rgst, index=react_id, columns=allGenes)

    # Fill rgst based on termNG
    for i in gpr1iso.index:
        if gpr1iso.loc[i, 'termNG'] == 1:
            if gpr1iso.loc[i, 'stoich']=='':#is.na in R
                rgst.loc[gpr1iso.loc[i, 'rxn'], gpr1iso.loc[i, 'genes'].split(',')] = 1.0
            else:
                rgst.loc[gpr1iso.loc[i, 'rxn'], gpr1iso.loc[i, 'genes'].split(',')] = float(gpr1iso.loc[i, 'stoich'])
        else:#TermNG>1
            if gpr1iso.loc[i, 'stoich']!='':
                rgst.loc[gpr1iso.loc[i, 'rxn'], gpr1iso.loc[i, 'genes'].split(',')] = \
                    [float(s) for s in str(gpr1iso.loc[i, 'stoich']).split(',')]
            else:#NA and termNG>1
                rgst.loc[gpr1iso.loc[i, 'rxn'], gpr1iso.loc[i, 'genes'].split(',')] = 1.0
                #gpr1iso[i,'stoich']=paste(rep(1,gpr1iso[i,'termNG']),collapse=',')
                gpr1iso.loc[i, 'stoich']=','.join(['1' for x in range(gpr1iso.loc[i, 'termNG'])])
                
    return({'gpr1iso':gpr1iso,'mod_cpx_mw':mod_cpx_mw,'rgst':rgst})
