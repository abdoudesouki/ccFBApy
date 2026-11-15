# 08/05/2021
#  gets fluxes from matrix

getccFBA_sol <- function (model, mod2=NULL, Kcat,MW=NULL,selected_rxns=NULL,verboseMode=2,objVal=NULL,
									RHS=NULL,solver=SYBIL_SETTINGS("SOLVER"),medval=NULL,
		geneConcLB=rep(0,length(allGenes(model))),geneConcUB=rep(1000,length(allGenes(model))),				
		cpx_stoich = NULL,C_mu_coef=0,mrMat=NULL,CSList=NULL,cs=1){
		#Calls function getccFBA_mat to build matrix then solves the problem 
		#Can take matrix mrMat then solve the problem
		
	if(is.null(mrMat)){
		mrMat = getccFBA_mat(model=model,mod2=mod2,Kcat=Kcat,MW=MW,verboseMode=verboseMode,RHS=RHS,cpx_stoich=cpx_stoich,
								geneConcLB=geneConcLB,geneConcUB=geneConcUB,medval=medval,selected_rxns=selected_rxns);
		
	}else{#apply geneConc bounds  
		  mrMat$clb[mrMat$geneCol[,'Col']] = geneConcLB[match(mrMat$geneCol[,'gene'],allGenes(model))]  # 
		  mrMat$cub[mrMat$geneCol[,'Col']] = geneConcUB[match(mrMat$geneCol[,'gene'],allGenes(model))]  # 

	}

    rgm=mrMat$rxnGeneCol
    cnames=mrMat$cnames
    gpr1iso=mrMat$gpr1iso

##Solve #########################
    prob=sysBiolAlg(mod2,LHS=mrMat$LHS,rlb=mrMat$rlb,rub=unlist(mrMat$rub),rtype=mrMat$rtype,lb=mrMat$clb,ub=mrMat$cub,obj=mrMat$obj_cf,
			lpdir='max',cnames=mrMat$cnames,solver=solver,  algorithm="ccFBA",#solverParm=solverParm,
					pname=sprintf("ccFBA mfe %s: %s,cpxst,kcatupd25",solver,mod_name(model)),
					writeProbToFileName=ifelse(verboseMode==3,sprintf('EC_ccFBA_Mat_%s.lp',solver),NULL));
	sol=optimizeProb(prob);
	str(sol)
	print(getMeanStatus(sol$stat,solver))
   ###
   all_flx=data.frame(stringsAsFactors=F,cs,csname=CSList[cs],rxn_id=react_id(mod2),flx=sol$fluxes[1:nr],kcat=mrMat$kcat_ir,kcat_src=mrMat$kcat_src_ir,
					enzConc=sol$fluxes[mrMat$rxnExprCol[match(1:nr,mrMat$rxnExprCol[,'flxCol']),'exprCol']],obj=obj_coef(mod2),
					ub=uppbnd(mod2))
    flx_MC = data.frame(stringsAsFactors=FALSE,rgm, #[,c('rsn','rxn_id','gene','gMW','gpr.MW')],
      flx=sol$fluxes[rgm[,'rxnCol']],
      grConc=ifelse(is.na(rgm[,'grCol']),0,sol$fluxes[rgm[,'grCol']]),
	enzConc=sol$fluxes[rgm[,'gpr.exprCol']]/as.numeric(rgm[,'stoich']),#shouldn't be NA # stoich shouldn't be 0:in case of r1g
	geneConc=ifelse(is.na(rgm[,'gCol']),0,sol$fluxes[rgm[,'gCol']]) )
	flx_MC=cbind(flx_MC,rgMC=rgm[,'gMW'] * flx_MC[,'grConc'] ,#/ as.numeric(flx_MC[,'stoich'])
			rxnMC=rgm[,'gpr.MW'] * flx_MC[,'enzConc'] ) 
	# rg_MC=flx_MC[flx_MC[,'rgMC']>1e-7 & flx_MC[,'flx']>1e-7 ,]
	rg_MC=flx_MC#25/5/2017
	rg_MC=rg_MC[order(rg_MC[,'rgMC'],decreasing=T),]
	rg_MC=cbind(rgMCorder=1:nrow(rg_MC),rg_MC)
	rg_MC=rg_MC[order(rg_MC[,'rxnMC'],decreasing=T),]
	rg_MC=cbind(rxnMCorder=1:nrow(rg_MC),rg_MC)
	
	geneConc=data.frame(gene=mrMat$geneCol[,'gene'],geneconc=sol$fluxes[mrMat$geneCol[,'Col']])
	##
	metflx_ir=sol$fluxes[1:react_num(mod2)]#Fluxes of reactions in irreversible model
	tmp_metflx=getRevFlux(model,react_id(mod2),metflx_ir)
	metflx=tmp_metflx[,"fwd"] - tmp_metflx[,"bwd"]  # fluxes of model
	
   return(list(sol=sol,metflx=metflx,metflx_ir=metflx_ir,geneConc=geneConc,rg_MC=rg_MC,flx_MC=flx_MC,mrMat=mrMat,prop=prob  ));
}