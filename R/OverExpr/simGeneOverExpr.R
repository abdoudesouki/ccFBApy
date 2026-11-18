# 14/4/2021
# implemtation of over expression of one gene using ccFBA (Tin Y Pang).
# By Abdelmonem Amer Desouki

simGeneOverExpr<-function(gx,y_pct,model,mod2=NULL, Kcat,MW=NULL,selected_rxns=NULL,verboseMode=2,objVal=NULL,
		RHS=NULL,solver=SYBIL_SETTINGS("SOLVER"),medval=NULL,runFVA=FALSE,fvaRxn=NULL,FUN=cfba_moment){
		#Input two more parameters + parameters to cfba_moment
		# gx is the gene to be over expressed (gene symbol e.g b2323
		# y_pct percent of C devoted to gene x, less than 1
		# FUN: ccFBA function: cfba_moment, cfba_moment_mr, etc 
		
		##Return List
		#Sol: solution returned from minimum total fulx solution,
		#OExSol: solution returned from cfba_moment call with Over Expression of gene gx,
		#WTSol: Wild type solution before over expression,
		#OverExprEff: flag indicating the over expression affect the Objective value (e.e Biomass)
		#objVal: the value of Objective reaction of the end solution.
		#
		
		if (! gx %in% allGenes(model)){
		 stop(sprintf('simGeneOverExpr: Error: gene %s is not in the model',gx)) 
		}
		
		if(y_pct > 1)
		 stop(sprintf('simGeneOverExpr: Error: y_pct = %.3f more than 1',y_pct)) 
		 
		####
		 WTSol <- FUN(model=model,mod2=mod2,Kcat=Kcat,MW=MW,verboseMode = verboseMode,RHS=RHS,solver=solver,medval=medval)
		sol <- WTSol$sol
		geneCol <- WTSol$geneCol
		# WTSol[[1]]$ok
		gene_conc <- sol$fluxes[as.integer(geneCol[,'Col'])]
		rxn_flux <- sol$fluxes[1:react_num(mod2)];
		gxMC=gene_conc[match(gx,geneCol[,'gene'])]*MW[match(gx,MW[,1]),2]
		if(verboseMode  > 0) 
			print(sprintf('simGeneOverExpr: MC of gene %s =%.9f , %.3f %% of C by default',gx,gxMC,gxMC/RHS))
		# browser()
		if( (gxMC/RHS) >= y_pct){# gxMC > 0 &&
		  # return(list(OExSol=NULL,WTSol=WTSol,OverExprEff=FALSE))
		  OExSol=WTSol
		  OverExprEff=FALSE
		}else{
		##optimal Solution contains less than y_pct of MC of gene X 
		
			geneConcLB=rep(0,length(allGenes(model)))
			geneConcUB=rep(1000,length(allGenes(model)))
			gindx=which(allGenes(model)==gx)
			# geneConcLB[gindx]=geneConcUB[gindx]=(RHS*y_pct)/MW[match(gx,MW[,1]),2]
			geneConcLB[gindx]=(RHS*y_pct)/MW[match(gx,MW[,1]),2]#Only Lower bound Update 20/4/2021
			# RHSx=RHS*(1-y_pct)
			mom_ret2 <- FUN(model,mod2,Kcat,MW=MW,verboseMode = verboseMode,RHS=RHS,solver=solver,medval=medval,
			                             geneConcLB=geneConcLB,geneConcUB=geneConcUB)
			if(verboseMode  > 0) 
				print(sprintf('simGeneOverExpr: MC of gene %s =%.6f , %.3f %% of C by default',gx,gxMC,gxMC/RHS))
			OExSol=mom_ret2
			OverExprEff=TRUE
		}
		####### MinTotal flux : Hint: from ATMFBA  ######
		###-----------------MTF step over nonmeasured expression ----------
		if(verboseMode  > 0) print('simGeneOverExpr: calculate min total flux')
		msrdGenes=allGenes(mod2) %in% MW[MW[,2]!=0,1]
		hasGpr=(rowSums(rxnGeneMat(mod2)[,msrdGenes])>0 )  # consider GPR with at least one measured gene
		# hasGpr[!gpr_rxn_ind]=FALSE;
		if(verboseMode  > 0) 
		    print(sprintf('simGeneOverExpr: min total flux #Rxns:%d, to be affected=%d ',length(hasGpr),sum(!hasGpr)))
	
		OExFlux=OExSol$sol$fluxes[1:react_num(mod2)]
		mod2_1=mod2
		lowbnd(mod2_1)[ hasGpr ] = OExFlux[ hasGpr ] 
		uppbnd(mod2_1)[ hasGpr ] = OExFlux[ hasGpr ]

		lowbnd(mod2_1)[obj_coef(mod2_1)!=0] = OExFlux[obj_coef(mod2_1)!=0] ;
		uppbnd(mod2_1)[obj_coef(mod2_1)!=0] = OExFlux[obj_coef(mod2_1)!=0] ;
		
		
		mtf=optimizeProb(mod2_1,algorithm="mtf", solver=solver)#,solverParm=solverParm

		return(list(Sol=mtf,OExSol=OExSol,WTSol=WTSol,OverExprEff=OverExprEff,objVal=OExFlux[obj_coef(mod2_1)!=0]))
		
 }
 
#####################

 
# Test units: test_simGeneOverExpr.R