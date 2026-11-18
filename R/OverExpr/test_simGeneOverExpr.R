#16/4/2021
# Test units

library(sybilccFBA)
data(kcat)
data(mw)
data(iAF1260)
model= iAF1260
mod2=mod2irrev(model)
setwd('D:\\Dropbox\\pprs\\ppr2_moment\\')
source('Tin\\simGeneOverExpr.R')
source('sybilccFBA\\R\\cfba_moment.R')
source('sybilccFBA\\R\\cfba_moment_mr.R')
source('sybilccFBA\\R\\getccFBA_mat.R')
source('sybilccFBA\\R\\getccFBA_sol.R')

# read the new kcat file
kcat2_tmp = read.table('Tin\\kcat2')
names(kcat2_tmp) <- c("rxn_id","val","dirxn");
for (qqq in 1:length(kcat$rxn_id)) {
	ii2 = which(kcat2_tmp$rxn_id == kcat$rxn_id[qqq]);
	ii3 = which(kcat2_tmp$dirxn == kcat$dirxn[qqq]);
	ii = intersect(ii2,ii3);
	if (length(ii)==0) {
		tmpp <- data.frame(kcat$rxn_id[qqq],kcat$val[qqq],kcat$dirxn[qqq]);
		names(tmpp) <- c("rxn_id","val","dirxn");
		kcat2_tmp <- rbind(kcat2_tmp,tmpp);
	}
}



uppbnd(mod2)[react_id(mod2)=="R_EX_glc_e__b"]=1000
uppbnd(mod2)[react_id(mod2)=="R_EX_glyc_e__b"]=0
uppbnd(mod2)[react_id(mod2)=="R_EX_ac_e__b"]=0
uppbnd(mod2)[react_id(mod2)=="R_EX_o2_e__b"]=1000
lowbnd(mod2)[react_id(mod2)=="R_ATPM"]=0

C_mu=0.27
# gx='b0929' ; pct=1e-04
gx='b2323'# 15% , 'b0728'
v_Fun=cfba_moment
allRes=NULL
for(pct in c(0.05,0.1,0.15,0.25,0.5)){#,1.1 different percentages
	print(pct)
	# RR=simGeneOverExpr(gx,pct,model,mod2,kcat,MW=mw,verboseMode = 1,RHS=C_mu,solver="glpkAPI",medval=3600*22.6)
	RR=simGeneOverExpr(gx,pct,model,mod2,kcat2_tmp,MW=mw,verboseMode = 1,RHS=C_mu,solver="glpkAPI",medval=3600*22.6,FUN=v_Fun)

	mtfFlx=getFluxDist(RR$Sol)
	OExFlx=RR$OExSol[[1]]$fluxes[1:react_num(mod2)]
	WTFlx=RR$WTSol[[1]]$fluxes[1:react_num(mod2)]
	
	sol <- RR$OExSol[[1]]
	geneCol <- RR$OExSol[[2]]
	gene_conc_OEx <- sol$fluxes[as.integer(geneCol[,'Col'])]
	gxMWOEx = gene_conc_OEx[match(gx,geneCol[,'gene'])]*mw[match(gx,mw[,1]),2]
	WTSol <- RR$WTSol[[1]]
	geneCol <- RR$WTSol[[2]]
	 gene_conc <- WTSol$fluxes[as.integer(geneCol[,'Col'])]
	gxMWWT=gene_conc[match(gx,geneCol[,'gene'])]*mw[match(gx,mw[,1]),2]

	allRes=rbind(allRes,
       data.frame(pct=pct,absMtfFlx=sum(abs(mtfFlx)), absOExFlx=sum(abs(OExFlx)),
             absWTFlx=sum(abs(WTFlx)), OverExprEff=RR$OverExprEff,Bm=RR$objVal,gxMWWT=gxMWWT,gxMWOEx=gxMWOEx,MCpctOExUsed=gxMWOEx/C_mu))
}

print(allRes)


######################## Display list of genes' Molecular Crawding (MC) and percentage in decending order
MCv=(gene_conc*mw[match(geneCol[,'gene'],mw[,1]),2])
 ix=order(-MCv)
 Res=data.frame(mw[match(geneCol[,'gene'],mw[,1]),],MCv,mcvpct=MCv/C_mu)
head( Res[ix,])
############ Gene Conc
# RR=simGeneOverExpr(gx,0.20,model,mod2,kcat,MW=mw,verboseMode = 1,RHS=0.27,solver="glpkAPI",medval=3600*22.6)


library(sybilccFBA)
data(kcat)
setwd('D:\\Dropbox\\pprs\\ppr2_moment\\')
source('Tin\\simGeneOverExpr.R')
source('sybilccFBA\\R\\cfba_moment.R')
source('sybilccFBA\\R\\cfba_moment_mr.R')
source('sybilccFBA\\R\\getccFBA_mat.R')
source('sybilccFBA\\R\\getccFBA_sol.R')

# read the new kcat file
kcat2_tmp = read.table('Tin\\kcat2')
names(kcat2_tmp) <- c("rxn_id","val","dirxn");
for (qqq in 1:length(kcat$rxn_id)) {
	ii2 = which(kcat2_tmp$rxn_id == kcat$rxn_id[qqq]);
	ii3 = which(kcat2_tmp$dirxn == kcat$dirxn[qqq]);
	ii = intersect(ii2,ii3);
	if (length(ii)==0) {
		tmpp <- data.frame(kcat$rxn_id[qqq],kcat$val[qqq],kcat$dirxn[qqq]);
		names(tmpp) <- c("rxn_id","val","dirxn");
		kcat2_tmp <- rbind(kcat2_tmp,tmpp);
	}
}

data(iAF1260)
model= iAF1260
data(mw)
data(kcat)
mod2=mod2irrev(model)

uppbnd(mod2)[react_id(mod2)=="R_EX_glc_e__b"]=1000
uppbnd(mod2)[react_id(mod2)=="R_EX_glyc_e__b"]=0
uppbnd(mod2)[react_id(mod2)=="R_EX_ac_e__b"]=0
uppbnd(mod2)[react_id(mod2)=="R_EX_o2_e__b"]=1000
lowbnd(mod2)[react_id(mod2)=="R_ATPM"]=0
kcat2_tmp$src='Tin'
C_mu=0.27
gx='b0929' ; pct=1e-04
v_Fun=cfba_moment_mr#
for(v_Fun in c("cfba_moment_mr","cfba_moment","getccFBA_sol")){#
print(paste(v_Fun,pct))
	# RR=simGeneOverExpr(gx,pct,model,mod2,kcat,MW=mw,verboseMode = 1,RHS=C_mu,solver="glpkAPI",medval=3600*22.6)
	RR=simGeneOverExpr(gx,pct,model,mod2,kcat2_tmp,MW=mw,verboseMode = 1,RHS=C_mu,solver="glpkAPI",medval=3600*22.6,FUN=get(v_Fun))
	mtfFlx=getFluxDist(RR$Sol)
	OExFlx=RR$OExSol$sol$fluxes[1:react_num(mod2)]
	WTFlx=RR$WTSol$sol$fluxes[1:react_num(mod2)]
	
	sol <- RR$OExSol$sol
	geneCol <- RR$OExSol$geneCol
	gene_conc_OEx <- sol$fluxes[as.integer(geneCol[,'Col'])]
	gxMWOEx = gene_conc_OEx[match(gx,geneCol[,'gene'])]*mw[match(gx,mw[,1]),2]
	WTSol <- RR$WTSol$sol
	geneCol <- RR$WTSol$geneCol
	 gene_conc <- WTSol$fluxes[as.integer(geneCol[,'Col'])]
	gxMWWT=gene_conc[match(gx,geneCol[,'gene'])]*mw[match(gx,mw[,1]),2]

	print(
       data.frame(pct=pct,v_Fun=v_Fun,absMtfFlx=sum(abs(mtfFlx)), absOExFlx=sum(abs(OExFlx)),
             absWTFlx=sum(abs(WTFlx)), OverExprEff=RR$OverExprEff,Bmwt=WTFlx[obj_coef(mod2)!=0],
			 Bm=RR$objVal,gxMWWT=gxMWWT,gxMWOEx=gxMWOEx,MCpctOExUsed=gxMWOEx/C_mu))
}