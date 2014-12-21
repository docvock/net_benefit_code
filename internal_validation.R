#####################################
#####################################
### model internal validation
####################################
####################################


########################################
### Plot transformed residual survival of person transplanted to those 
### that have not been transplanted. These curves should overlap with the same 
### recipient characteristics. Here we only adjuste for LAS and not recipient age
### native disease grouping, or center volume
### In these plots, use full follow-up information
######################################### 
coef.beta.early <-  matrix(c(coef.beta[1:length(change.comp),]), ncol=1)
coef.beta.diff <-  matrix(c(coef.beta[-(1:length(change.comp)),])[change.comp], ncol=1)	
coef.beta.const <- matrix(c(coef.beta[-(1:length(change.comp)),])[-change.comp], ncol=1)	

LC_END_MATCH_LAS <- c(X.AT.LTX.CONST%*%coef.beta.const)
LC_END_MATCH_LAS.1 <- c(X.AT.LTX.TV%*%coef.beta.early)
LC_END_MATCH_LAS.2 <- c(X.AT.LTX.TV%*%coef.beta.diff)

LC_PLUS_MATCH_LAS <- c(X.MAX.FUTURE.CONST%*%coef.beta.const)
LC_PLUS_MATCH_LAS.1 <- c(X.MAX.FUTURE.TV%*%coef.beta.early)
LC_PLUS_MATCH_LAS.2 <- c(X.MAX.FUTURE.TV%*%coef.beta.diff)

LC_PLUS_MATCH_LAS1_ALT <- c(X.ALT.CONST%*%coef.beta.const)
LC_PLUS_MATCH_LAS1_ALT.1 <- c(X.ALT.TV%*%coef.beta.early)
LC_PLUS_MATCH_LAS1_ALT.2 <- c(X.ALT.TV%*%coef.beta.diff)

RESID_LIFE_SCALE <- RESID_TO_TRANS+RESID_AFTER_TRANS*exp(LC_END_MATCH_LAS)*exp(LC_END_MATCH_LAS.1)*
	exp((LC_END_MATCH_LAS.2)*RESID_AFTER_TRANS/(RESID_AFTER_TRANS+180))
POT_CENSOR_C <- pmin(POT_CENSOR*exp(LC_PLUS_MATCH_LAS)*exp(LC_PLUS_MATCH_LAS.1)*exp((LC_PLUS_MATCH_LAS.2)*POT_CENSOR/(POT_CENSOR+180)),
	(POT_CENSOR))

POT_CENSOR_1 <- RESID_LIFE_SCALE
TEST1 <- ifelse(POT_CENSOR_1 < POT_CENSOR_C & EVENT_INDICATOR_USE==0,1,0)

pdf("graphics/internal_validation.pdf",width=10,height=6)

par(mfrow=c(2,3))
for (j in 1:6) {
X.matrix.plus <- cbind(1,0,(LAS.mean[j]+5-30),(LAS.mean[j]+5-30)^2,0,0,0,0,0,0,0,0)

LC_PLUS_MATCH_LAS.MAX <- c(X.matrix.plus[,3:12]%*%coef.beta.phase.2[3:12,]) 
TIME_MAX1 <- pmin(exp(LC_PLUS_MATCH_LAS.MAX)*exp(coef.beta[1])*
			exp(coef.beta[3]*2340/(2340+180))*2340,2340)


T1A <- EVENT_INDICATOR_USE[TRANSPLANTED==1&MATCH_LAS>=LAS.lower[j]&MATCH_LAS<LAS.upper[j]&TEST1==0]
T1B <- RESID_LIFE_SCALE[TRANSPLANTED==1&MATCH_LAS>=LAS.lower[j]&MATCH_LAS<LAS.upper[j] &TEST1==0]

T2A <- EVENT_INDICATOR_USE[TRANSPLANTED==0&MATCH_LAS>=LAS.lower[j]&MATCH_LAS<LAS.upper[j]&TEST1==0]
T2B <- RESID_LIFE_SCALE[TRANSPLANTED==0&MATCH_LAS>=LAS.lower[j]&MATCH_LAS<LAS.upper[j]&TEST1==0]

f1A <- survfit(Surv(T1B, T1A) ~ 1)
f2A <- survfit(Surv(T2B, T2A) ~ 1)

plot(f1A,xlab="Days",ylab="Survival Proportion",main=paste("LAS ",LAS.lower[j],"-",LAS.upper[j],sep=""),mark.time=F,conf.int=F,col="blue",xlim=c(0,TIME_MAX1),
ylim=c(0.55,1.0),lwd=3,cex.lab=1.5,cex.main=2)
lines(f2A,mark.time=F,conf.int=F,col="red",lty=2,lwd=3)
}

dev.off()

