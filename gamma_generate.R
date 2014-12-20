#########################
### obtain starting values for scale accelerated failure time model
### without using artificial censoring
#########################


#########################
### obtain starting value estimates assuming
### that treatment had no effect
#########################
print("Obtain Starting Values for Gamma Parameters")
survreg.data <- as.data.frame(cbind(RESID_LIFE,EVENT_INDICATOR_USE,X.CURRENT.RECIP[,-1])[RESID_LIFE>0,])
t1 <- survreg(Surv(RESID_LIFE,EVENT_INDICATOR_USE)~ . ,data=survreg.data,dist="exponential")
sv.gamma <- t1$coef
names(sv.gamma) <- names.gamma.parm
print(sv.gamma)

#########################
### obtain estimates of gamma by solving estimating function 
### assuming a scale AFT for the SNFTM
#########################
print("Obtain Estimates of Gamma Parameters")
beta.parm <- ncol(X.CURRENT)
gamma.parm <- ncol(X.CURRENT.RECIP)

coef.vec <- c(rep(0,beta.parm),-t1$coef)
coef.vec <- matrix(coef.vec,ncol=1)
rownames(coef.vec) <- c(names.beta.parm, names.gamma.parm)
coef.total <- coef.vec

for (q in 1:10) {
LC_END_MATCH_LAS <- c(X.AT.LTX%*%coef.vec[1:beta.parm,1])
LC_MATCH_LAS <- c(X.CURRENT.RECIP%*%coef.vec[(beta.parm+1):(beta.parm+gamma.parm),1])
###
mean.func <- EVENT_INDICATOR_USE-exp(LC_MATCH_LAS)*
		(RESID_TO_TRANS+exp(LC_END_MATCH_LAS)*RESID_AFTER_TRANS)
f.func.1 <- exp(LC_MATCH_LAS)*exp(LC_END_MATCH_LAS)*RESID_AFTER_TRANS
f.func.2 <- exp(LC_MATCH_LAS)*(RESID_TO_TRANS+exp(LC_END_MATCH_LAS)*RESID_AFTER_TRANS)
estimating.eq <- c(apply(mean.func*TPT*X.CURRENT,2,sum),apply(mean.func*X.CURRENT.RECIP,2,sum))

gradient <- rbind( cbind(t(f.func.1*TPT*X.CURRENT)%*%X.AT.LTX, t(f.func.2*TPT*X.CURRENT)%*%X.CURRENT.RECIP),
	cbind(t(f.func.1*X.CURRENT.RECIP)%*%X.AT.LTX, t(f.func.2*X.CURRENT.RECIP)%*%X.CURRENT.RECIP) )

coef.vec <- solve(gradient)%*%matrix(estimating.eq,ncol=1)+coef.vec
coef.total <- cbind(coef.total,coef.vec) 
print(q)
} # end q loop
