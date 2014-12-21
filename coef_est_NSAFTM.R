#####################################################
#####################################################
### non-scale (time-varying) accelerated failure time model
### estimation using artifical censoring
#####################################################
#####################################################

coef.beta.init <- c(coef.beta[change.comp,1],coef.beta[,1])
##### after several iteration re-define second phase coefficients
coef.beta.init[c(3,4)] <- coef.beta.init[c(3,4)]-coef.beta.init[c(1,2)]
coef.beta <- matrix(coef.beta.init,ncol=1)


coef.total.SNNSAFTM <- coef.beta
rownames(coef.beta.SNNSAFTM) <- names.beta.parm.tv

for (q in 1:35) {
start <- proc.time()
# values of coefficients of the structural model initially for those that are time-varying
coef.beta.early <-  matrix(c(coef.beta[1:length(change.comp),]), ncol=1)
# difference in values of coefficients of the structural model long-term - initial
# for those that are time-varying
coef.beta.diff <-  matrix(c(coef.beta[-(1:length(change.comp)),])[change.comp], ncol=1)	
# values of coefficient in the structural model which are not time-varying
coef.beta.const <- matrix(c(coef.beta[-(1:length(change.comp)),])[-change.comp], ncol=1)	

####

#################
### get needed linear combinations of data
#################
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

POT_CENSOR_C <- pmin(POT_CENSOR*exp(LC_PLUS_MATCH_LAS)*exp(LC_PLUS_MATCH_LAS.1)*exp((LC_PLUS_MATCH_LAS.2)*
		POT_CENSOR/(POT_CENSOR+180)),(POT_CENSOR))
POT_CENSOR_C_IND <- ifelse(POT_CENSOR*exp(LC_PLUS_MATCH_LAS)*exp(LC_PLUS_MATCH_LAS.1)*exp((LC_PLUS_MATCH_LAS.2)*
		POT_CENSOR/(POT_CENSOR+180)) < (POT_CENSOR),1,0)

POT_CENSOR_ALT_C <- pmin(c.adj*exp(LC_PLUS_MATCH_LAS1_ALT)*exp(LC_PLUS_MATCH_LAS1_ALT.1)*exp((LC_PLUS_MATCH_LAS1_ALT.2)*
		c.adj/(c.adj+180)),POT_CENSOR_C) 
POT_CENSOR_ALT_C_IND <- ifelse(c.adj*exp(LC_PLUS_MATCH_LAS1_ALT)*exp(LC_PLUS_MATCH_LAS1_ALT.1)*exp((LC_PLUS_MATCH_LAS1_ALT.2)*
		c.adj/(c.adj+180))<c.adj,1,0)
POT_CENSOR_ALT_C_IND2 <- ifelse(POT_CENSOR_C<(c.adj*exp(LC_PLUS_MATCH_LAS1_ALT)*exp(LC_PLUS_MATCH_LAS1_ALT.1)*
		exp((LC_PLUS_MATCH_LAS1_ALT.2)*c.adj/(c.adj+180))),1,0)

POT_CENSOR_C_CENSOR <- ifelse((POT_CENSOR_C>(RESID_LIFE_SCALE)&(EVENT_INDICATOR_USE==0)),1,0)
POT_CENSOR_ALT_C_CENSOR <- ifelse((POT_CENSOR_ALT_C>(RESID_LIFE_SCALE)&(EVENT_INDICATOR_USE==0)),1,0)
POT_CENSOR_C <- ifelse((POT_CENSOR_C>(RESID_LIFE_SCALE)&(EVENT_INDICATOR_USE==0)),RESID_LIFE_SCALE,POT_CENSOR_C)
POT_CENSOR_ALT_C <- ifelse((POT_CENSOR_ALT_C>(RESID_LIFE_SCALE)&(EVENT_INDICATOR_USE==0)),RESID_LIFE_SCALE,POT_CENSOR_ALT_C)


mean.func.1 <- smooth.indicator(RESID_LIFE_SCALE,POT_CENSOR_ALT_C,eps)-exp(LC_MATCH_LAS)*
        smooth.min(RESID_LIFE_SCALE,POT_CENSOR_ALT_C,eps)
mean.func.2 <- smooth.indicator(RESID_LIFE_SCALE,POT_CENSOR_C,eps2)-exp(LC_MATCH_LAS)*
				smooth.min(RESID_LIFE_SCALE,POT_CENSOR_C,eps2)

estimating.eq <- c(apply(mean.func.1*TPT*X.CURRENT.TV,2,sum),apply(mean.func.2*TPT*X.CURRENT,2,sum))

f.func.2 <- RESID_AFTER_TRANS*exp(LC_END_MATCH_LAS)*exp(LC_END_MATCH_LAS.1)*
		exp((LC_END_MATCH_LAS.2)*RESID_AFTER_TRANS/(RESID_AFTER_TRANS+180))

gradient.x.1 <- t((smooth.indicator.deriv.x(RESID_LIFE_SCALE,POT_CENSOR_ALT_C,eps)-
  exp(LC_MATCH_LAS)*smooth.min.deriv.x(RESID_LIFE_SCALE,POT_CENSOR_ALT_C,eps))*TPT*X.CURRENT.TV)%*%
		cbind((f.func.2)*X.AT.LTX.TV,
				(f.func.2)*RESID_AFTER_TRANS/(RESID_AFTER_TRANS+180)*X.AT.LTX.TV,
				(f.func.2)*X.AT.LTX.CONST)
gradient.x.2 <- t((smooth.indicator.deriv.x(RESID_LIFE_SCALE,POT_CENSOR_C,eps2)-
	exp(LC_MATCH_LAS)*smooth.min.deriv.x(RESID_LIFE_SCALE,POT_CENSOR_C,eps2))*TPT*X.CURRENT)%*%
		cbind((f.func.2)*X.AT.LTX.TV,
				(f.func.2)*RESID_AFTER_TRANS/(RESID_AFTER_TRANS+180)*X.AT.LTX.TV,
				(f.func.2)*X.AT.LTX.CONST)
gradient.x <- rbind(gradient.x.1,gradient.x.2)

gradient.c.2 <- t((smooth.indicator.deriv.c(RESID_LIFE_SCALE,POT_CENSOR_C,eps2)-
  exp(LC_MATCH_LAS)*smooth.min.deriv.c(RESID_LIFE_SCALE,POT_CENSOR_C,eps2))*TPT*X.CURRENT)%*%
	(cbind(POT_CENSOR*exp(LC_PLUS_MATCH_LAS)*exp(LC_PLUS_MATCH_LAS.1)*exp((LC_PLUS_MATCH_LAS.2)*POT_CENSOR/(POT_CENSOR+180))*X.MAX.FUTURE.TV,
		POT_CENSOR*exp(LC_PLUS_MATCH_LAS)*exp(LC_PLUS_MATCH_LAS.1)*exp((LC_PLUS_MATCH_LAS.2)*POT_CENSOR/(POT_CENSOR+180))*
			POT_CENSOR/(POT_CENSOR+180)*X.MAX.FUTURE.TV,
		POT_CENSOR*exp(LC_PLUS_MATCH_LAS)*exp(LC_PLUS_MATCH_LAS.1)*exp((LC_PLUS_MATCH_LAS.2)*POT_CENSOR/(POT_CENSOR+180))*
			X.MAX.FUTURE.CONST)*POT_CENSOR_C_IND*(1-POT_CENSOR_C_CENSOR)+	 
	cbind((f.func.2)*X.AT.LTX.TV,(f.func.2)*RESID_AFTER_TRANS/(RESID_AFTER_TRANS+180)*X.AT.LTX.TV,(f.func.2)*X.AT.LTX.CONST)*POT_CENSOR_C_CENSOR)

gradient.c.1 <- t((smooth.indicator.deriv.c(RESID_LIFE_SCALE,POT_CENSOR_ALT_C,eps)-
	exp(LC_MATCH_LAS)*smooth.min.deriv.c(RESID_LIFE_SCALE,POT_CENSOR_ALT_C,eps))*TPT*X.CURRENT.TV)%*%
	((cbind(c.adj*exp(LC_PLUS_MATCH_LAS1_ALT)*exp(LC_PLUS_MATCH_LAS1_ALT.1)*exp((LC_PLUS_MATCH_LAS1_ALT.2)*c.adj/(c.adj+180))*X.ALT.TV,
		c.adj*exp(LC_PLUS_MATCH_LAS1_ALT)*exp(LC_PLUS_MATCH_LAS1_ALT.1)*exp((LC_PLUS_MATCH_LAS1_ALT.2)*c.adj/(c.adj+180))*c.adj/(c.adj+180)*X.ALT.TV,
		c.adj*exp(LC_PLUS_MATCH_LAS1_ALT)*exp(LC_PLUS_MATCH_LAS1_ALT.1)*exp((LC_PLUS_MATCH_LAS1_ALT.2)*c.adj/(c.adj+180))*X.ALT.CONST)*
			(1-POT_CENSOR_ALT_C_CENSOR)*(1-POT_CENSOR_ALT_C_IND2)+

     (cbind(POT_CENSOR*exp(LC_PLUS_MATCH_LAS)*exp(LC_PLUS_MATCH_LAS.1)*exp((LC_PLUS_MATCH_LAS.2)*POT_CENSOR/(POT_CENSOR+180))*X.MAX.FUTURE.TV,
		POT_CENSOR*exp(LC_PLUS_MATCH_LAS)*exp(LC_PLUS_MATCH_LAS.1)*exp((LC_PLUS_MATCH_LAS.2)*POT_CENSOR/(POT_CENSOR+180))*
     		POT_CENSOR/(POT_CENSOR+180)*X.MAX.FUTURE.TV,
		POT_CENSOR*exp(LC_PLUS_MATCH_LAS)*exp(LC_PLUS_MATCH_LAS.1)*exp((LC_PLUS_MATCH_LAS.2)*POT_CENSOR/(POT_CENSOR+180))*X.MAX.FUTURE.CONST))*
			POT_CENSOR_C_IND*POT_CENSOR_ALT_C_IND2)*(1-POT_CENSOR_ALT_C_CENSOR)+

	cbind((f.func.2)*X.AT.LTX.TV,(f.func.2)*RESID_AFTER_TRANS/(RESID_AFTER_TRANS+180)*X.AT.LTX.TV,(f.func.2)*X.AT.LTX.CONST)*POT_CENSOR_ALT_C_CENSOR)
gradient.c <- rbind(gradient.c.1,gradient.c.2)

gradient <-(gradient.x+gradient.c)

if(q <= 20) {
coef.beta <- -solve(gradient)%*%matrix(estimating.eq,ncol=1)*0.10+coef.beta
}

if(q > 20 & q <=30) {
coef.beta <- -solve(gradient)%*%matrix(estimating.eq,ncol=1)*0.25+coef.beta
}

if(q > 30) {
coef.beta <- -solve(gradient)%*%matrix(estimating.eq,ncol=1)+coef.beta
}


print(proc.time()-start)
rownames(coef.beta) <- names.beta.parm.tv

coef.total.SNNSAFTM <- cbind(coef.total.SNNSAFTM,coef.beta)
print(q)
print(c(estimating.eq) )
print(c(coef.beta) )
}


########################################################
### add in standard errors
########################################################

h.func.full <- cbind(mean.func.1*X.CURRENT.TV,mean.func.2*X.CURRENT)

part.1 <- t(h.func.full*PROB_TRANS)%*%h.func.full
h.func.prob <- h.func.full*PROB_TRANS
h.func.prob.sum <- NULL
for (i in 1:ncol(h.func.prob)) {
h.func.prob.sum <- cbind(h.func.prob.sum,tapply(h.func.prob[,i],T_TRR_ID_CODE,sum))
}
part.2 <- t(h.func.prob.sum)%*%h.func.prob.sum
B.matrix <- (part.1-part.2)
A.matrix <- gradient
v.cov <- solve(A.matrix)%*%B.matrix%*%t(solve(A.matrix))
se <- sqrt(diag(v.cov))


print("Convergence of Parameters in NSAFTM: Relative Change Between Iterations")
rel.change <- max(abs((coef.total.SNNSAFTM[,36]-coef.total.SNNSAFTM[,35])/coef.total.SNNSAFTM[,35]))
print(rel.change)


