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

for (q in 1:35) {
start <- proc.time()
# values of coefficients of the structural model initially for those that are time-varying
coef.beta.early <-  matrix(c(coef.beta[1:length(change.comp),]), ncol=1)
# difference in values of coefficients of the structural model long-term - initial
# for those that are time-varying
coef.beta.early <-  matrix(c(coef.beta[1:length(change.comp),]), ncol=1)	
coef.beta.diff <-  matrix(c(coef.beta[-(1:length(change.comp)),])[change.comp], ncol=1)	
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

estimating.eq <- c(apply(mean.func.1*TPT*X.matrix.12,2,sum)[change.comp],apply(mean.func.2*TPT*X.matrix.12[,1:14],2,sum))

f.func.2 <- RESID_AFTER_TRANS*exp(LC_END_MATCH_LAS)*exp(LC_END_MATCH_LAS.1)*exp((LC_END_MATCH_LAS.2)*RESID_AFTER_TRANS/(RESID_AFTER_TRANS+180))

gradient.x.1 <- t((smooth.indicator.deriv.x(RESID_LIFE_SCALE,POT_CENSOR_ALT_C,eps)-
  exp(LC_MATCH_LAS)*smooth.min.deriv.x(RESID_LIFE_SCALE,POT_CENSOR_ALT_C,eps))*TPT*X.matrix.12[,1:2])%*%
		cbind((f.func.2)*X.matrix.11[,1:2],
				(f.func.2)*RESID_AFTER_TRANS/(RESID_AFTER_TRANS+180)*X.matrix.11[,1:2],
				(f.func.2)*X.matrix.11[,3:14])
gradient.x.2 <- t((smooth.indicator.deriv.x(RESID_LIFE_SCALE,POT_CENSOR_C,eps2)-
	exp(LC_MATCH_LAS)*smooth.min.deriv.x(RESID_LIFE_SCALE,POT_CENSOR_C,eps2))*TPT*X.matrix.12[,1:14])%*%
		cbind((f.func.2)*X.matrix.11[,1:2],
				(f.func.2)*RESID_AFTER_TRANS/(RESID_AFTER_TRANS+180)*X.matrix.11[,1:2],
				(f.func.2)*X.matrix.11[,3:14])
gradient.x <- rbind(gradient.x.1,gradient.x.2)

gradient.c.2 <- t((smooth.indicator.deriv.c(RESID_LIFE_SCALE,POT_CENSOR_C,eps2)-
  exp(LC_MATCH_LAS)*smooth.min.deriv.c(RESID_LIFE_SCALE,POT_CENSOR_C,eps2))*TPT*X.matrix.12[,1:14])%*%
	(cbind(POT_CENSOR*exp(LC_PLUS_MATCH_LAS)*exp(LC_PLUS_MATCH_LAS.1)*exp((LC_PLUS_MATCH_LAS.2)*POT_CENSOR/(POT_CENSOR+180))*X.matrix.11A[,1:2],
		POT_CENSOR*exp(LC_PLUS_MATCH_LAS)*exp(LC_PLUS_MATCH_LAS.1)*exp((LC_PLUS_MATCH_LAS.2)*POT_CENSOR/(POT_CENSOR+180))*POT_CENSOR/(POT_CENSOR+180)*X.matrix.11A[,1:2],
		POT_CENSOR*exp(LC_PLUS_MATCH_LAS)*exp(LC_PLUS_MATCH_LAS.1)*exp((LC_PLUS_MATCH_LAS.2)*POT_CENSOR/(POT_CENSOR+180))*X.matrix.11A[,3:14])*POT_CENSOR_C_IND*(1-POT_CENSOR_C_CENSOR)+	 
	cbind((f.func.2)*X.matrix.11[,1:2],(f.func.2)*RESID_AFTER_TRANS/(RESID_AFTER_TRANS+180)*X.matrix.11[,1:2],(f.func.2)*X.matrix.11[,3:14])*POT_CENSOR_C_CENSOR)

gradient.c.1 <- t((smooth.indicator.deriv.c(RESID_LIFE_SCALE,POT_CENSOR_ALT_C,eps)-
	exp(LC_MATCH_LAS)*smooth.min.deriv.c(RESID_LIFE_SCALE,POT_CENSOR_ALT_C,eps))*TPT*X.matrix.12[,1:2])%*%
	((cbind(c.adj*exp(LC_PLUS_MATCH_LAS1_ALT)*exp(LC_PLUS_MATCH_LAS1_ALT.1)*exp((LC_PLUS_MATCH_LAS1_ALT.2)*c.adj/(c.adj+180))*X.matrix.11ALT[,1:2],
		c.adj*exp(LC_PLUS_MATCH_LAS1_ALT)*exp(LC_PLUS_MATCH_LAS1_ALT.1)*exp((LC_PLUS_MATCH_LAS1_ALT.2)*c.adj/(c.adj+180))*c.adj/(c.adj+180)*X.matrix.11ALT[,1:2],
		c.adj*exp(LC_PLUS_MATCH_LAS1_ALT)*exp(LC_PLUS_MATCH_LAS1_ALT.1)*exp((LC_PLUS_MATCH_LAS1_ALT.2)*c.adj/(c.adj+180))*X.matrix.11ALT[,3:14])*(1-POT_CENSOR_ALT_C_CENSOR)*(1-POT_CENSOR_ALT_C_IND2)+

     (cbind(POT_CENSOR*exp(LC_PLUS_MATCH_LAS)*exp(LC_PLUS_MATCH_LAS.1)*exp((LC_PLUS_MATCH_LAS.2)*POT_CENSOR/(POT_CENSOR+180))*X.matrix.11A[,1:2],
		POT_CENSOR*exp(LC_PLUS_MATCH_LAS)*exp(LC_PLUS_MATCH_LAS.1)*exp((LC_PLUS_MATCH_LAS.2)*POT_CENSOR/(POT_CENSOR+180))*POT_CENSOR/(POT_CENSOR+180)*X.matrix.11A[,1:2],
		POT_CENSOR*exp(LC_PLUS_MATCH_LAS)*exp(LC_PLUS_MATCH_LAS.1)*exp((LC_PLUS_MATCH_LAS.2)*POT_CENSOR/(POT_CENSOR+180))*X.matrix.11A[,3:14]))*POT_CENSOR_C_IND*POT_CENSOR_ALT_C_IND2)*(1-POT_CENSOR_ALT_C_CENSOR)+

	cbind((f.func.2)*X.matrix.11[,1:2],(f.func.2)*RESID_AFTER_TRANS/(RESID_AFTER_TRANS+180)*X.matrix.11[,1:2],(f.func.2)*X.matrix.11[,3:14])*POT_CENSOR_ALT_C_CENSOR)
gradient.c <- rbind(gradient.c.1,gradient.c.2)

gradient <-(gradient.x+gradient.c)

if(q <= 20) {
coef.beta[-c(13,14),] <- -solve(gradient[c(1:12,15:16),c(1:12,15:16)])%*%matrix(estimating.eq[-c(13,14)],ncol=1)*0.10+coef.beta[-c(13,14),]
}

if(q > 20 & q <=30) {
coef.beta[-c(13,14),] <- -solve(gradient[c(1:12,15:16),c(1:12,15:16)])%*%matrix(estimating.eq[-c(13,14)],ncol=1)*0.25+coef.beta[-c(13,14),]
}

if(q > 30) {
coef.beta[-c(13,14),] <- -solve(gradient[c(1:12,15:16),c(1:12,15:16)])%*%matrix(estimating.eq[-c(13,14)],ncol=1)+coef.beta[-c(13,14),]
}


print(proc.time()-start)


coef.total.SNNSAFTM <- cbind(coef.total.SNNSAFTM,coef.beta)
print(q)
print(c(estimating.eq) )
print(c(coef.beta) )
}


########################################################
### add in standard errors
########################################################

h.func.full <- cbind(mean.func.1*X.matrix.12[,1:2],mean.func.2*X.matrix.12)

part.1 <- t(h.func.full*PROB_TRANS)%*%h.func.full
h.func.prob <- h.func.full*PROB_TRANS
h.func.prob.sum <- NULL
for (i in 1:ncol(h.func.prob)) {
h.func.prob.sum <- cbind(h.func.prob.sum,tapply(h.func.prob[,i],T_TRR_ID_CODE,sum))
}
part.2 <- t(h.func.prob.sum)%*%h.func.prob.sum
B.matrix <- (part.1-part.2)[c(1:12,15:16),c(1:12,15:16)]
A.matrix <- gradient[c(1:12,15:16),c(1:12,15:16)]
v.cov <- solve(A.matrix)%*%B.matrix%*%t(solve(A.matrix))
se <- sqrt(diag(v.cov))

#sqrt(v.cov[2,2]+v.cov[4,4]+2*v.cov[2,4])
#coef.beta[2,1]+coef.beta[4,1]
#(coef.beta[2,1]+coef.beta[4,1])/sqrt(v.cov[2,2]+v.cov[4,4]+2*v.cov[2,4])

###########################################################
### create table for model results
###########################################################

results <- cbind(coef.beta[-c(13,14)],se)
colnames(results) <- c("Estimate","StdError")
results[3,1] <- coef.beta[1,]+coef.beta[3,]
results[4,1] <- coef.beta[2,]+coef.beta[4,]

results[3,2] <- sqrt(matrix(c(1,1),nrow=1)%*%v.cov[c(1,3),c(1,3)]%*%matrix(c(1,1),ncol=1)) 
results[4,2] <- sqrt(matrix(c(1,1),nrow=1)%*%v.cov[c(2,4),c(2,4)]%*%matrix(c(1,1),ncol=1)) 
rownames(results) <-  c("Intercept Early","Sinlge Lung Early","Intercept Late", "Single Lung Late", 
"LAS","LAS Squarred", "Group A (ref: Group D)","Group B (ref: Group D)", "Group C (ref: Group D)",                      
"Age (per 5 increase)","Age Donor (>55)", "Cigarette Donor","2 year Center Volume (per 10 increase)","Size Difference (per 5 cm increase)")

results <- as.data.frame(results)
results$zscore <- results$Estimate/results$StdError


mult.beta.parm.nonscale <- c(1,1,1,1,1,1,1,1,1,5,1,1,10,5)
results[,1:2] <- results[,1:2]*mult.beta.parm.nonscale 
results <- cbind(results[,1:2],exp(results[,1]),results[,3])
colnames(results) <- c("Estimate","StdError","exp(Estimate)","Z-score")



print("Non-Scale AFTM Results")
print(results)
print(xtable(results,caption="Parameter estimates for accelerated failure time model. Note that we centered the LAS score at 30
and recipient age at the average age of patients in the same native disease group, 2-year center volume at 50, and other characteristics
as indicated in the chart below. Therefore, the interpretation of the intercept terms is the effect of transplantation at these reference levels.
",digits=c(1,5,5,3,2),align=c("l","|","c","c","c","c")),
include.rownames=TRUE,append=TRUE,caption.placement="top", sanitize.text.function = function(x){x})

print("Convergence of Parameters in NSAFTM: Relative Change Between Iterations")
rel.change <- max(abs((coef.total.SNNSAFTM[-c(13,14),36]-coef.total.SNNSAFTM[-c(13,14),35])/coef.total.SNNSAFTM[-c(13,14),35]))
print(rel.change)


