#####################################################
#####################################################
### scale accelerated failure time model
### estimation using artifical censoring
#####################################################
#####################################################

coef.beta <- matrix(coef.beta.init,ncol=1)
coef.total.SNAFTM <- coef.beta

for (q in 1:25) {
start <- proc.time()

####

#################
### get needed linear combinations of data
#################
LC_END_MATCH_LAS <- c(X.AT.LTX%*%coef.beta)
LC_PLUS_MATCH_LAS <- c(X.MAX.FUTURE%*%coef.beta)

#################
## define transformed residual lifetime
#################
RESID_LIFE_SCALE <- RESID_TO_TRANS+(RESID_AFTER_TRANS)*exp(LC_END_MATCH_LAS)

#################
### define transformed artificial (potential) censoring time 
#################

	# minimum of transformed artificial (potential) censoring time and observed (potential) censoring time
POT_CENSOR_C <- pmin((POT_CENSOR)*exp(LC_PLUS_MATCH_LAS),(POT_CENSOR))
	# indicator for transformed (potential) censoring time as opposed to observed (potential) censoring time
POT_CENSOR_C_IND <- ifelse((POT_CENSOR)*exp(LC_PLUS_MATCH_LAS) < (POT_CENSOR),1,0)
	# if the subject is actually censored and the transformed artificial censoring is longer than transformed
	# residual follow-up, set the transformed artificial censoring equal to transformed residual follow-up
	# and set POT_CENSOR_C_CENSOR to be indicator if this is true 
POT_CENSOR_C_CENSOR <- ifelse((POT_CENSOR_C>(RESID_LIFE_SCALE)&(EVENT_INDICATOR_USE==0)),1,0)
POT_CENSOR_C <- ifelse((POT_CENSOR_C>(RESID_LIFE_SCALE)&(EVENT_INDICATOR_USE==0)),
	RESID_LIFE_SCALE,POT_CENSOR_C)

#################
### evaluate the estimating function
#################
mean.func <- smooth.indicator(RESID_LIFE_SCALE,POT_CENSOR_C,eps2)-exp(LC_MATCH_LAS)*
				smooth.min(RESID_LIFE_SCALE,POT_CENSOR_C,eps2)
estimating.eq <- apply(mean.func*TPT*X.CURRENT,2,sum)

#################
### obtain gradient of estimating function
#################
f.func <- exp(LC_END_MATCH_LAS)*RESID_AFTER_TRANS
gradient.x <- t((smooth.indicator.deriv.x(RESID_LIFE_SCALE,POT_CENSOR_C,eps2)-
	exp(LC_MATCH_LAS)*smooth.min.deriv.x(RESID_LIFE_SCALE,POT_CENSOR_C,eps2))*TPT*X.CURRENT)%*%
		cbind(f.func*X.AT.LTX)
gradient.c <- t((smooth.indicator.deriv.c(RESID_LIFE_SCALE,POT_CENSOR_C,eps2)-
  exp(LC_MATCH_LAS)*smooth.min.deriv.c(RESID_LIFE_SCALE,POT_CENSOR_C,eps2))*TPT*X.CURRENT)%*%
	(cbind(POT_CENSOR*exp(LC_PLUS_MATCH_LAS)*X.MAX.FUTURE)*POT_CENSOR_C_IND*(1-POT_CENSOR_C_CENSOR)+
	 cbind(f.func*X.AT.LTX)*POT_CENSOR_C_CENSOR)
gradient <- (gradient.x+gradient.c)

if(q <= 20) {
coef.beta <- -solve(gradient)%*%matrix(estimating.eq,ncol=1)*0.25+coef.beta}

if(q>20) {
coef.beta <- -solve(gradient)%*%matrix(estimating.eq,ncol=1)+coef.beta}

print(proc.time()-start)

coef.total.SNAFTM <- cbind(coef.total.SNAFTM,coef.beta)
print(q)
print(c(estimating.eq))
print(c(coef.beta))
}

########################################################
### add in standard errors
########################################################

h.func.full <- mean.func*X.CURRENT

part.1 <- t(h.func.full*PROB_TRANS)%*%h.func.full
h.func.prob <- h.func.full*PROB_TRANS
h.func.prob.sum <- NULL
for (i in 1:ncol(h.func.prob)) {
h.func.prob.sum <- cbind(h.func.prob.sum,tapply(h.func.prob[,i],T_TRR_ID_CODE,sum))
}
part.2 <- t(h.func.prob.sum)%*%h.func.prob.sum
B.matrix <- (part.1-part.2)
A.matrix <- (gradient)
v.cov <- solve(A.matrix)%*%B.matrix%*%t(solve(A.matrix))
se <- sqrt(diag(v.cov))


###########################################################
### create table for model results
###########################################################

results <- cbind(coef.beta,se)
colnames(results) <- c("Estimate","StdError")
results <- as.data.frame(results)
results$zscore <- results$Estimate/results$StdError
rownames(results) <- names.beta.parm.table
results[,1:2] <- results[,1:2]*mult.beta.parm 
results <- cbind(results[,1:2],exp(results[,1]),results[,3])
colnames(results) <- c("Estimate","StdError","exp(Estimate)","Z-score")

print("Scale AFTM Results")
print(results)

print(xtable(results,caption="Parameter estimates for scale accelerated failure time model. Note that we centered the LAS score at 30 and recipient age at the average age of patients in the same native disease group, 2-year center volume at 50, and other characteristics as indicated in the chart below. 
Therefore, the interpretation of the intercept terms is the effect of transplantation at these reference levels.",
	digits=c(1,5,5,3,2),align=c("l","|","c","c","c","c")),
	include.rownames=TRUE,append=TRUE,caption.placement="top", sanitize.text.function = function(x){x},
	file="results/SAFTM_results.txt")

