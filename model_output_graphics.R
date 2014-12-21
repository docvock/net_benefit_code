###########################################################
### create table for model results
###########################################################

results <- cbind(coef.beta,se)
colnames(results) <- c("Estimate","StdError")
med.LAS <- 38.86
mean.age.A <- 55.99137
mean.age.B <- 45.7975
mean.age.C <- 30.08858
mean.age.D <- 54.65367
med.age <- 58
med.age <- 45
med.center.vol <- 69

intercept.adj <- coef.beta[5]*(med.LAS-30)+coef.beta[6]*(med.LAS-30)^2+coef.beta[10]*(med.age-mean.age.D)+
		coef.beta[13]*(med.center.vol-50) 
results[,1] <- c(coef.beta[1]+intercept.adj,
coef.beta[1]+coef.beta[3]+intercept.adj,
coef.beta[2]+coef.beta[1]+intercept.adj,
coef.beta[2]+coef.beta[4]+coef.beta[1]+coef.beta[3]+intercept.adj,	
coef.beta[5]+2*coef.beta[6]*(med.LAS-30),
coef.beta[6],
coef.beta[7]+coef.beta[10]*(mean.age.D-mean.age.A),
coef.beta[8]+coef.beta[10]*(mean.age.D-mean.age.B),
coef.beta[9]+coef.beta[10]*(mean.age.D-mean.age.C),
coef.beta[10],
coef.beta[11],
coef.beta[12],
coef.beta[13],
coef.beta[14])
results[,1]

coef.beta.matrix <- matrix(0,nrow=14,ncol=14)
coef.beta.matrix[1,c(1,5,6,10,13)] <- c(1,(med.LAS-30),(med.LAS-30)^2,(med.age-mean.age.D),(med.center.vol-50))
coef.beta.matrix[2,c(1,3,5,6,10,13)] <- c(1,1,(med.LAS-30),(med.LAS-30)^2,(med.age-mean.age.D),(med.center.vol-50))
coef.beta.matrix[3,c(1,2,5,6,10,13)] <- c(1,1,(med.LAS-30),(med.LAS-30)^2,(med.age-mean.age.D),(med.center.vol-50))
coef.beta.matrix[4,c(1,2,3,4,5,6,10,13)] <- c(1,1,1,1,(med.LAS-30),(med.LAS-30)^2,(med.age-mean.age.D),(med.center.vol-50))
coef.beta.matrix[6,6] <- coef.beta.matrix[10,10] <- coef.beta.matrix[11,11] <- 
	coef.beta.matrix[12,12] <- coef.beta.matrix[13,13] <- coef.beta.matrix[14,14] <-1
coef.beta.matrix[5,c(5,6)] <- c(1,2*(med.LAS-30))
coef.beta.matrix[7,c(7,10)] <- c(1,(mean.age.D-mean.age.A))
coef.beta.matrix[8,c(8,10)] <- c(1,(mean.age.D-mean.age.B))
coef.beta.matrix[9,c(9,10)] <- c(1,(mean.age.D-mean.age.C))
c(coef.beta.matrix%*%coef.beta)
results[,1]



results[,2] <- sqrt(diag(coef.beta.matrix%*%v.cov%*%t(coef.beta.matrix))) 
rownames(results) <-  c("Intercept Bilateral Early","Intercept Bilateral Late", "Intercept Single Early", "Intercept Single Late",
"LAS","LAS Squarred", "Group A (ref: Group D)","Group B (ref: Group D)", "Group C (ref: Group D)",                      
"Age (per 5 increase)","Age Donor (>55)", "Donor Smoking History",
"2 year Center Volume (per 30 increase)","Size Difference (per 5 cm increase)")

results <- as.data.frame(results)
results$zscore <- results$Estimate/results$StdError

results[,1] <- results[,1]*mult.beta.parm.tv 
results[,2] <- results[,2]*abs(mult.beta.parm.tv) 
results <- cbind(results[,1:2],exp(results[,1]),results[,3])
colnames(results) <- c("Estimate","StdError","exp(Estimate)","Z-score")


print("Non-Scale AFTM Results")
print(results)
print(xtable(results,caption="Parameter estimates for structural nested accelerated failure time model. Note that we centered the continuous
covariates at the median value and categorical covariates at the most common value among transplant recipients given in Table 1. Therefore, 
the interpretation of the intercept terms is the effect of transplantation at these reference levels.
",digits=c(1,5,5,3,2),align=c("l","|","c","c","c","c")),
include.rownames=TRUE,append=TRUE,caption.placement="top", sanitize.text.function = function(x){x})

results.2 <- cbind(results[,1:3],exp(results[,1]-qnorm(0.975)*results[,2]),exp(results[,1]+qnorm(0.975)*results[,2]),2*pnorm(-abs(results[,4])))
colnames(results.2) <- c("Estimate","StdError","exp(Estimate)","Lower CI","Upper CI","p-value")

print(xtable(results.2,caption="Parameter estimates for structural nested accelerated failure time model. Note that we centered the continuous
covariates at the median value and categorical covariates at the most common value among transplant recipients given in Table 1. Therefore, 
the interpretation of the intercept terms is the effect of transplantation at these reference levels.
",digits=c(1,5,5,3,3,3,3),align=c("l","|","c","c","c","c","c","c")),
include.rownames=TRUE,append=TRUE,caption.placement="top", sanitize.text.function = function(x){x},
file="results/NSAFTM_results_r2.txt")


########################################################
### proportion of transplanted subjects that would benefit at 1, 2, and 3- years from transplantation 
########################################################

LC_END_MATCH_LAS.T <- LC_END_MATCH_LAS[TRANSPLANTED==1]
LC_END_MATCH_LAS.1.T <- LC_END_MATCH_LAS.1[TRANSPLANTED==1]
LC_END_MATCH_LAS.2.T <- LC_END_MATCH_LAS.2[TRANSPLANTED==1]
END_MATCH_LAS.T <- END_MATCH_LAS[TRANSPLANTED==1]
GROUPING.T <- GROUPING[TRANSPLANTED==1]
MOV_SUM_USE.T <- MOV_SUM_TRUE[TRANSPLANTED==1]
HIST_CIG_DON_Y.T <- HIST_CIG_DON_Y[TRANSPLANTED==1]
SINGLE_TRANS.T <- SINGLE_TRANS_TRUE[TRANSPLANTED==1]
HGT_DIFF.T <- HGT_DIFF_TRUE[TRANSPLANTED==1]
HGT_DIFF_ALT <- T_HGT_DIFF[TRANSPLANTED==1]
AGE_DON_B.T <- AGE_DON_B[TRANSPLANTED==1]

lower.LAS <- c(0,32.5,35,40,45,50)
upper.LAS <- c(32.5,35,40,45,50,100)
lower.CV <- c(0,25,50,75,100)
upper.CV <- c(25,50,75,100,300)
lower.HD <- c(-50,0,5,10)
upper.HD <- c(0,5,10,40)

GROUP.LEVEL <- c("A","B","C","D")
BINARY.LEVEL <- c(0,1)


lifetime_1yr_notrans <- 365*exp(LC_END_MATCH_LAS.T)*exp(LC_END_MATCH_LAS.1.T)*exp((LC_END_MATCH_LAS.2.T[1])*365/(365+180))
lifetime_2yr_notrans <- 730*exp(LC_END_MATCH_LAS.T)*exp(LC_END_MATCH_LAS.1.T)*exp((LC_END_MATCH_LAS.2.T)*730/(730+180))
lifetime_3yr_notrans <- 1095*exp(LC_END_MATCH_LAS.T)*exp(LC_END_MATCH_LAS.1.T)*exp((LC_END_MATCH_LAS.2.T)*1095/(1095+180))
lifetime_notrans <- cbind(lifetime_1yr_notrans,lifetime_2yr_notrans,lifetime_3yr_notrans)

print("Proportion of subjects with a 1-year survival benefit")
print(prop.table(table( lifetime_1yr_notrans<365)))
print("Proportion of subjects with a 2-year survival benefit")
print(prop.table(table( lifetime_2yr_notrans<730)))
print("Proportion of subjects with a 3-year survival benefit")
print(prop.table(table( lifetime_3yr_notrans<1095)))

########################
### Among all transplants the proportion that could be considered futile
### by patient and donor characteristics
########################
time_level <- c(365,730,1095)
prop.succ.total <- NULL
for(q in 1:3) {
lifetime_tyr_notrans <- lifetime_notrans[,q]
prop.succ <- prop.table(table( lifetime_tyr_notrans< time_level[q]))[2]
prop.succ <- c(prop.succ,NA)
for (j in 1:4) {
prop_table <- prop.table(table( lifetime_tyr_notrans[GROUPING.T==GROUP.LEVEL[j]] < time_level[q]))
prop.succ <- c(prop.succ,prop_table[2]) }
prop.succ <- c(prop.succ,NA)
for(j in 1:6) {
prop_table <- prop.table(table( lifetime_tyr_notrans[END_MATCH_LAS.T>=lower.LAS[j]&END_MATCH_LAS.T<upper.LAS[j]] <time_level[q]))
if(length(prop_table)>1) {
prop.succ <- c(prop.succ,prop_table[2]) 
} else{
prop.succ <- c(prop.succ,prop_table) } }
prop.succ <- c(prop.succ,NA)
for(j in 1:5) {
prop_table <- prop.table(table( lifetime_tyr_notrans[MOV_SUM_USE.T>=lower.CV[j]&MOV_SUM_USE.T<upper.CV[j]] < time_level[q]))
if(length(prop_table)>1) {
prop.succ <- c(prop.succ,prop_table[2]) 
} else{
prop.succ <- c(prop.succ,prop_table) } }
prop.succ <- c(prop.succ,NA)
for (j in 1:2) {
prop_table <- prop.table(table( lifetime_tyr_notrans[AGE_DON_B.T ==BINARY.LEVEL[j]] < time_level[q]))
prop.succ <- c(prop.succ,prop_table[2]) }
prop.succ <- c(prop.succ,NA)
for (j in 1:2) {
prop_table <- prop.table(table( lifetime_tyr_notrans[HIST_CIG_DON_Y.T==BINARY.LEVEL[j]] < time_level[q]))
prop.succ <- c(prop.succ,prop_table[2]) }
prop.succ <- c(prop.succ,NA)
for (j in 1:2) {
prop_table <- prop.table(table( lifetime_tyr_notrans[SINGLE_TRANS.T ==BINARY.LEVEL[j]] < time_level[q]))
prop.succ <- c(prop.succ,prop_table[2]) }
prop.succ <- c(prop.succ,NA)
for(j in 1:4) {
prop_table <- prop.table(table( lifetime_tyr_notrans[HGT_DIFF.T>=lower.HD[j]&HGT_DIFF.T<upper.HD[j]] < time_level[q]))
if(length(prop_table)>1) {
prop.succ <- c(prop.succ,prop_table[2]) 
} else{
prop.succ <- c(prop.succ,prop_table) } }
prop.succ.total <- cbind(prop.succ.total,prop.succ)
}


row.names.use <- c("Overall","Native Disease Grouping",
"\\hspace{1em} Obstructive","\\hspace{1em} Vascular","\\hspace{1em} Cystic/bronchiestasis","\\hspace{1em} Restrictive",
"LAS at Transplantation",paste("\\hspace{1em}",lower.LAS,"-",upper.LAS,sep=" "),
"2-year Center Volume",paste("\\hspace{1em}",lower.CV,"-",upper.CV,sep=" "),
"Age Donor", "\\hspace{1em} $<$55 years", "\\hspace{1em} $>$55 years",
"Donor cigarette history", "\\hspace{1em} $<$20 pack years","\\hspace{1em} $>$20 pack years and recent use",
"Type of Transplant","\\hspace{1em} Bilateral Transpalnt","\\hspace{1em} Single Transplant",
"Height Difference (recipient-donor)", "\\hspace{1em} $<$0 cm","\\hspace{1em} 0-5 cm", "\\hspace{1em} 5-10 cm","\\hspace{1em} $>$10 cm")

rownames(prop.succ.total ) <- row.names.use
colnames(prop.succ.total ) <- c("1-year","2-year","3-year")
print(prop.succ.total)

prop.succ.total <- prop.succ.total*100

print(xtable(prop.succ.total,caption="Percentage of transplants with 1-, 2-, and 3-year survival benefit by patient and donor characteristics.",
digits=c(1,1,1,1),align=c("l","|","c","c","c")),
include.rownames=TRUE,append=TRUE,caption.placement="top", sanitize.text.function = function(x){x},
file="results/NSAFTM_results_r2.txt"
)

######################################################################
######################################################################
### Plot relative survival benefit by LAS among those transplanted
######################################################################
######################################################################

rsb_1yr <- 1/(exp(LC_END_MATCH_LAS.T)*exp(LC_END_MATCH_LAS.1.T)*exp((LC_END_MATCH_LAS.2.T[1])*365/(365+180)))
rsb_2yr <- 1/(exp(LC_END_MATCH_LAS.T)*exp(LC_END_MATCH_LAS.1.T)*exp((LC_END_MATCH_LAS.2.T)*730/(730+180)))
rsb_3yr <- 1/(exp(LC_END_MATCH_LAS.T)*exp(LC_END_MATCH_LAS.1.T)*exp((LC_END_MATCH_LAS.2.T)*1095/(1095+180)))


n.trans <- length(END_MATCH_LAS.T)
rsb_data <- cbind(c(rsb_1yr, rsb_2yr, rsb_3yr), rep(END_MATCH_LAS.T, 3))
rsb_data <- as.data.frame(rsb_data)
rsb_data$Year <- c(rep("1-Year", n.trans), rep("2-Year", n.trans), rep("3-Year", n.trans))
colnames(rsb_data) <- c("SurvBen", "LAS_T", "Year")
pdf("graphics/trt_effect_LAS_rsb.pdf",width=10,height=6)
xyplot(log(SurvBen) ~ LAS_T |Year, data = rsb_data, xlab = "LAS at Transplantation", ylab = "log(Relative Survival Benefit)", layout = c(3,1))
dev.off()
#xyplot(log(SurvBen) ~ LAS_T |Year, data = rsb_data, panel = panel.hexbinplot)

cor(END_MATCH_LAS.T, log(rsb_1yr), method="spearman")
cor(END_MATCH_LAS.T, log(rsb_2yr), method="spearman")
cor(END_MATCH_LAS.T, log(rsb_3yr), method="spearman")

#######################################################################
#######################################################################
### Simple Criteria for assessing futile transplants
######################################################################
########################################################################

q <- 2
lifetime_tyr_notrans <- lifetime_notrans[,q]
#criteria <- ifelse(GROUPING.T=="A"&END_MATCH_LAS.T<35,1,0)
criteria <- ifelse(GROUPING.T=="A"&END_MATCH_LAS.T<35&MOV_SUM_USE.T<50,1,0)
#criteria <- ifelse(GROUPING.T=="A"&END_MATCH_LAS.T<35&(MOV_SUM_USE.T<50|
# (HGT_DIFF.T>5 & HIST_CIG_DON_Y.T==0)|
# (HGT_DIFF.T>5&AGE_DON_B.T==0)|
# (HGT_DIFF.T>5&MOV_SUM_USE.T<65)|
# (MOV_SUM_USE.T<65&HIST_CIG_DON_Y.T==0)|
# (MOV_SUM_USE.T<65&AGE_DON_B.T==0)|
# (HIST_CIG_DON_Y.T==0&AGE_DON_B.T==1)),1,0)

n <- 9091
n.fail <- table( lifetime_2yr_notrans<730)[1]
prop_table <- prop.table(table( lifetime_tyr_notrans[criteria==1] < time_level[q]))
print(prop_table)
table.qual <- table( lifetime_tyr_notrans[criteria==1] < time_level[q])
print(table.qual[1]/n)
print(table.qual[1]/n.fail)

########################################################
### initial and long-term ratio and equity point for variety of LAS-T
### and reference values for other covariates (including bilateral transplant)
########################################################
coef.beta.center <- c(coef.beta.matrix%*%coef.beta)
equity.func <- function(time,LAS) {
equity.func <- exp(coef.beta.center[1])*exp((-coef.beta.center[1]+coef.beta.center[2])*time*365/(time*365+180))*exp(coef.beta.center[5]*LAS+coef.beta.center[6]*LAS^2)-1
return(equity.func)
}

LAS.seq <- seq(from=30,to=70,by=5)-med.LAS
init.ratio <- exp(coef.beta.center[1])*exp(coef.beta.center[5]*LAS.seq+coef.beta.center[6]*LAS.seq^2)
long.term.ratio <- exp(coef.beta.center[1])*exp(-coef.beta.center[1]+coef.beta.center[2])*exp(coef.beta.center[5]*LAS.seq+coef.beta.center[6]*LAS.seq^2)
equity.point <- NULL
for (i in 1:length(LAS.seq)) {
if (init.ratio[i]>1) {
equity.point.i <- uniroot(equity.func,interval=c(0,7),LAS=LAS.seq[i])$root } else {
equity.point.i <- 0 }
equity.point <- c(equity.point,equity.point.i)
}

LAS.effect <- as.data.frame(cbind(LAS.seq+med.LAS,init.ratio,long.term.ratio,equity.point))
colnames(LAS.effect) <- c("LAS-T","Initial Ratio","Long-term Ratio","Equity Point(years)")
rownames(LAS.effect) <- NULL


print(xtable(LAS.effect,caption="The initial and long-term ratio of the amount of time a patient would have to remain on the waiting list to observe the same proportion of
patients dying post-transplantation as a function of LAS at transplantation. In addition, we provide the equity point, the time (years) when the survival proprotion post-transplant
is equal to the survival proportion among patients remaining on the waiting list. All other covariates were set at the reference level given in Table 1",
digits=c(0,0,2,2,2),align=c("l","l","|","c","c","c")),
include.rownames=FALSE,append=TRUE,caption.placement="top", sanitize.text.function = function(x){x},
file="results/NSAFTM_results_r2.txt"
)

####################################################
#### refined and group hypothesis tests
####################################################


###########################################################
### Test if the acceleration factor is different over native disease grouping
###########################################################

test.coef.matrix <- coef.beta.matrix[7:9,]
test.vcov.matrix <- test.coef.matrix%*%v.cov%*%t(test.coef.matrix)
test.effect <- test.coef.matrix%*%coef.beta
chi.sq <- t(test.effect)%*%solve(test.vcov.matrix)%*%test.effect
p.value <- 1-pchisq(chi.sq,nrow(test.coef.matrix))
test.result <- c(chi.sq,nrow(test.coef.matrix),p.value)
names(test.result) <- c("Chi.sq","df","p.value")
print("Chi-square test if acceleration factor varies between native disease")
print(test.result)

###########################################################
### Difference in acceleration factor between Group A and Group C 
###########################################################

test.coef.matrix <- matrix((coef.beta.matrix[7,]-coef.beta.matrix[9,]),nrow=1)
test.vcov.matrix <- test.coef.matrix%*%v.cov%*%matrix(test.coef.matrix,ncol=1)
test.effect <- matrix(test.coef.matrix,nrow=1)%*%coef.beta
chi.sq <- t(test.effect)%*%solve(test.vcov.matrix)%*%test.effect
p.value <- 1-pchisq(chi.sq,1)
test.result <- c(chi.sq,1,p.value)
print("Chi-square test if acceleration factor varies between native Group A and Group C")
names(test.result) <- c("Chi.sq","df","p.value")
print(test.result)

print("Effect of Group C relative to Group A and 95% CI")
print(c(exp(test.effect),exp(test.effect-qnorm(0.975)*sqrt(test.vcov.matrix)),exp(test.effect+qnorm(0.975)*sqrt(test.vcov.matrix))))


###########################################################
### Long term bilateral transplant
###########################################################

test.coef.matrix <- matrix((coef.beta.matrix[4,]-coef.beta.matrix[2,]),nrow=1)
test.vcov.matrix <- test.coef.matrix%*%v.cov%*%matrix(test.coef.matrix,ncol=1)
test.effect <- matrix(test.coef.matrix,nrow=1)%*%coef.beta
chi.sq <- t(test.effect)%*%solve(test.vcov.matrix)%*%test.effect
p.value <- 1-pchisq(chi.sq,1)
test.result <- c(chi.sq,1,p.value)
print("Chi-square test if acceleration factor differs early compared to Late for single lung transplant")
names(test.result) <- c("Chi.sq","df","p.value")
print(test.result)

print("Effect of Single to Bilateral and 95% CI")
print(c(exp(test.effect),exp(test.effect-qnorm(0.975)*sqrt(test.vcov.matrix)),exp(test.effect+qnorm(0.975)*sqrt(test.vcov.matrix))))


###########################################################
### Difference in exp(acceleration) factor for LAS-T of 30 versus 35
###########################################################

test.coef.vec  <- rep(0,14)
test.coef.vec[5] <- 5
test.coef.vec[6] <- 25
test.coef.matrix <- matrix(test.coef.vec,ncol=1)
test.vcov.matrix <- t(test.coef.matrix)%*%v.cov%*%test.coef.matrix
test.effect <- matrix(apply(test.coef.matrix*c(coef.beta),2,sum),ncol=1)
print("Effect of LAS-T 35 versus LAS-T 30 and 95% CI")
print(c(exp(test.effect),exp(test.effect-qnorm(0.975)*sqrt(test.vcov.matrix)),exp(test.effect+qnorm(0.975)*sqrt(test.vcov.matrix))))
print(c(1/exp(test.effect),1/exp(test.effect-qnorm(0.975)*sqrt(test.vcov.matrix)),1/exp(test.effect+qnorm(0.975)*sqrt(test.vcov.matrix))))
print(c(1-exp(test.effect),1-exp(test.effect-qnorm(0.975)*sqrt(test.vcov.matrix)),1-exp(test.effect+qnorm(0.975)*sqrt(test.vcov.matrix))))

###########################################################
### Difference in exp(acceleration) factor for LAS-T of 50 versus 55
###########################################################

test.coef.vec  <- rep(0,14)
test.coef.vec[5] <- 5
test.coef.vec[6] <- 25^2-20^2
test.coef.matrix <- matrix(test.coef.vec,ncol=1)
test.vcov.matrix <- t(test.coef.matrix)%*%v.cov%*%test.coef.matrix
test.effect <- matrix(apply(test.coef.matrix*c(coef.beta),2,sum),ncol=1)
print("Effect of LAS-T 55 versus LAS-T 50 and 95% CI")
print(c(exp(test.effect),exp(test.effect-qnorm(0.975)*sqrt(test.vcov.matrix)),exp(test.effect+qnorm(0.975)*sqrt(test.vcov.matrix))))
print(c(1/exp(test.effect),1/exp(test.effect-qnorm(0.975)*sqrt(test.vcov.matrix)),1/exp(test.effect+qnorm(0.975)*sqrt(test.vcov.matrix))))
print(c(1-exp(test.effect),1-exp(test.effect-qnorm(0.975)*sqrt(test.vcov.matrix)),1-exp(test.effect+qnorm(0.975)*sqrt(test.vcov.matrix))))


###########################################################
### graphics of the cummulative incidence by LAS strata
###########################################################
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

LAS.lower <- c(0,32.5,35,40,45,50)
LAS.upper <- c(32.5,35,40,45,50,100)
LAS.mean <- c(31.25,33.75,37.5,42.5,47.5,55)


pdf("graphics/trt_effect_LAS.pdf",width=10,height=6)
m <- matrix(c(1,2,3,4,5,6,7,7,7),nrow = 3,ncol = 3,byrow = TRUE)
layout(mat = m,heights = c(5,5,1))

for(j in 1:6) {
par(mar=c(c(5, 4, 4, 2) + 0.1))

T1A <- EVENT_INDICATOR_USE[TRANSPLANTED==1&MATCH_LAS>=LAS.lower[j]&MATCH_LAS<LAS.upper[j]]
T1B <- RESID_LIFE_SCALE[TRANSPLANTED==1&MATCH_LAS>=LAS.lower[j]&MATCH_LAS<LAS.upper[j]]

f1A <- survfit(Surv(T1B, T1A) ~ 1)

surv.1 <- f1A[[6]]
inc.1 <- 1-surv.1
time.1 <- f1A[[2]]
f1A[[6]] <- time.1
f1A[[2]] <- inc.1

T1T <- RESID_LIFE[TRANSPLANTED==1&MATCH_LAS>=LAS.lower[j]&MATCH_LAS<LAS.upper[j]]
T1E <- EVENT_INDICATOR_USE[TRANSPLANTED==1&MATCH_LAS>=LAS.lower[j]&MATCH_LAS<LAS.upper[j]]
f1T <- survfit(Surv(T1T, T1E) ~ 1)

#TIME_MAX2 <- f1T[[2]][which((1-f1T[[6]])<CI.max)]
#TIME_MAX2 <- TIME_MAX2[length(TIME_MAX2)]

surv.1T <- f1T[[6]]
inc.1T <- 1-surv.1T
time.1T <- f1T[[2]]
f1T[[6]] <- time.1T
f1T[[2]] <- inc.1T

plot(f1T,mark.time=F,conf.int=F,xlim=c(0,0.35),ylim=c(0,1095),ylab="Days",xlab="Cumulative Incidence",
	main=paste("LAS ",LAS.lower[j],"-",LAS.upper[j],sep=""),col="blue",lwd=3,cex.lab=1.5,cex.main=2)
lines(f1A,mark.time=F,conf.int=F,col="red",lwd=3,lty=2)
}

par(mar=c(c(1, 1, 1, 1) + 0.1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top",inset = 0,legend=c("Transplanted","Non-Transplanted"),col=c("blue","red"),
	lty=1:2,lwd=2,cex=1.4, horiz = TRUE)

dev.off()

###########################################################
### create table for length of follow-up on transplant scale
###########################################################

#####################
### Assume double lung transplant and other donor characteristics set as follows in the SNFTM 
### 	1) no smoking history
###		2) gender mismatch
###		3) death due to head trauma
###		4) recipient shorter than donor
### Assume that all other recipient characteristics remain constant over time. This is true for 
###	native disease group and age at registration but necessarily for center volume.
#####################

uni.function <- function(x,followup,covariates){
uni.function <- exp(covariates)*exp(coef.beta[1])*exp(coef.beta[3]*x*365/(x*365+180))*365*x -followup
return(uni.function)
}

LC_END_MATCH_LAS_FOLLOWUP <- c(X.matrix.followup[,3:12]%*%matrix(coef.beta[5:14],ncol=1))
LC_END_MATCH_LAS_FOLLOWUP.P1 <- c(X.matrix.followup.p[,3:12]%*%matrix(coef.beta[5:14],ncol=1))

pot.censor.p <- exp(LC_END_MATCH_LAS_FOLLOWUP.P1)*exp(coef.beta[1])*
			exp(coef.beta[3]*6.41*365/(6.41*365+180))*365*6.41
time.trans.full <- NULL
for (i in 1:length(LC_END_MATCH_LAS_FOLLOWUP)) {
time.trans.full <- c(time.trans.full,
	uniroot(uni.function,interval=c(0,6.41),followup=pot.censor.p[i],covariates=LC_END_MATCH_LAS_FOLLOWUP[i])$root)
}

pot.censor.p <- exp(LC_END_MATCH_LAS_FOLLOWUP.P1)*exp(coef.beta[1])*
			exp(coef.beta[3]*4*365/(4*365+180))*365*4
time.trans.4yr <- NULL
for (i in 1:length(LC_END_MATCH_LAS_FOLLOWUP)) {
time.trans.4yr <- c(time.trans.4yr,
	uniroot(uni.function,interval=c(0,4),followup=pot.censor.p[i],covariates=LC_END_MATCH_LAS_FOLLOWUP[i])$root)
}

pot.censor.p <- exp(LC_END_MATCH_LAS_FOLLOWUP.P1)*exp(coef.beta[1])*
			exp(coef.beta[3]*2*365/(2*365+180))*365*2
time.trans.2yr <- NULL
for (i in 1:length(LC_END_MATCH_LAS_FOLLOWUP)) {
time.trans.2yr <- c(time.trans.2yr,
	uniroot(uni.function,interval=c(0,2),followup=pot.censor.p[i],covariates=LC_END_MATCH_LAS_FOLLOWUP[i])$root)
}
follow_up <- cbind(LAS.seq,time.trans.full,time.trans.4yr,time.trans.2yr)
print("Length of follow-up on transplant scale with 6.41, 4, and 2 year follow-up on non-transformed follow-up")
colnames(follow_up) <- c("Current LAS","6.41 Years","4 Years","2 Years")
print(follow_up,digits=3)

########################################################
### proportion of subjects with ad hoc approach for censoring
########################################################
POT_CENSOR_C <- pmin(POT_CENSOR*exp(LC_PLUS_MATCH_LAS)*exp(LC_PLUS_MATCH_LAS.1)*exp((LC_PLUS_MATCH_LAS.2)*POT_CENSOR/(POT_CENSOR+180)),(POT_CENSOR))


print("Proportion of (i,j) observations with ad hoc appraoch for censoring")
print(prop.table(table((POT_CENSOR_C<=(RESID_LIFE_SCALE)+1.0),(EVENT_INDICATOR_USE==1)))[1,1])
ind.fail <- which((POT_CENSOR_C > (RESID_LIFE_SCALE)+1.0) & (EVENT_INDICATOR_USE==0))
print("Proportion of (i,j) subjects with at least one instance of ad hoc appraoch for censoring")
print(length(unique(WL_ID_CODE[ind.fail]))/length(unique(WL_ID_CODE)))

GROUP.LEVEL <- c("A","B","C","D")
for (j in 1:4) {
print(paste("Proportion of (i,j) observations with ad hoc appraoch for censoring in Group ",GROUP.LEVEL[j],sep=""))
print(prop.table(table((POT_CENSOR_C[GROUPING==GROUP.LEVEL[j]]<=(RESID_LIFE_SCALE[GROUPING==GROUP.LEVEL[j]])+1.0),(EVENT_INDICATOR_USE[GROUPING==GROUP.LEVEL[j]]==1)))[1,1])
ind.fail <- which((POT_CENSOR_C[GROUPING==GROUP.LEVEL[j]] > (RESID_LIFE_SCALE[GROUPING==GROUP.LEVEL[j]])+1.0) & (EVENT_INDICATOR_USE[GROUPING==GROUP.LEVEL[j]]==0))
print(paste("Proportion of (i,j) subjects with at least one instance of ad hoc appraoch for censoring in Group ",GROUP.LEVEL[j],sep=""))
print(length(unique((WL_ID_CODE[GROUPING==GROUP.LEVEL[j]])[ind.fail]))/length(unique(WL_ID_CODE[GROUPING==GROUP.LEVEL[j]])))
}
