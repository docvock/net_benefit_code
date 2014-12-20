#######################################################
### R code to generate the KM cumulative incidence plots for WL dynamics
### and post transplant survival as well as numeric estimates
###
###
### PROGRAMS RELIED ON: combined_cohort_development.sas
### DATA FILES USED: combined_patient_record.txt
### ADDITIONAL PROGRAMS USED: none
### DATA FILES GENERATED: none
### PDF FILES GENERATED: WL_dynamics_overall.pdf and post_LTx_overall.pdf 
#######################################################


###################
### Functions and libraries used
###################
library(survival)
library(cmprsk)
library(date)

#######################
### read-in data
#######################
LAS.data <- read.table("data/combined_patient_record.txt",sep="\t",header=T,quote="\"")

#######################
### generate cumulative incidence curves for waitlist dynamics
### 	both overall and by LAS strata
#######################

WL_TIME <- ifelse(LAS.data$TRANSPLANTED==1,
as.Date(LAS.data$TX_DATE1,format="%m/%d/%Y")-as.Date(LAS.data$INIT_DATE_LAS1, format="%m/%d/%Y"),
as.Date(LAS.data$END_DATE_USE,origin="1960-01-01")-as.Date(LAS.data$INIT_DATE_LAS1, format="%m/%d/%Y"))
WL_TIME <- ifelse(WL_TIME<(-100),2,WL_TIME)
WL_FAIL <- ifelse(LAS.data$TRANSPLANTED==1,2,LAS.data$FAILURE)

TRANS_WL <- cuminc(ftime=WL_TIME, fstatus=WL_FAIL, cencode=0, na.action=na.omit)
TRANS_WL1 <- cuminc(ftime=WL_TIME/365, fstatus=WL_FAIL, cencode=0, na.action=na.omit)
print("Cumulative Incidence of Transplantation and Death")
print("1=Death,2=Transplantation")
print(TRANS_WL1,maxtime=2.5)

pdf("graphics/WL_dynamics_overall.pdf",width=8,height=6)
m <- matrix(c(1,2),nrow = 2,ncol = 1,byrow = TRUE)
layout(mat = m,heights = c(5,1))

plot(TRANS_WL[[1]]$time,TRANS_WL[[1]]$est,xlim=c(0,730),ylim=c(0,1),
	xlab="Days Since Initial Registration",ylab="Proportion",lty=1,
	col="blue",lwd=3,type="l",main="Waitlist Dynamics")
lines(TRANS_WL[[2]]$time,TRANS_WL[[2]]$est,col="red",lwd=3,lty=2)
surv.use <- NULL
for (i in 0:730) {
t1 <- TRANS_WL[[1]]$est[which(TRANS_WL[[1]]$time<=i)]
t1 <- t1[length(t1)]
t2 <- TRANS_WL[[2]]$est[which(TRANS_WL[[1]]$time<=i)]
t2 <- t2[length(t2)]
surv.use <- c(surv.use,1-t1-t2)
}
lines(0:730,surv.use,lwd=3) 
par(mar=c(c(1, 1, 1, 1) + 0.1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top",inset = 0,c("Died","Transplanted","Alive and Not Transplanted"),
	col=c("blue","red","black"),lty=c(1,2,1),lwd=3,cex=0.9, horiz = TRUE)
dev.off()

lower.LAS <- c(30,32.5,35,40,45,50)
upper.LAS <- c(32.5,35,40,45,50,60)

pdf("graphics/WL_dynamics_LAS_strata.pdf",width=8,height=6)
m <- matrix(c(1,2,3,4,5,6,7,7,7),nrow = 3,ncol = 3,byrow = TRUE)
layout(mat = m,heights = c(3,3,1))

for (j in 1:6) {
TRANS_WL <- cuminc(ftime=WL_TIME[LAS.data$INIT_MATCH_LAS<upper.LAS[j]&
		LAS.data$INIT_MATCH_LAS>=lower.LAS[j]], 
	fstatus=WL_FAIL[LAS.data$INIT_MATCH_LAS<upper.LAS[j]&
			LAS.data$INIT_MATCH_LAS>=lower.LAS[j]], cencode=0, na.action=na.omit)
plot(TRANS_WL[[1]]$time,TRANS_WL[[1]]$est,xlim=c(0,730),ylim=c(0,1),
	xlab="Days Since Initial Registration",ylab="Proportion",lty=1,
	col="blue",lwd=3,type="l",main=paste("Initial LAS ",lower.LAS[j],"-",upper.LAS[j],sep=""))
lines(TRANS_WL[[2]]$time,TRANS_WL[[2]]$est,col="red",lwd=3,lty=2)
surv.use <- NULL
for (i in 0:730) {
t1 <- TRANS_WL[[1]]$est[which(TRANS_WL[[1]]$time<=i)]
t1 <- t1[length(t1)]
t2 <- TRANS_WL[[2]]$est[which(TRANS_WL[[2]]$time<=i)]
t2 <- t2[length(t2)]
surv.use <- c(surv.use,1-t1-t2)
}
lines(0:730,surv.use,lwd=3) 
print(j)
}
par(mar=c(c(1, 1, 1, 1) + 0.1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top",inset = 0,c("Died","Transplanted","Alive and Not Transplanted"),
	col=c("blue","red","black"),lty=c(3,2,1),lwd=3,cex=1.15, horiz = TRUE)
dev.off()

#######################
### generate survival curves post-transplantation
#######################

TX_TIME <- as.numeric(as.Date(LAS.data$EVENT_DATE_USE[LAS.data$TRANSPLANTED==1],
	origin="1960-01-01")-as.Date(LAS.data$TX_DATE1[LAS.data$TRANSPLANTED==1],format="%m/%d/%Y"))
TX_FAIL <- LAS.data$EVENT_INDICATOR[LAS.data$TRANSPLANTED==1]

POST_TX <- survfit(Surv(TX_TIME, TX_FAIL) ~ 1)
print("Post-transplant graft survival")
print(summary(POST_TX,times=c(365/2,365,730,1095)))

pdf("graphics/post_LTx_overall.pdf",width=6,height=6)
plot(POST_TX,mark.time=F,conf.int=F,xlab="Days Since Trnasplantation",
	ylab="Survival",xlim=c(0,1095),main="Post-Transplantation Graft Survival")
dev.off()





