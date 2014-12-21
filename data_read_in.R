##############################################
### read-in data - GROUP A
##############################################
print("Loading in Data for Group A")
data.LAS.A <- fread("data/combined_longitudinal_listing_prob_trans_GROUP_A_r2.txt")
data.LAS.A <- as.data.frame(data.LAS.A)
data.LAS.A <- data.LAS.A[is.na(data.LAS.A$PROB_TRANS)==F,]
n.subj.A <- nrow(data.LAS.A)

##############################################
### read-in data - GROUP B
##############################################
print("Loading in Data for Group B")
data.LAS.B <- fread("data/combined_longitudinal_listing_prob_trans_GROUP_B_r2.txt")
data.LAS.B <- as.data.frame(data.LAS.B)
data.LAS.B <- data.LAS.B[is.na(data.LAS.B$PROB_TRANS)==F,]
n.subj.B <- nrow(data.LAS.B)

##############################################
### read-in data - GROUP C
##############################################
print("Loading in Data for Group C")
data.LAS.C <- fread("data/combined_longitudinal_listing_prob_trans_GROUP_C_r2.txt")
data.LAS.C <- as.data.frame(data.LAS.C)
data.LAS.C <- data.LAS.C[is.na(data.LAS.C$PROB_TRANS)==F,]
n.subj.C <- nrow(data.LAS.C)


##############################################
### read-in data - GROUP D
##############################################
print("Loading in Data for Group D")
data.LAS.D <- fread("data/combined_longitudinal_listing_prob_trans_GROUP_D_r2.txt")
data.LAS.D <- as.data.frame(data.LAS.D)
data.LAS.D <- data.LAS.D[is.na(data.LAS.D$PROB_TRANS)==F,]
n.subj.D <- nrow(data.LAS.D)

##############################################
### read-in data - combine all native disease groupings
##############################################
print("Manipulating Data")
data.LAS <- rbind(data.LAS.A,data.LAS.B,data.LAS.C,data.LAS.D)
n.subj <- nrow(data.LAS)
rm(data.LAS.A,data.LAS.B,data.LAS.C,data.LAS.D)

###############################################
### get variables that we will need as separate vectors 
### from data frame and then remove data frame
################################################

##############
### follow-up and at risk time; failure indicator
##############
POT_CENSOR <- data.LAS$POT_CENSOR
EVENT_INDICATOR_USE <- data.LAS$EVENT_INDICATOR
RESID_LIFE <- data.LAS$RESID_LIFE
RESID_TO_TRANS <- data.LAS$RESID_TO_TRANS
RESID_AFTER_TRANS <- data.LAS$RESID_AFTER_TRANS

##############
### transplantation indicator and probability of transplantation
##############
TRANSPLANTED <- data.LAS$TRANSPLANTED
EVER_TRANSPLANTED <- data.LAS$EVER_TRANSPLANTED  
PROB_TRANS <- data.LAS$PROB_TRANS
T_TRR_ID_CODE <- data.LAS$T_TRR_ID_CODE
TPT <- TRANSPLANTED-PROB_TRANS

##############
### recipient characteristics
##############
WL_ID_CODE <- data.LAS$WL_ID_CODE
GROUPING <- data.LAS$GROUPING
GROUP_A <- ifelse(GROUPING=="A",1,0)
GROUP_B <- ifelse(GROUPING=="B",1,0)
GROUP_C <- ifelse(GROUPING=="C",1,0)
GROUP_D <- ifelse(GROUPING=="D",1,0)
MATCH_LAS <- data.LAS$MATCH_LAS 
END_MATCH_LAS <- data.LAS$T_MATCH_LAS1
END_MATCH_LAS <- ifelse(is.na(END_MATCH_LAS)==T,1,END_MATCH_LAS) 
INIT_AGE_USE <- data.LAS$INIT_AGE_USE
INIT_AGE_USE[GROUP_A==1] <- INIT_AGE_USE[GROUP_A==1]-mean(INIT_AGE_USE[GROUP_A==1])
INIT_AGE_USE[GROUP_B==1] <- INIT_AGE_USE[GROUP_B==1]-mean(INIT_AGE_USE[GROUP_B==1])
INIT_AGE_USE[GROUP_C==1] <- INIT_AGE_USE[GROUP_C==1]-mean(INIT_AGE_USE[GROUP_C==1])
INIT_AGE_USE[GROUP_D==1] <- INIT_AGE_USE[GROUP_D==1]-mean(INIT_AGE_USE[GROUP_D==1])
MOV_SUM_TRUE <- data.LAS$mov_sum
MOV_SUM_USE <- pmin(data.LAS$mov_sum,100)-50

##############
### donor, surgery, and matching characteristics
##############
T_AGE_DON_B <- data.LAS$T_AGE_DON_B
AGE_DON_B <- data.LAS$AGE_DON_B
AGE_DON_B <- ifelse(is.na(AGE_DON_B)==T,0,AGE_DON_B)

T_HIST_CIG_DON_Y <- data.LAS$T_HIST_CIG_DON_Y
HIST_CIG_DON_Y <- data.LAS$HIST_CIG_DON_Y
HIST_CIG_DON_Y <- ifelse(is.na(HIST_CIG_DON_Y)==T,0,HIST_CIG_DON_Y)

T_SINGLE_TRANS <- 1-data.LAS$T_BILATERAL
SINGLE_TRANS <- 1-data.LAS$BILATERAL
SINGLE_TRANS <- ifelse(is.na(SINGLE_TRANS)==T,0,SINGLE_TRANS)

T_HGT_DIFF <- data.LAS$T_HGT_DIFF 
T_HGT_DIFF <- pmax(T_HGT_DIFF,0)
HGT_DIFF <- data.LAS$HGT_CM_TCR-data.LAS$HGT_CM_DON_CALC1
HGT_DIFF <- ifelse(is.na(HGT_DIFF)==T,0,HGT_DIFF)
HGT_DIFF_TRUE <- HGT_DIFF
HGT_DIFF <- pmax(HGT_DIFF,0)

rm(data.LAS)

########################################################
### define names of parameters for tables and within R environment
########################################################

## Parameters for scale accelerated failure time model
names.beta.parm.table <- c("Intercept", "Single Lung", "LAS","LAS Squarred", "Group A (ref: Group D)",
	"Group B (ref: Group D)", "Group C (ref: Group D)", "Age Recipient (per 5 year increase)",
	"Age Donor (>55)", "Cigarette Donor", "2 year Center Volume (per 10 increase)",
	"Size Difference (per 5 cm increase)")
names.beta.parm <- c("Intercept", "Single_LTX", "LAS","LAS_SQ", "GROUP_A",
	"GROUP_B", "GROUP_C", "AGE", "AGE_DON", "CIG_DON", "CENTER_VOL", "HEIGHT_DIFF")
mult.beta.parm <- c(1,1,1,1,1,1,1,5,1,1,10,5)

## Parameters for time-varying accelerated failure time model
names.beta.parm.tv.table <- c("Intercept Bilateral Lung Early", "Intercept Bilateral Lung Late",
	"Intercept Sinlgle Lung Early", "Intercept Single Lung Late", "LAS","LAS Squarred", 
	"Group A (ref: Group D)", "Group B (ref: Group D)", "Group C (ref: Group D)", 
	"Age Recipient (per 5 year increase)", "Age Donor (>55)", "Cigarette Donor",
	"2 year Center Volume (per 10 increase)", "Size Difference (per 5 cm increase)")
names.beta.tv.parm <- c("INT_EARLY", "Single_LTX_EARLY", "INT_DIFF", "Single_LTX_DIFF", "LAS","LAS_SQ",
	"GROUP_A", "GROUP_B", "GROUP_C", "AGE", "AGE_DON", "CIG_DON", "CENTER_VOL", "HEIGHT_DIFF")
mult.beta.parm.tv <- c(1,1,1,1,1,1,1,1,1,5,1,1,10,5)

names.gamma.parm.table <- c("Intercept","LAS","LAS Squared","Group A","Group B","Group C",
	"Age","Center Volume")
names.gamma.parm <- c("Intercept","LAS","LAS_SQ","GROUP_A","GROUP_B","GROUP_C",
	"AGE","CENTER_VOL")

########################################################
### get design matrices used in the subsequent analysis
########################################################
MATCH_LAS_DELTA <- pmin(MATCH_LAS+delta.use,100)

## Design matrix at time of transplant for each subject
X.AT.LTX <- cbind(1,SINGLE_TRANS,(END_MATCH_LAS-30),(END_MATCH_LAS-30)^2,GROUP_A,GROUP_B,
	GROUP_C,INIT_AGE_USE,AGE_DON_B,HIST_CIG_DON_Y,MOV_SUM_USE,HGT_DIFF)
## Design matrix at current transplant (i.e. current value of patient characteristic and characteristics
## of the current organ transplant)
X.CURRENT <- cbind(1,T_SINGLE_TRANS,MATCH_LAS-30,(MATCH_LAS-30)^2,GROUP_A,GROUP_B,GROUP_C,
	INIT_AGE_USE,T_AGE_DON_B,T_HIST_CIG_DON_Y,MOV_SUM_USE,T_HGT_DIFF)
## Design matrix at current transplant for recipient characteristics
X.CURRENT.RECIP <- cbind(1,MATCH_LAS-30,(MATCH_LAS-30)^2,GROUP_A,GROUP_B,GROUP_C,INIT_AGE_USE,
	MOV_SUM_USE)
## Design matrix at maximum future value (used for artificial censoring)
X.MAX.FUTURE <- cbind(1,rep(0,n.subj),MATCH_LAS_DELTA-30,(MATCH_LAS_DELTA-30)^2,GROUP_A,GROUP_B,
	GROUP_C,INIT_AGE_USE,rep(0,n.subj),rep(0,n.subj),MOV_SUM_USE,rep(0,n.subj))
## Design matrix used in calculating gradients of estimating functions
X.ALT <- cbind(1,rep(0,n.subj),MATCH_LAS-30,(MATCH_LAS-30)^2,GROUP_A,GROUP_B,GROUP_C,
	INIT_AGE_USE, rep(0,n.subj),rep(0,n.subj),MOV_SUM_USE,rep(0,n.subj))

## Design matrices for the time-varying components
change.com <- 1:2
X.AT.LTX.TV <- X.AT.LTX[,change.comp]
X.CURRENT.TV <- X.CURRENT[,change.comp] 
X.MAX.FUTURE.TV <- X.MAX.FUTURE[,change.comp]
X.ALT.TV <- X.ALT[,change.comp] 

X.AT.LTX.CONST <- X.AT.LTX[,-change.comp]
X.CURRENT.CONST <- X.CURRENT[,-change.comp] 
X.MAX.FUTURE.CONST <- X.MAX.FUTURE[,-change.comp]
X.ALT.CONST <- X.ALT[,-change.comp] 

################################
### matrices to assess follow-up on transplant scale 
################################

LAS.seq <- seq(from=30,to=70,by=5)
LAS.seq.p <-LAS.seq+delta.use

X.matrix.followup <- cbind(1,0,(LAS.seq-30),(LAS.seq-30)^2,0,0,0,0,0,0,0,0,0,0)
X.matrix.followup.p <- cbind(1,0,(LAS.seq.p-30),(LAS.seq.p-30)^2,0,0,0,0,0,0,0,0,0,0)
