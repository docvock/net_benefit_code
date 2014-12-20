
### R code to generate the tables for characteristics of recipeints (at registration and 
###  transplantation), donors, surgeries, and types of transplantation
###
### PROGRAMS RELIED ON: combined_cohort_development.sas
### DATA FILES USED: combined_patient_record.txt
### ADDITIONAL PROGRAMS USED: none
### DATA FILES GENERATED: none
### TXT FILES GENERATED: characteristics_table.txt (file with Latex code for tables)


###################
### Functions and libraries used
###################
library(xtable)
library(date)

cont.summary <- function(variable) {
cont.summary <- summary(variable)
if(length(cont.summary)==6) {
cont.summary <- c(cont.summary,0) }
cont.summary <- c(cont.summary[3],paste("(",cont.summary[2],", ",cont.summary[5],")",sep=""),
	round(cont.summary[7]/length(variable)*100,1))
return(cont.summary) }

multi.summary <- function(variable,describe) {
table.miss <- prop.table(table(is.na(variable)))*100
if(length(table.miss)==1) {
table.miss <- c(table.miss,0) } 
table.multi <- table(variable[is.na(variable)==F])
prop.multi <- round(prop.table(table.multi)*100,1)
multi.summary <- rbind(c(NA,NA,round(table.miss[2],1)),cbind(table.multi,prop.multi,NA))
rownames(multi.summary) <- describe
return(multi.summary) }

#######################
### read-in data
#######################
LAS.data <- read.table("data/combined_patient_record.txt",sep="\t",header=T,quote="\"")

#######################
### Define categorical variables - listing
#######################

### Ethnic Category
LAS.data$ETHCAT_G <- ifelse(LAS.data$ETHCAT==6,10,LAS.data$ETHCAT)
LAS.data$ETHCAT_G <- ifelse(LAS.data$ETHCAT==7,10,LAS.data$ETHCAT_G)
LAS.data$ETHCAT_G <- ifelse(LAS.data$ETHCAT==9,10,LAS.data$ETHCAT_G)
LAS.data$ETHCAT_G <- factor(LAS.data$ETHCAT_G)
r.ETHCAT_G <- c("\\hspace{1em} White","\\hspace{1em} Black","\\hspace{1em} Hispanic",
"\\hspace{1em} Asian","\\hspace{1em} Other")
### Vent Use
r.VENT_USE <- c("\\hspace{1em} BiPAP","\\hspace{1em} Continuous mechanical",
	"\\hspace{1em} Intermittent mechanical","\\hspace{1em} No assisted ventilation",
	"\\hspace{1em} CPAP")
### Diabetes
LAS.data$DIABETES1_G <- ifelse(LAS.data$DIABETES1=="I",1,
				 ifelse(LAS.data$DIABETES1=="N",2,
				 ifelse(LAS.data$DIABETES1=="U",NA,
				 ifelse(LAS.data$DIABETES1=="X",3,
				 ifelse(LAS.data$DIABETES1=="",NA,5)))))
r.DIABETES <- c("\\hspace{1em} Insulin dependent","\\hspace{1em} Non-insulin dependent or 
	dependence unknown","\\hspace{1em} Non-diabetic")
### NYHA class and functional class
LAS.data$FUNC_CLASS_USE <- ifelse(is.na(LAS.data$CALC_FUNC_STAT1)==F, LAS.data$CALC_FUNC_STAT1,
				ifelse(LAS.data$CALC_NYHA1<=7,1,
				ifelse(LAS.data$CALC_NYHA1==8,2,
				ifelse(LAS.data$CALC_NYHA1==9,3,4))))
r.FUNC_CLASS <- c("\\hspace{1em} No assistance","\\hspace{1em} Some assistance",
	"\\hspace{1em} Total assistance")


#######################
### Define categorical variables - at transplantation
#######################

LAS.data.trans <- LAS.data[LAS.data$TRANSPLANTED==1,]
### Ethnic Category Donor
LAS.data.trans$ETHCAT_DON_G <- ifelse(LAS.data.trans$ETHCAT_DON==6,10,LAS.data.trans$ETHCAT_DON)
LAS.data.trans$ETHCAT_DON_G <- ifelse(LAS.data.trans$ETHCAT_DON==7,10,LAS.data.trans$ETHCAT_DON_G)
LAS.data.trans$ETHCAT_DON_G <- ifelse(LAS.data.trans$ETHCAT_DON==9,10,LAS.data.trans$ETHCAT_DON_G)
LAS.data.trans$ETHCAT_DON_G <- factor(LAS.data.trans$ETHCAT_DON_G)

### Donor smoking history
LAS.data.trans$CIG <- ifelse(LAS.data.trans$HIST_CIG_DON1=="U",NA,
				ifelse(LAS.data.trans$HIST_CIG_DON1=="",NA,LAS.data.trans$HIST_CIG_DON1))
### CMV donor
LAS.data.trans$CMV <- ifelse(LAS.data.trans$CMV_DON1=="ND",NA,
				ifelse(LAS.data.trans$CMV_DON1=="U",NA,
				ifelse(LAS.data.trans$CMV_DON1=="",NA,LAS.data.trans$CMV_DON1)))
r.CMV <- c("\\hspace{1em} Indeterminant","\\hspace{1em} Negative","\\hspace{1em} Positive")

### COD donor
#LAS.data.trans$COD_CAD_DON1 <- ifelse(is.na(LAS.data.trans$COD_CAD_DON1)==T,999,LAS.data.trans$COD_CAD_DON1)
r.COD <- c("\\hspace{1em} Anoxia","\\hspace{1em} Cerebrovascular/Stroke",
	"\\hspace{1em} Head Trauma","\\hspace{1em} CNS Tumor","\\hspace{1em} Other COD")
### Diabetes
LAS.data.trans$T_DIABETES_G <- ifelse(LAS.data.trans$T_DIABETES=="I",1,
				 ifelse(LAS.data.trans$T_DIABETES=="N",2,
				 ifelse(LAS.data.trans$T_DIABETES=="U",NA,
				 ifelse(LAS.data.trans$T_DIABETES=="X",3,
				 ifelse(LAS.data.trans$T_DIABETES=="",NA,5)))))
### NYHA class and functional class
LAS.data.trans$T_FUNC_CLASS_USE <- ifelse(is.na(LAS.data.trans$T_CALC_FUNC_STAT)==F, 
	LAS.data.trans$T_CALC_FUNC_STAT,
				ifelse(LAS.data.trans$T_CALC_NYHA<=7,1,
				ifelse(LAS.data.trans$T_CALC_NYHA==8,2,
				ifelse(LAS.data.trans$T_CALC_NYHA==9,3,4))))


######################################
### Characteristics recipient at listing and transplantation
######################################

results.recipient.list <- rbind(cont.summary(LAS.data$INIT_AGE_LAS),
multi.summary(LAS.data$GENDER,c("Gender","\\hspace{1em}  Female","\\hspace{1em}  Male")),
multi.summary(LAS.data$ETHCAT_G,c("Race and ethnicity",r.ETHCAT_G)),
multi.summary(LAS.data$GROUPING,c("Native disease grouping","\\hspace{1em} Obstructive",
	"\\hspace{1em} Vascular","\\hspace{1em} Cystic/bronchiestasis","\\hspace{1em} Restrictive")),
cont.summary(LAS.data$INIT_MATCH_LAS1),
cont.summary(LAS.data$CALC_SIX_MIN_WALK1),
cont.summary(LAS.data$CALC_FVC_PRE1),
cont.summary(LAS.data$CALC_O21),
multi.summary(LAS.data$FUNC_CLASS_USE,c("Functional class",r.FUNC_CLASS)),
multi.summary(LAS.data$CALC_VENT_USE1,c("Ventilator use",r.VENT_USE)),
multi.summary(LAS.data$DIABETES1_G,c("Diabetes",r.DIABETES)),
cont.summary(LAS.data$CALC_BMI1),
cont.summary(LAS.data$CENTER_VOL1)
)
colnames(results.recipient.list) <-c("Median or n","(25$^{th}$, 75$^{th}$ centile) or \\%",
	"\\% Missing")
rownames(results.recipient.list)[c(1,16:19,34:35)]  <- c("Age","LAS",
	"Six min. walk distance (ft.)","FVC (% predicted)","O$_2$ required (L/day)","BMI (kg/m$^2$)",
	"Center Volume Preceding 2 years")

results.recipient.trans <- rbind(cont.summary(LAS.data.trans$AGE1),
multi.summary(LAS.data.trans$GENDER,c("Gender","\\hspace{1em}  Female","\\hspace{1em}  Male")),
multi.summary(LAS.data.trans$ETHCAT_G,c("Race and ethnicity",r.ETHCAT_G)),
multi.summary(LAS.data.trans$GROUPING,c("Native disease grouping","\\hspace{1em} Obstructive",
	"\\hspace{1em} Vascular","\\hspace{1em} Cystic/bronchiestasis","\\hspace{1em} Restrictive")),
cont.summary(LAS.data.trans$T_MATCH_LAS),
cont.summary(LAS.data.trans$T_CALC_SIX_MIN_WALK),
cont.summary(LAS.data.trans$T_CALC_FVC_PRE),
cont.summary(LAS.data.trans$T_CALC_O2),
multi.summary(LAS.data.trans$T_FUNC_CLASS_USE,c("Functional class",r.FUNC_CLASS)),
multi.summary(LAS.data.trans$T_CALC_VENT_USE,c("Ventilator use",r.VENT_USE)),
multi.summary(LAS.data.trans$T_DIABETES_G,c("Diabetes",r.DIABETES)),
cont.summary(LAS.data.trans$T_CALC_BMI),
cont.summary(LAS.data.trans$T_CENTER_VOL)
)

colnames(results.recipient.trans) <-c("Median or n","(25$^{th}$, 75$^{th}$ centile) or \\%","\\% Missing")
rownames(results.recipient.trans)[c(1,16:19,34:35)]  <- c("Age (yrs.)","LAS",
	"Six min. walk distance (ft.)","FVC (% predicted)","O$_2$ required (L/day)","BMI (kg/m$^2$)",
	"Center Volume Preceding 2 years")

recipient.characteristics <- cbind(results.recipient.list,results.recipient.trans)

print(xtable(recipient.characteristics,label=c("t:recipient"), caption=c("Characteristics of 13,040
unique patients included in the study cohort at time of first registration (or first non-zero LAS 
for patients listed prior to LAS era) and at transplantation for 9,091 patients transplanted during study period. Missing includes those observations categorized as unknown or test not performed."),
	align=c("p{1.75in}",rep(c("p{0.5in}","p{1in}","p{0.5in}"),2))),
	sanitize.rownames.function = identity,
	sanitize.colnames.function = identity,caption.placement="top",
	file="results/listing_characteristics_table.txt",append=F)	

######################################
### Donor characteristics, surgery characteristics, mathcing characteristics and transplant type 
######################################
LAS.data.trans$HGT_DIFF <- LAS.data.trans$END_HGT_CM-LAS.data.trans$HGT_CM_DON_CALC
LAS.data.trans$GENDER_MISS <- ifelse(LAS.data.trans$GENDER_DON1=="F" & LAS.data.trans$GENDER=="M",1,
						ifelse(LAS.data.trans$GENDER_DON1=="M" & LAS.data.trans$GENDER=="F",1,0))
LAS.data.trans$TX_TYPE <- ifelse(LAS.data.trans$TX_PROCEDUR_TY1==601,1,
					ifelse(LAS.data.trans$TX_PROCEDUR_TY1==602,1,0))
		

results.cont.donor <- rbind(cont.summary(LAS.data.trans$AGE_DON1),
multi.summary(LAS.data.trans$GENDER_DON1,c("Donor gender","\\hspace{1em} Unknown",
	"\\hspace{1em} Female","\\hspace{1em} Male")),
multi.summary(LAS.data.trans$ETHCAT_DON_G,c("Race and ethnicity",r.ETHCAT_G)),
multi.summary(LAS.data.trans$COD_CAD_DON1,c("Donor cause of death",r.COD)),
multi.summary(LAS.data.trans$T_DIABETES_G,c("Donor diabetes",r.DIABETES)),
multi.summary(LAS.data.trans$ECD_DONOR1,c("Expanded Criterion Donor","\\hspace{1em} No",
	"\\hspace{1em} Yes")),
multi.summary(LAS.data.trans$CMV,c("Donor CMV Histology",r.CMV)),
multi.summary(LAS.data.trans$CIG,c("Donor cigarette history", 
	"\\hspace{1em} $>$20 pack years and recent use","\\hspace{1em} $<$20 pack years")),
cont.summary(LAS.data.trans$PO21),
cont.summary(LAS.data.trans$BMI_DON_CALC1),
cont.summary(LAS.data.trans$ISCHTIME1),
cont.summary(LAS.data.trans$DISTANCE),
cont.summary(LAS.data.trans$HGT_DIFF),
multi.summary(LAS.data.trans$GENDER_MISS ,c("Gender mismatch","\\hspace{1em} No mismatch",
	"\\hspace{1em} Mismatch")),
multi.summary(LAS.data.trans$TX_TYPE ,c("Transplant type","\\hspace{1em} Bilateral",
	"\\hspace{1em} Single")))
colnames(results.cont.donor) <-c("Median or n","(25$^{th}$, 75$^{th}$ centile) or \\%","\\% Missing")

rownames(results.cont.donor)[c(1,32:36)] <- c("Donor age (yrs.)","Donor PO2","BMI Donor",
	"Ischemic time (hrs.)","Distance for graft retrieval (miles)",
	"Difference in Height (recipient-donor,cm)")

print(xtable(results.cont.donor,label=c("t:donor"),caption=c("Characteristics of donor, surgery, transplantation type, and organ match for the 9,091 transplantations included in the study 
	cohort."), align=c("p{1.75in}", rep(c("p{0.5in}","p{1in}","p{0.5in}"),1))),
	sanitize.rownames.function = identity,
	sanitize.colnames.function = identity,caption.placement="top",
	file="results/transplant_characteristics_table.txt",append=F)


##############################

