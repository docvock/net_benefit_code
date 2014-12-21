/*
#########################################################
### function to generate the estimated probabilities of receiving 
### an organ transplant
### this will be included in the main function
##########################################################
*/

proc sort data=total_final1; by T_TRR_ID_CODE;
run;

/* slim down data set for phreg use */
data phreg_data;
set total_final1 (keep=dummy_time transplanted match_las T_TRR_ID_CODE WL_ID_CODE region  T_TX_PROCEDUR_TY  T_AGE_DON
 HGT_DIFF INIT_DATE mov_sum INIT_AGE_USE GENDER T_GENDER_DON GROUPING T_HIST_CIG_DON_Y LIST_DAYS T_CMV_DON_P T_DIABETES_DON_Y T_COD_CAD_DON_HT);
if mov_sum=. then mov_sum=0;
BILATERAL = 0;
IF T_TX_PROCEDUR_TY=603 then BILATERAL=1;
IF T_TX_PROCEDUR_TY=604 then BILATERAL=1;

MATCH_LAS1 = max((MATCH_LAS - &match_las_k1)**3/(&match_las_k4-&match_las_k1)**2, 0)+
	(&match_las_k3-&match_las_k1)*max((MATCH_LAS - &match_las_k4)**3/(&match_las_k4-&match_las_k1)**2, 0)/(&match_las_k4-&match_las_k3)-
	(&match_las_k4-&match_las_k1)*max((MATCH_LAS - &match_las_k3)**3/(&match_las_k4-&match_las_k1)**2, 0)/(&match_las_k4-&match_las_k3);
MATCH_LAS2 = max((MATCH_LAS - &match_las_k2)**3/(&match_las_k4-&match_las_k1)**2, 0)+
	(&match_las_k3-&match_las_k2)*max((MATCH_LAS - &match_las_k4)**3/(&match_las_k4-&match_las_k1)**2, 0)/(&match_las_k4-&match_las_k3)-
	(&match_las_k4-&match_las_k2)*max((MATCH_LAS - &match_las_k3)**3/(&match_las_k4-&match_las_k1)**2, 0)/(&match_las_k4-&match_las_k3);
MOV_SUM1 = max((MOV_SUM - &mov_sum_k1)**3/(&mov_sum_k4-&mov_sum_k1)**2, 0)+
	(&mov_sum_k3-&mov_sum_k1)*max((MOV_SUM - &mov_sum_k4)**3/(&mov_sum_k4-&mov_sum_k1)**2, 0)/(&mov_sum_k4-&mov_sum_k3)-
	(&mov_sum_k4-&mov_sum_k1)*max((MOV_SUM - &mov_sum_k3)**3/(&mov_sum_k4-&mov_sum_k1)**2, 0)/(&mov_sum_k4-&mov_sum_k3);
MOV_SUM2 = max((MOV_SUM - &mov_sum_k2)**3/(&mov_sum_k4-&mov_sum_k1)**2, 0)+
	(&mov_sum_k3-&mov_sum_k2)*max((MOV_SUM - &mov_sum_k4)**3/(&mov_sum_k4-&mov_sum_k1)**2, 0)/(&mov_sum_k4-&mov_sum_k3)-
	(&mov_sum_k4-&mov_sum_k2)*max((MOV_SUM - &mov_sum_k3)**3/(&mov_sum_k4-&mov_sum_k1)**2, 0)/(&mov_sum_k4-&mov_sum_k3);
LIST_DAYS1 = max((LIST_DAYS - &LIST_DAYS_k1)**3/(&LIST_DAYS_k4-&LIST_DAYS_k1)**2, 0)+
	(&LIST_DAYS_k3-&LIST_DAYS_k1)*max((LIST_DAYS - &LIST_DAYS_k4)**3/(&LIST_DAYS_k4-&LIST_DAYS_k1)**2, 0)/(&LIST_DAYS_k4-&LIST_DAYS_k3)-
	(&LIST_DAYS_k4-&LIST_DAYS_k1)*max((LIST_DAYS - &LIST_DAYS_k3)**3/(&LIST_DAYS_k4-&LIST_DAYS_k1)**2, 0)/(&LIST_DAYS_k4-&LIST_DAYS_k3);
LIST_DAYS2 = max((LIST_DAYS - &LIST_DAYS_k2)**3/(&LIST_DAYS_k4-&LIST_DAYS_k1)**2, 0)+
	(&LIST_DAYS_k3-&LIST_DAYS_k2)*max((LIST_DAYS - &LIST_DAYS_k4)**3/(&LIST_DAYS_k4-&LIST_DAYS_k1)**2, 0)/(&LIST_DAYS_k4-&LIST_DAYS_k3)-
	(&LIST_DAYS_k4-&LIST_DAYS_k2)*max((LIST_DAYS - &LIST_DAYS_k3)**3/(&LIST_DAYS_k4-&LIST_DAYS_k1)**2, 0)/(&LIST_DAYS_k4-&LIST_DAYS_k3);
HGT_DIFF1 = max((HGT_DIFF - &HGT_DIFF_k1)**3/(&HGT_DIFF_k4-&HGT_DIFF_k1)**2, 0)+
	(&HGT_DIFF_k3-&HGT_DIFF_k1)*max((HGT_DIFF - &HGT_DIFF_k4)**3/(&HGT_DIFF_k4-&HGT_DIFF_k1)**2, 0)/(&HGT_DIFF_k4-&HGT_DIFF_k3)-
	(&HGT_DIFF_k4-&HGT_DIFF_k1)*max((HGT_DIFF - &HGT_DIFF_k3)**3/(&HGT_DIFF_k4-&HGT_DIFF_k1)**2, 0)/(&HGT_DIFF_k4-&HGT_DIFF_k3);
HGT_DIFF2 = max((HGT_DIFF - &HGT_DIFF_k2)**3/(&HGT_DIFF_k4-&HGT_DIFF_k1)**2, 0)+
	(&HGT_DIFF_k3-&HGT_DIFF_k2)*max((HGT_DIFF - &HGT_DIFF_k4)**3/(&HGT_DIFF_k4-&HGT_DIFF_k1)**2, 0)/(&HGT_DIFF_k4-&HGT_DIFF_k3)-
	(&HGT_DIFF_k4-&HGT_DIFF_k2)*max((HGT_DIFF - &HGT_DIFF_k3)**3/(&HGT_DIFF_k4-&HGT_DIFF_k1)**2, 0)/(&HGT_DIFF_k4-&HGT_DIFF_k3);
INIT_AGE_USE1 = max((INIT_AGE_USE - &INIT_AGE_USE_k1)**3/(&INIT_AGE_USE_k4-&INIT_AGE_USE_k1)**2, 0)+
	(&INIT_AGE_USE_k3-&INIT_AGE_USE_k1)*max((INIT_AGE_USE - &INIT_AGE_USE_k4)**3/(&INIT_AGE_USE_k4-&INIT_AGE_USE_k1)**2, 0)/(&INIT_AGE_USE_k4-&INIT_AGE_USE_k3)-
	(&INIT_AGE_USE_k4-&INIT_AGE_USE_k1)*max((INIT_AGE_USE - &INIT_AGE_USE_k3)**3/(&INIT_AGE_USE_k4-&INIT_AGE_USE_k1)**2, 0)/(&INIT_AGE_USE_k4-&INIT_AGE_USE_k3);
INIT_AGE_USE2 = max((INIT_AGE_USE - &INIT_AGE_USE_k2)**3/(&INIT_AGE_USE_k4-&INIT_AGE_USE_k1)**2, 0)+
	(&INIT_AGE_USE_k3-&INIT_AGE_USE_k2)*max((INIT_AGE_USE - &INIT_AGE_USE_k4)**3/(&INIT_AGE_USE_k4-&INIT_AGE_USE_k1)**2, 0)/(&INIT_AGE_USE_k4-&INIT_AGE_USE_k3)-
	(&INIT_AGE_USE_k4-&INIT_AGE_USE_k2)*max((INIT_AGE_USE - &INIT_AGE_USE_k3)**3/(&INIT_AGE_USE_k4-&INIT_AGE_USE_k1)**2, 0)/(&INIT_AGE_USE_k4-&INIT_AGE_USE_k3);

GENDER_MATCH=0;
if GENDER=T_GENDER_DON then GENDER_MATCH=1;
T_AGE_DON_B=0;
if T_AGE_DON>=55 then T_AGE_DON_B=1;
run;

%macro prob_trans;
%if (&group_var="A" or &group_var="D" ) %then %do;
proc phreg data=phreg_data NOSUMMARY ;
class T_TRR_ID_CODE;
model dummy_time*transplanted(0) = 
				match_LAS match_LAS1 match_LAS2 /*baseline effect of LAS*/
				mov_sum mov_sum1 mov_sum2 /*baseline effect of center volume */
				list_days list_days1 list_days2 /*baseline effect of time on waiting list */
				hgt_diff hgt_diff1 hgt_diff2 /*baseline effect of size mismatch */
				init_age_use  /*baseline effect of recipeient age */
				bilateral*init_age_use bilateral*match_las /*interaction of bilateral with age and match_las*/
				T_AGE_DON_B*init_age_use /*age of donor and recipeint */
				T_HIST_CIG_DON_Y*init_age_use /*donor cigarette and recipient age interaction */
				/ ties=breslow;
strata T_TRR_ID_CODE; 
ods output ParameterEstimates=parm_est;
output  out=final_data_pred xbeta=linear_pred;
run;
%end;

%if (&group_var="C" or &group_var="B") %then %do;
proc phreg data=phreg_data NOSUMMARY ;
class T_TRR_ID_CODE;
model dummy_time*transplanted(0) = 
				match_LAS match_LAS1 match_LAS2 /*baseline effect of LAS*/
				mov_sum mov_sum1 mov_sum2 /*baseline effect of center volume */
				list_days list_days1 list_days2 /*baseline effect of time on waiting list */
				hgt_diff hgt_diff1 hgt_diff2 /*baseline effect of size mismatch */
				init_age_use  /*baseline effect of recipeient age */
				T_AGE_DON_B*init_age_use /*age of donor and recipeint */
				T_HIST_CIG_DON_Y*init_age_use /*donor cigarette and recipient age interaction */
				/ ties=breslow;
strata T_TRR_ID_CODE; 
*id WL_ID_CODE;
ods output ParameterEstimates=parm_est;
output  out=final_data_pred xbeta=linear_pred;
run;
%end;
%mend;

%prob_trans;

data parm_est; set parm_est;
drop df HazardRatio Label;
zscore = sqrt(Chisq); run;

data parm_est; set parm_est;
drop ChiSq; run;

proc export data=parm_est
 outfile="xtable/treatment_allocation_model_GROUP_&group_var1._r2.txt" dbms=tab replace; run; 

data final_data_pred1; set final_data_pred (keep=T_TRR_ID_CODE WL_ID_CODE linear_pred);
e_linear_pred = exp(linear_pred); run;

proc means data=final_data_pred1 sum noprint;
var e_linear_pred;
by T_TRR_ID_CODE; 
output out=sum_linear_pred sum=s_e_linear_pred;  run;

data sum_linear_pred; set sum_linear_pred (keep=T_TRR_ID_CODE s_e_linear_pred);
run;

data final_data_pred1;
merge final_data_pred1 sum_linear_pred;
by T_TRR_ID_CODE; 
PROB_TRANS = e_linear_pred/s_e_linear_pred;
run;

data final_data_pred1; 
set final_data_pred1 (keep= T_TRR_ID_CODE WL_ID_CODE PROB_TRANS);
run;

proc sort data=final_data_pred1;
by T_TRR_ID_CODE WL_ID_CODE; run;
proc sort data=total_final1;
by T_TRR_ID_CODE WL_ID_CODE; run;

data total_final_prob; merge total_final1 final_data_pred1; 
by T_TRR_ID_CODE WL_ID_CODE; run;

/*proc contents data=total_final_prob; run;*/

data total_final_prob; set total_final_prob;
if mov_sum=. then mov_sum=0;
T_BILATERAL = 0;
IF T_TX_PROCEDUR_TY=603 then T_BILATERAL=1;
IF T_TX_PROCEDUR_TY=604 then T_BILATERAL=1;
IF TX_PROCEDUR_TY1=603 or TX_PROCEDUR_TY1=604 then BILATERAL=1;
else BILATERAL=0;
MATCH_LAS1 = max((MATCH_LAS - &match_las_k1)**3/(&match_las_k4-&match_las_k1)**2, 0)+
	(&match_las_k3-&match_las_k1)*max((MATCH_LAS - &match_las_k4)**3/(&match_las_k4-&match_las_k1)**2, 0)/(&match_las_k4-&match_las_k3)-
	(&match_las_k4-&match_las_k1)*max((MATCH_LAS - &match_las_k3)**3/(&match_las_k4-&match_las_k1)**2, 0)/(&match_las_k4-&match_las_k3);
MATCH_LAS2 = max((MATCH_LAS - &match_las_k2)**3/(&match_las_k4-&match_las_k1)**2, 0)+
	(&match_las_k3-&match_las_k2)*max((MATCH_LAS - &match_las_k4)**3/(&match_las_k4-&match_las_k1)**2, 0)/(&match_las_k4-&match_las_k3)-
	(&match_las_k4-&match_las_k2)*max((MATCH_LAS - &match_las_k3)**3/(&match_las_k4-&match_las_k1)**2, 0)/(&match_las_k4-&match_las_k3);
MOV_SUM1 = max((MOV_SUM - &mov_sum_k1)**3/(&mov_sum_k4-&mov_sum_k1)**2, 0)+
	(&mov_sum_k3-&mov_sum_k1)*max((MOV_SUM - &mov_sum_k4)**3/(&mov_sum_k4-&mov_sum_k1)**2, 0)/(&mov_sum_k4-&mov_sum_k3)-
	(&mov_sum_k4-&mov_sum_k1)*max((MOV_SUM - &mov_sum_k3)**3/(&mov_sum_k4-&mov_sum_k1)**2, 0)/(&mov_sum_k4-&mov_sum_k3);
MOV_SUM2 = max((MOV_SUM - &mov_sum_k2)**3/(&mov_sum_k4-&mov_sum_k1)**2, 0)+
	(&mov_sum_k3-&mov_sum_k2)*max((MOV_SUM - &mov_sum_k4)**3/(&mov_sum_k4-&mov_sum_k1)**2, 0)/(&mov_sum_k4-&mov_sum_k3)-
	(&mov_sum_k4-&mov_sum_k2)*max((MOV_SUM - &mov_sum_k3)**3/(&mov_sum_k4-&mov_sum_k1)**2, 0)/(&mov_sum_k4-&mov_sum_k3);
LIST_DAYS1 = max((LIST_DAYS - &LIST_DAYS_k1)**3/(&LIST_DAYS_k4-&LIST_DAYS_k1)**2, 0)+
	(&LIST_DAYS_k3-&LIST_DAYS_k1)*max((LIST_DAYS - &LIST_DAYS_k4)**3/(&LIST_DAYS_k4-&LIST_DAYS_k1)**2, 0)/(&LIST_DAYS_k4-&LIST_DAYS_k3)-
	(&LIST_DAYS_k4-&LIST_DAYS_k1)*max((LIST_DAYS - &LIST_DAYS_k3)**3/(&LIST_DAYS_k4-&LIST_DAYS_k1)**2, 0)/(&LIST_DAYS_k4-&LIST_DAYS_k3);
LIST_DAYS2 = max((LIST_DAYS - &LIST_DAYS_k2)**3/(&LIST_DAYS_k4-&LIST_DAYS_k1)**2, 0)+
	(&LIST_DAYS_k3-&LIST_DAYS_k2)*max((LIST_DAYS - &LIST_DAYS_k4)**3/(&LIST_DAYS_k4-&LIST_DAYS_k1)**2, 0)/(&LIST_DAYS_k4-&LIST_DAYS_k3)-
	(&LIST_DAYS_k4-&LIST_DAYS_k2)*max((LIST_DAYS - &LIST_DAYS_k3)**3/(&LIST_DAYS_k4-&LIST_DAYS_k1)**2, 0)/(&LIST_DAYS_k4-&LIST_DAYS_k3);
HGT_DIFF1 = max((HGT_DIFF - &HGT_DIFF_k1)**3/(&HGT_DIFF_k4-&HGT_DIFF_k1)**2, 0)+
	(&HGT_DIFF_k3-&HGT_DIFF_k1)*max((HGT_DIFF - &HGT_DIFF_k4)**3/(&HGT_DIFF_k4-&HGT_DIFF_k1)**2, 0)/(&HGT_DIFF_k4-&HGT_DIFF_k3)-
	(&HGT_DIFF_k4-&HGT_DIFF_k1)*max((HGT_DIFF - &HGT_DIFF_k3)**3/(&HGT_DIFF_k4-&HGT_DIFF_k1)**2, 0)/(&HGT_DIFF_k4-&HGT_DIFF_k3);
HGT_DIFF2 = max((HGT_DIFF - &HGT_DIFF_k2)**3/(&HGT_DIFF_k4-&HGT_DIFF_k1)**2, 0)+
	(&HGT_DIFF_k3-&HGT_DIFF_k2)*max((HGT_DIFF - &HGT_DIFF_k4)**3/(&HGT_DIFF_k4-&HGT_DIFF_k1)**2, 0)/(&HGT_DIFF_k4-&HGT_DIFF_k3)-
	(&HGT_DIFF_k4-&HGT_DIFF_k2)*max((HGT_DIFF - &HGT_DIFF_k3)**3/(&HGT_DIFF_k4-&HGT_DIFF_k1)**2, 0)/(&HGT_DIFF_k4-&HGT_DIFF_k3);
INIT_AGE_USE1 = max((INIT_AGE_USE - &INIT_AGE_USE_k1)**3/(&INIT_AGE_USE_k4-&INIT_AGE_USE_k1)**2, 0)+
	(&INIT_AGE_USE_k3-&INIT_AGE_USE_k1)*max((INIT_AGE_USE - &INIT_AGE_USE_k4)**3/(&INIT_AGE_USE_k4-&INIT_AGE_USE_k1)**2, 0)/(&INIT_AGE_USE_k4-&INIT_AGE_USE_k3)-
	(&INIT_AGE_USE_k4-&INIT_AGE_USE_k1)*max((INIT_AGE_USE - &INIT_AGE_USE_k3)**3/(&INIT_AGE_USE_k4-&INIT_AGE_USE_k1)**2, 0)/(&INIT_AGE_USE_k4-&INIT_AGE_USE_k3);
INIT_AGE_USE2 = max((INIT_AGE_USE - &INIT_AGE_USE_k2)**3/(&INIT_AGE_USE_k4-&INIT_AGE_USE_k1)**2, 0)+
	(&INIT_AGE_USE_k3-&INIT_AGE_USE_k2)*max((INIT_AGE_USE - &INIT_AGE_USE_k4)**3/(&INIT_AGE_USE_k4-&INIT_AGE_USE_k1)**2, 0)/(&INIT_AGE_USE_k4-&INIT_AGE_USE_k3)-
	(&INIT_AGE_USE_k4-&INIT_AGE_USE_k2)*max((INIT_AGE_USE - &INIT_AGE_USE_k3)**3/(&INIT_AGE_USE_k4-&INIT_AGE_USE_k1)**2, 0)/(&INIT_AGE_USE_k4-&INIT_AGE_USE_k3);
T_GENDER_MATCH=0;
if GENDER=T_GENDER_DON then T_GENDER_MATCH=1;
if GENDER=GENDER_DON1 then GENDER_MATCH=1;
else GENDER_MATCH=0;
T_AGE_DON_B=0;
if T_AGE_DON>=55 then T_AGE_DON_B=1;
if AGE_DON1>=55 then AGE_DON_B=1;
else AGE_DON_B=0;
if COD_CAD_DON1=3 then COD_CAD_DON_HT=1;
else COD_CAD_DON_HT=0;
if HIST_CIG_DON1="Y" then HIST_CIG_DON_Y=1;
else HIST_CIG_DON_Y=0;
rename HGT_DIFF = T_HGT_DIFF;
HGT_DIFF1 = HGT_CM_TCR-HGT_CM_DON_CALC1;
run;


data total_final_prob1; set total_final_prob;
drop ABO ABO_DON AGE BEGIN_DT BMI_TCR CALC_LAS CITIZENSHIP CMV_DON COD_CAD_DON DEATH_DATE DEATH_DATE_COMBINED DIABETES_DON DUMMY_TIME ECD_DONOR END_DATE END_DT ETHCAT 
FAILURE FINAL_HGT GENDER GENDER_DON HGT_DIFF1 HGT_DIFF2 HIST_CIG_DON INIT_AGE INIT_AGE_USE1 INIT_AGE_USE2 LAST_AGE LISTING_CTR_CODE LIST_DAYS1 LIST_DAYS2 PO2
PRI_PAYMENT_TCR PX_STAT PX_STAT_DATE RACE RACE_DON REGION REM_CD RETXDATE RETXDATE1 SHARE_TY SSDMF_DEATH_DATE TCR_DGN TX_DATE TX_PROCEDUR_TY T_ABO_DON T_CMV_DON T_CMV_DON_P
T_COD_CAD_DON T_DIABETES_DON T_DIABETES_DON_Y T_ECD_DONOR T_GENDER_DON T_HGT_CM_DON_CALC T_HIST_CIG_DON T_PO2 T_RACE_DON cohort trans
TCR_DGN_OSTXT; run;

data total_final_prob1; set total_final_prob1;
drop AGE1 AGE_DON BMI_DON_CALC BMI_DON_CALC1 CMV_DON1 DEATH_DATE_COMBINED_USE DEATH_DATE_SOURCE DISTANCE ECD_DONOR1 ETHCAT_DON ETHCAT_DON1 HGT_CM_DON_CALC ISCHTIME1 MATCH_LAS1 MATCH_LAS2
MOV_SUM1 MOV_SUM2 PO21 RACE_DON1 RETXDATE_USE; run;

proc export data=total_final_prob1 outfile="data/combined_longitudinal_listing_prob_trans_GROUP_&group_var1._r2.txt" dbms=tab replace; run;   

