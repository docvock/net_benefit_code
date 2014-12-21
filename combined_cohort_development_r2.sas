/*
#######################################################
### SAS code to generate the groups of patients eligeable for 
### at each organ transplant and their estimated probability of transplant.
### In this analysis we will subdivide by native disease grouping
###
### PROGRAMS RELIED ON: none
### SAS DATA FILES USED: thoracic_public_use_data 
###                      thoracic_las_audit_data
###                      thoracic_las_history_data
### ADDITIONAL PROGRAMS USED: experiment_generation.sas
###                           prob_generation.sas
### SAS DATA FILES GENERATED: public_nodup_final (one observation per patient in the cohort for future analysis - 
###								generally the last listing for each patient)
### TXT DATA FILES GENERATED: combined_all_listings (one observation per listing in the cohort)
###                           combined_patient_record (one observation per patient in the cohort)
###                           combined_longitudinal_listing_prob_trans_GROUP_&group_var
#######################################################
*/

options formdlim="_" nodate nonumber;
x 'cd C:\\Users\\bstvock\\Documents\\research\\net_benefit_CF\\project_files';

Libname newlib 'C:\\Users\\bstvock\\Documents\\UNOS data\\updated_data_12_31_2011';  
Libname mylib 'data/'; 

/*************************************************************/
/* read-in data */
/*************************************************************/
data public;
set newlib.thoracic_public_use_data (keep=
PT_CODE WL_ID_CODE TRR_ID_CODE /*patient identifiers*/
citizenship init_age ABO BMI_TCR PRI_PAYMENT_TCR ETHCAT RACE REGION GENDER HGT_CM_TCR TCR_DGN TCR_DGN_OSTXT /*baseline patient characteristics*/
init_date init_stat organ WL_ORG WLHL WLHR WLIN WLKI WLKP WLLI WLPA WLPI listing_ctr_code init_blu_flg init_llu_flg init_rlu_flg/*Listing information*/
AGE end_match_las TX_DATE TXED TXHRT TXINT TXKID TXLIV TXPAN TX_PROCEDUR_TY ABO_DON ABO_MAT AGE AGE_DON ECD_DONOR PO2 COD_CAD_DON 
CONTIN_CIG_DON HIST_CIG_DON ISCHTIME SHARE_TY HGT_CM_DON_CALC GENDER_DON CMV_DON RACE_DON ETHCAT_DON DIABETES_DON DISTANCE END_HGT_CM ECD_DONOR BMI_DON_CALC
/*transplant and donor characteristics*/
PX_STAT PX_STAT_DATE  death_date end_date REM_CD RETXDATE SSDMF_DEATH_DATE /*status information*/ 
); run;

data LAS_audit; 
set newlib.thoracic_las_audit_data; run;

data LAS_hist_comp;
set newlib.thoracic_las_history_data; run;

/************************************************************/
/* Combine LAS grouping covariate with audit data */
/************************************************************/
proc sort data=LAS_hist_comp; by WL_ID_CODE; run;

data LAS_hist; set LAS_hist_comp (keep=WL_ID_CODE GROUPING); 
by WL_ID_CODE; f_WL_ID = first.WL_ID_CODE; run;

data LAS_hist; set LAS_hist; where f_WL_ID=1; run;

data LAS_hist; set LAS_hist; drop f_WL_ID; run;

proc sort data=LAS_audit; by WL_ID_CODE; run; 

data LAS_audit; merge LAS_audit LAS_hist;
by WL_ID_CODE; run;

/*************************************************************/
/* Get grouping variable to the pulbic dataset and initial non-zero LAS*/
/*************************************************************/
data group_var; set LAS_audit; where MATCH_LAS>0; 
INIT_DATE_LAS = datepart(BEGIN_DT); run;

proc sort data=group_var; by WL_ID_CODE INIT_DATE_LAS descending BEGIN_DT; run;

data group_var; set group_var; by WL_ID_CODE INIT_DATE_LAS descending BEGIN_DT;
first_WL = first.WL_ID_CODE; run; 

data group_var;
set group_var (keep=WL_ID_CODE GROUPING MATCH_LAS INIT_DATE_LAS first_WL); where first_WL=1; run;

data group_var;
set group_var; 
rename MATCH_LAS=INIT_MATCH_LAS; run;

proc sort data=public; by WL_ID_CODE; run;

proc sort data=group_var; by WL_ID_CODE; run;

data public; merge public group_var; by WL_ID_CODE; run;

data public; set public; drop first_WL; run;

/************************************************************/
/* Combined death date for all subjects based on lung and heart/lung transplants data */
/************************************************************/

data public_death_date; set public (keep=PT_CODE WL_ID_CODE PX_stat PX_stat_date death_date SSDMF_death_date init_stat);
where INIT_STAT=7010 or INIT_STAT=7999 or INIT_STAT=1010 or INIT_STAT=1020 or INIT_STAT=1030 or INIT_STAT=1090 or INIT_STAT=1999;    
PX_death_date = .;
if PX_stat="D" then PX_death_date=PX_stat_date;
death_date_avail=0;
if death_date ne . then death_date_avail=1;
PX_death_date_avail =0;
if PX_death_date ne . then PX_death_date_avail=1;
SSDMF_death_date_avail = 0;
if SSDMF_death_date ne . then SSDMF_death_date_avail=1; 
run;

proc sort data=public_death_date out=public_death_date;
by PT_CODE descending death_date_avail death_date
descending PX_death_date_avail PX_death_date
descending SSDMF_death_date_avail SSDMF_death_date;
run;

data public_death_date; set public_death_date; by PT_CODE;
first_patient = first.PT_CODE;  run;

data public_death_date_s; set public_death_date;
where first_patient=1; run;

data public_death_date_s; set public_death_date_s;
DEATH_DATE_SOURCE = 0;
DEATH_DATE_COMBINED = death_date;
if death_date ne . then death_date_source=1;
if death_date= . then death_date_combined=PX_death_date;
if death_date= . & PX_death_date ne . then death_date_source=2;
if death_date= . & PX_death_date=. then death_date_combined=SSDMF_death_date;
if death_date= . & PX_death_date=. & SSDMF_death_date ne . then death_date_source=3; run;

proc sort data=public; by PT_CODE; run;

data public; merge public_death_date_s(keep=PT_CODE death_date_combined death_date_source) public;
by PT_CODE; run;

/*************************************************************/
/* Determine the moving average center volume of all lung (only) adult transplants */
/*************************************************************/
data public_trans; set public (keep=PT_CODE AGE listing_ctr_code tx_date trr_id_code organ);
where tx_date ne . & organ="LU"; run;

data public_trans; set public_trans;
where age ge 18 & tx_date ge '04MAY2003'd & tx_date le '30SEP2011'd; 
trans=1; run; 

proc sort data=public_trans; by listing_ctr_code tx_date;  run;

proc means data=public_trans N noprint;
var trans;
by listing_ctr_code tx_date; 
output out=public_trans_tot; run;

data public_trans_tot; set public_trans_tot;
where _STAT_="N"; run;

data listing_center; set public_trans_tot; by listing_ctr_code;
list_first = first.listing_ctr_code; run;

data listing_center; set listing_center;
where list_first=1; run;

data listing_center; set listing_center;
keep listing_ctr_code; run;

data listing_center; set listing_center;
num_list_ctr=_N_; run;

data list_dates;
do num_list_ctr=1 to 71;
do tx_date='04MAY2003'd to '30SEP2011'd; output; end; end; run;

data list_dates; merge list_dates listing_center;
by num_list_ctr; run;

data list_dates; set list_dates; drop num_list_ctr; run;

data public_trans_tot; set public_trans_tot;
drop _TYPE_ _FREQ_ _STAT_; run;

data public_trans_tot; merge public_trans_tot list_dates;
by listing_ctr_code tx_date; run;

data public_trans_tot; set public_trans_tot; 
if trans=. then trans=0; run;

%let n = 731;

data public_trans_tot;
  set public_trans_tot;
  by listing_ctr_code;
  retain num_sum 0;
  if first.listing_ctr_code then do;
    count=0;
    num_sum=0;
  end;
  count+1;
  last&n=lag&n(trans);
  if count gt &n then num_sum=sum(num_sum,trans,-last&n);
  else num_sum=sum(num_sum,trans);
  if count ge &n then mov_sum=num_sum-trans;
  else mov_sum=.;
run;

data public_trans_tot_s; set public_trans_tot;
where tx_date ge '04MAY2005'd; drop num_sum count last731; run; 

data public_trans_tot_s; set public_trans_tot_s;
where LISTING_CTR_CODE ne ""; run;

/*************************************************************/
/* Limit Analysis to Patients that have not had a prior lung transplant */
/*************************************************************/

data public_retrans; set public (keep=PT_CODE tx_date trr_id_code organ);
where tx_date ne . & (organ="LU" or organ="HL"); run;

data public_retrans; set public_retrans; 
where TX_DATE le '30SEP2011'd; run; 

proc sort data=public_retrans;
by PT_CODE tx_date; run;

data public_retrans;
set public_retrans;
by PT_CODE;
first_organ = first.PT_CODE; run;

data LAS_final_organ1;
set public_retrans;
where first_organ=1;
rename tx_date = TX_DATE1; 
run;

data LAS_final_organ1; 
set LAS_final_organ1 (keep=PT_CODE tx_date1); run;

data public;
merge public LAS_final_organ1; 
by PT_CODE; run;

/*************************************************************/
/* restrict analysis using the inclusion/exclusion criteria described */
/*************************************************************/

data public1; set public;
where init_stat=7010 or init_stat=7999; 
LIST_TRANS = 0;
if TRR_ID_CODE ne "" & TX_DATE le '30SEP2011'd then LIST_TRANS=1; run;

data public1; set public1;
where end_date ge '4MAY2005'd & init_date le '30SEP2011'd ; run;

/*proc freq data=public1; tables LIST_TRANS; run;*/

data public1; set public1;
PRIOR_TX = 0;
if TX_DATE1 ne . & INIT_DATE > TX_DATE1 then PRIOR_TX=1; 
if TX_DATE1 ne . & INIT_DATE = TX_DATE1 & TX_DATE > TX_DATE1 then PRIOR_TX=1;
if TX_DATE1 ne . & INIT_DATE = TX_DATE1 & TX_DATE =. then PRIOR_TX=1;
run;

data public1; set public1;
where PRIOR_TX=0; run;

/*proc freq data=public1; tables LIST_TRANS; run;*/

data public1; set public1;
where WLHR="" & WLHL="" & WLIN="" & WLKI="" & WLKP="" & WLLI="" & WLPA="" & WLPI="" &
PT_CODE ne 125779 & PT_CODE ne 619743 & PT_CODE ne 632410 & PT_CODE ne 638611 & PT_CODE ne 644831 &
PT_CODE ne 663088 & PT_CODE ne 665398 & PT_CODE ne 734998 & PT_CODE ne 750685 & PT_CODE ne 778780 & 
PT_CODE ne 792102 & PT_CODE ne 871377 & PT_CODE ne 29495 & PT_CODE ne 70683 & PT_CODE ne 121742 & 
PT_CODE ne 474489 & PT_CODE ne 692711 & PT_CODE ne 786442; run;
/* PT_CODE 121742, 125779, 474489, 692711, 70683, and 786442 listed initially for lung only and then lung+liver subsequently */
/* PT_CODE 644831, 663088, 665398, 734998, 750685, 778780, 792102, and 871377 listed initially for lung only and then lung+hear subsequently */
/* PT_CODE 29495, 6197432, 632410, 638611 concurrently listed for heart/lung double transplant at same cetner */

/*proc freq data=public1; tables LIST_TRANS; run;*/

data public1; set public1;
where init_age ge 18 or ((init_date_las-init_date)/365+init_age ge 18); run;

/*proc freq data=public1; tables LIST_TRANS; run;*/

data public1; set public1;
nonliving = 1;
if TRR_ID_CODE ne "" & TX_DATE le '30SEP2011'd & TXED=0 then nonliving=0; 
if REM_CD=15 & END_DATE le '30SEP2011'd then nonliving=0; run;

data public1; set public1;
where nonliving=1; run;

/*proc freq data=public1; tables LIST_TRANS; run;*/ 

data public1; set public1;
where citizenship=1 & PT_CODE ne 635399 & PT_CODE ne 778539 & PT_CODE ne 862252 & PT_CODE ne 669266; run;
/* PT_CODE=669266, 635399, 778539, and 862252 was listed once as citizen and then once as permanent resident*/

/*proc freq data=public1; tables LIST_TRANS; run;*/ 

data public1; set public1;
where TRR_ID_CODE = "" or END_MATCH_LAS > 0; run;

data public1; set public1;
where INIT_DATE_LAS ne . & INIT_DATE_LAS le '30SEP2011'd; run;

data public1; set public1;
if GROUPING="" then GROUPING="D"; run;
/* PT_CODE=573640 has missing grouping but TCR_DGN=422 which is always assigned group D */

proc freq data=public1; tables LIST_TRANS; run; 
proc freq data=public1; tables LIST_TRANS*GROUPING; run; 

/**************************************************************************/
/* get the transplant and approximate match date for all transplanted listing */

data listings_all; set public1; run;

proc sort data=listings_all;
by WL_ID_CODE; run;

proc sort data=LAS_audit;
by WL_ID_CODE; run;

data LAS; merge LAS_audit listings_all; 
by WL_ID_CODE; run;

data LAS_trans; set LAS;
where LIST_TRANS=1; run;

/*Figure out which date corresponds to match listing */
data LAS_trans; 
set LAS_trans;
LAS_diff=end_match_las-match_las; run;

data LAS_final_sub_1a;
set LAS_trans; where LAS_diff=0; 
approx_match_date=datepart(begin_dt);
diff_match_tx=tx_date-approx_match_date;
after_tx = 0;
if approx_match_date > tx_date then after_tx=1; run;

proc sort data=LAS_final_sub_1a;
by WL_ID_CODE after_tx descending begin_dt; run;

data LAS_final_sub_1a;
set LAS_final_sub_1a; 
by WL_ID_CODE;
last_day=first.WL_ID_CODE; run;

data LAS_final_sub_1a; set LAS_final_sub_1a; 
where last_day=1; run;

proc sort data=LAS_final_sub_1a;
by PT_CODE approx_match_date; run;

data LAS_final_sub_1a;
set LAS_final_sub_1a;
by PT_CODE;
first_organ = first.PT_CODE; run;

data LAS_final_organ1;
set LAS_final_sub_1a;
where first_organ=1;
rename approx_match_date = APPROX_MATCH_DATE1;
rename tx_date = tx_date_alt; 
rename RETXDATE=RETXDATE1;
rename END_MATCH_LAS=T_MATCH_LAS1;
rename TX_PROCEDUR_TY=TX_PROCEDUR_TY1; 
rename AGE_DON=AGE_DON1;
rename COD_CAD_DON=COD_CAD_DON1;
rename HIST_CIG_DON=HIST_CIG_DON1;
rename HGT_CM_DON_CALC=HGT_CM_DON_CALC1;
rename GENDER_DON=GENDER_DON1;
rename DIABETES_DON = DIABETES_DON1;
rename PO2 = PO21;
rename RACE_DON=RACE_DON1;
rename ECD_DONOR = ECD_DONOR1;
rename CMV_DON = CMV_DON1;
rename BMI_DON_CALC = BMI_DON_CALC1;
rename ISCHTIME = ISCHTIME1;
rename AGE = AGE1;
rename ETHCAT_DON = ETHCAT_DON1;
run;

data LAS_final_organ_not1;
set LAS_final_sub_1a;
where first_organ=0; run;

data LAS_final_organ1; 
set LAS_final_organ1 (keep=PT_CODE approx_match_date1 tx_date_alt RETXDATE1 T_MATCH_LAS1 TX_PROCEDUR_TY1 AGE_DON1 COD_CAD_DON1 HIST_CIG_DON1 HGT_CM_DON_CALC1 GENDER_DON1
DIABETES_DON1 PO21 RACE_DON1 ECD_DONOR1 CMV_DON1 BMI_DON_CALC1 ISCHTIME1 AGE1 ETHCAT_DON1); run;

data public1;
merge public1 LAS_final_organ1; 
by PT_CODE; run;

data public1; set public1;
tx_date_diff = tx_date1-tx_date_alt; 
tx_avail1 = 0;
tx_avail_alt = 0;
if tx_date1 ne . then tx_avail1=1;
if tx_date_alt ne . then tx_avail_alt=1;
run;

proc univariate data=public1; var tx_date_diff; run;
proc freq data=public1; tables tx_avail1*tx_avail_alt; run;

data public1; set public1;
TRANSPLANTED = 0;
if tx_date1 ne . then TRANSPLANTED=1; run;

data public_nodup;
set public1; where (TRANSPLANTED=0) or (INIT_DATE le tx_date1); run;

proc sort data=public_nodup;
by PT_CODE INIT_DATE; run;

data public_nodup; set public_nodup;
END_DATE_USE = min(end_date,approx_match_date1,death_date_combined,'30SEP2011'd); 
REM_TRANS=0;
REM_DETER=0;
REM_HEALTHY=0;
REM_OTHER=0;
REM_DEATH=0;
if end_date le '30SEP2011'd then do;
if REM_CD=4 then REM_TRANS=1;
if REM_CD=8 then REM_DEATH=1;
if REM_CD=13 then REM_DETER=1;
if REM_CD=12 then REM_HEALTHY=1;
if REM_CD=6 or REM_CD=9 or REM_CD=16 or REM_CD=24 then REM_OTHER=1;
end; run;

proc sort data=public_nodup; 
by PT_CODE descending end_date_use descending REM_TRANS descending REM_DEATH descending REM_DETER descending REM_HEALTHY descending REM_OTHER descending END_MATCH_LAS descending INIT_DATE; run;

data public_nodup; set public_nodup; by PT_CODE;
first_PT = first.PT_CODE; run;

data public_nodup_final; set public_nodup;
where first_PT=1; run;

proc sort data=public_nodup; 
by PT_CODE init_date_las descending REM_TRANS descending REM_DEATH descending REM_DETER descending REM_HEALTHY descending REM_OTHER; run;

data public_nodup; set public_nodup; by PT_CODE;
first_PT = first.PT_CODE; run;

data public_nodup1; set public_nodup;
where first_PT=1; run;

/* get covariate information from first listing */
data public_nodup1; set public_nodup1 (keep=PT_CODE WL_ID_CODE INIT_DATE_LAS INIT_MATCH_LAS INIT_DATE INIT_AGE LISTING_CTR_CODE TX_DATE1 APPROX_MATCH_DATE1);
rename INIT_DATE = INIT_DATE1;
rename INIT_DATE_LAS=INIT_DATE_LAS1;
rename INIT_MATCH_LAS=INIT_MATCH_LAS1; 
INIT_AGE_LAS = INIT_AGE+(INIT_DATE_LAS-INIT_DATE)/365;
cohort=1;
run;

/* add in center volume at listing and transplantation */
data public_trans_tot_m1;
set public_trans_tot_s;
rename mov_sum=CENTER_VOL1;
drop trans;
rename TX_DATE = INIT_DATE_LAS1; run;

proc sort data=public_nodup1; by LISTING_CTR_CODE INIT_DATE_LAS1; run;

data public_nodup1; merge public_nodup1 public_trans_tot_m1;
by LISTING_CTR_CODE INIT_DATE_LAS1; run;

data public_nodup1; set public_nodup1;
where cohort=1; run;

proc sort data=public_nodup1; by LISTING_CTR_CODE TX_DATE1; run;

data public_trans_tot_m2;
set public_trans_tot_s;
rename mov_sum=T_CENTER_VOL;
drop trans;
rename TX_DATE = TX_DATE1; run;

data public_nodup1; merge public_nodup1 public_trans_tot_m2;
by LISTING_CTR_CODE TX_DATE1; run;

data public_nodup1; set public_nodup1;
where cohort=1; run;

proc print data=public_nodup1 (obs=20); run;

/* Add in clinical characteristics at listing and transplantation 
from LAS history datafile*/
proc sort data=public_nodup1; by WL_ID_CODE; run;

proc sort data=LAS_hist_comp; by WL_ID_CODE; run;

data public_nodup1s; set public_nodup1 (keep= WL_ID_CODE INIT_DATE_LAS1 cohort);

proc contents data=LAS_hist_comp;
run;

data LAS_hist_list; merge LAS_hist_comp(keep= WL_ID_CODE CHG_DATE CALC_BMI
CALC_FUNC_STAT
CALC_FVC_PRE
CALC_NYHA
CALC_O2
CALC_SIX_MIN_WALK
CALC_VENT_USE
DIABETES) public_nodup1s;
by WL_ID_CODE; run;

data LAS_hist_list; set LAS_hist_list; where cohort=1 & CHG_DATE le INIT_DATE_LAS1; run;

data LAS_hist_list; set LAS_hist_list; by WL_ID_CODE;
last_WL =last.WL_ID_CODE; run;

data LAS_hist_list; set LAS_hist_list; where last_WL=1; 
rename CALC_BMI= CALC_BMI1;
rename CALC_FUNC_STAT=CALC_FUNC_STAT1;
rename CALC_FVC_PRE=CALC_FVC_PRE1;
rename CALC_NYHA=CALC_NYHA1;
rename CALC_O2= CALC_O21;
rename CALC_SIX_MIN_WALK= CALC_SIX_MIN_WALK1;
rename CALC_VENT_USE= CALC_VENT_USE1;
rename DIABETES= DIABETES1;
drop CHG_DATE cohort INIT_DATE_LAS1;
run; 

data public_nodup1; merge public_nodup1 LAS_hist_list;
by WL_ID_CODE; run;

/********/
data public_nodup1s; set public_nodup1 (keep= WL_ID_CODE APPROX_MATCH_DATE1 cohort); run;

data LAS_hist_trans; merge LAS_hist_comp(keep= WL_ID_CODE CHG_DATE CALC_BMI
CALC_FUNC_STAT
CALC_FVC_PRE
CALC_NYHA
CALC_O2
CALC_SIX_MIN_WALK
CALC_VENT_USE
DIABETES) public_nodup1s;
by WL_ID_CODE; run;

data LAS_hist_trans; set LAS_hist_trans; where cohort=1 & CHG_DATE le APPROX_MATCH_DATE1 & APPROX_MATCH_DATE1 ne .; run;

data LAS_hist_trans; set LAS_hist_trans; by WL_ID_CODE;
last_WL =last.WL_ID_CODE; run;

data LAS_hist_trans; set LAS_hist_trans; where last_WL=1; 
rename CALC_BMI= T_CALC_BMI;
rename CALC_FUNC_STAT= T_CALC_FUNC_STAT;
rename CALC_FVC_PRE=T_CALC_FVC_PRE;
rename CALC_NYHA=T_CALC_NYHA;
rename CALC_O2= T_CALC_O2;
rename CALC_SIX_MIN_WALK= T_CALC_SIX_MIN_WALK;
rename CALC_VENT_USE= T_CALC_VENT_USE;
rename DIABETES= T_DIABETES;
drop CHG_DATE cohort APPROX_MATCH_DATE1;
run; 


data public_nodup1; merge public_nodup1 LAS_hist_trans;
by WL_ID_CODE; run;

proc sort data=public_nodup1; by PT_CODE; run;

data public_nodup1; set public_nodup1;
drop WL_ID_CODE INIT_AGE LISTING_CTR_CODE TX_DATE1 APPROX_MATCH_DATE1; run;

/**************/

data public_nodup_final; merge public_nodup_final public_nodup1; 
by PT_CODE; run;

data public_nodup_final; set public_nodup_final;
FAILURE = 0;
if DEATH_DATE_COMBINED ne . & DEATH_DATE_COMBINED le '30SEP2011'd & DEATH_DATE_COMBINED ge '01JAN2001'd then FAILURE=1; run;

data public_nodup_final; set public_nodup_final;
if REM_CD=21 then REM_CD=4;
if REM_CD=5 then REM_CD=13;
if REM_CD=14 then REM_CD=12;
if REM_CD=16 then REM_CD=12; run;

/* Creates dataset with one observation per listing with information from generally the last listing time if 
listed and delisted multiple times. If listed concurrently at more than one center, this uses center where
patient was transplanted. Need this dataset for analysis of quality of SSDMF done in another files */
 data mylib.public_nodup_final;
set public_nodup_final; run;



data public1; set public1;
END_DATE_USE = min(end_date,approx_match_date1,death_date_combined,'30SEP2011'd); 
FAILURE = 0;
if DEATH_DATE_COMBINED ne . & DEATH_DATE_COMBINED le '30SEP2011'd & DEATH_DATE_COMBINED ge '01JAN2001'd then FAILURE=1;
drop ABO_MAT INIT_STAT ISCHTIME ORGAN PRIOR_TX TXED TXHRT TXINT TXKID TXLIV TXPAN WLHL WLHR WLIN WLKI WLKP WLLI WLPA WLPI WL_ORG nonliving tx_avail1 tx_avail_alt
tx_date_alt tx_date_diff CONTIN_CIG_DON INIT_BLU_FLG INIT_LLU_FLG INIT_RLU_FLG; run;

data public_nodup_final; set public_nodup_final;
*END_DATE_USE = min(end_date,approx_match_date1,death_date_combined,'30SEP2011'd); 
drop ABO_MAT INIT_STAT ORGAN PRIOR_TX TXED TXHRT TXINT TXKID TXLIV TXPAN WLHL WLHR WLIN WLKI WLKP WLLI WLPA WLPI WL_ORG nonliving tx_avail1 tx_avail_alt
tx_date_alt tx_date_diff CONTIN_CIG_DON INIT_BLU_FLG INIT_LLU_FLG INIT_RLU_FLG list_trans REM_DEATH REM_HEALTHY REM_OTHER REM_TRANS FIRST_PT cohort last_WL; run;

data public_nodup_final; set public_nodup_final;
EVENT_DATE_USE = '30SEP2011'd;
EVENT_INDICATOR = 0;
DEATH_DATE_COMBINED_USE = DEATH_DATE_COMBINED;
RETXDATE_USE = RETXDATE1;
if DEATH_DATE_COMBINED>'30SEP2011'd then DEATH_DATE_COMBINED_USE = .;
if RETXDATE1>'30SEP2011'd then RETXDATE_USE = .;
if DEATH_DATE_COMBINED_USE ne . then EVENT_DATE_USE=DEATH_DATE_COMBINED_USE; 
if DEATH_DATE_COMBINED_USE ne . then EVENT_INDICATOR=1; 
if RETXDATE_USE ne . then EVENT_DATE_USE=RETXDATE_USE ; 
if RETXDATE_USE ne . then EVENT_INDICATOR=1;
run;

proc contents data=public_nodup_final; run;

proc export data=public1 outfile='data/combined_all_listings.txt' dbms=tab replace; run;   
proc export data=public_nodup_final outfile='data/combined_patient_record.txt' dbms=tab replace; run;   

/*************************************************************************/
/* Get the longitudinal data for the patients actually included in the cohort
of interest. We will merge all the covariate information that is of interest later*/
/**************************************************************************/

data public_trans_dates_s;
set public1 (keep=PT_CODE WL_ID_CODE TRR_ID_CODE END_DATE_USE);
cohort=1; run; 

proc sort data=public_trans_dates_s;
by WL_ID_CODE; run;

proc sort data=LAS_audit;
by WL_ID_CODE; run;

data longitudinal; merge LAS_audit public_trans_dates_s;
by WL_ID_CODE; run;

data longitudinal;
set longitudinal;
where cohort=1 & datepart(begin_dt) le end_date_use & MATCH_LAS>0; run;

data longitudinal; set longitudinal;
drop cohort; run;

/*************************************************************************/
/* Create transplant data set and merge with the longitudinal dataset to create  
final dataset of use */
/**************************************************************************/

/* organ dataset of use */

data organ_dates;
set public1 (keep=WL_ID_CODE TRR_ID_CODE approx_match_date1 END_MATCH_LAS list_trans TX_PROCEDUR_TY AGE_DON ABO_DON HGT_CM_DON_CALC PO2 ECD_DONOR
GENDER_DON COD_CAD_DON HIST_CIG_DON CMV_DON RACE_DON DIABETES_DON GROUPING);
where list_trans=1; 
rename TRR_ID_CODE=T_TRR_ID_CODE; 
rename TX_PROCEDUR_TY=T_TX_PROCEDUR_TY; 
rename AGE_DON=T_AGE_DON;
rename ABO_DON=T_ABO_DON;
rename HGT_CM_DON_CALC=T_HGT_CM_DON_CALC;
rename PO2=T_PO2;
rename ECD_DONOR=T_ECD_DONOR; 
rename GENDER_DON = T_GENDER_DON;
rename COD_CAD_DON=T_COD_CAD_DON;
rename HIST_CIG_DON=T_HIST_CIG_DON;
rename CMV_DON=T_CMV_DON;
rename RACE_DON=T_RACE_DON;
rename DIABETES_DON = T_DIABETES_DON ;
run;

/* merge with longitudinal dataset - prune longitudinal dataset to one observation
per day using the following established rules - note that we will need to merge the organ dates data set first
Rules
1) Take the earliest observation except
2) If the LAS score corresponds to the match score (for organ allocation) then use that observation
*/

data organ_dates_a;
set organ_dates (keep=WL_ID_CODE END_MATCH_LAS); run;

proc sort data=organ_dates_a;
by WL_ID_CODE; run;

proc sort data=longitudinal; 
by WL_ID_CODE; run;

data longitudinal_s;
merge longitudinal organ_dates_a;
by WL_ID_CODE; 
LAS_diff = END_MATCH_LAS-MATCH_LAS; run;

data longitudinal_s;
set longitudinal_s;
LAS_match_exact=0;
if LAS_diff=0 then LAS_match_exact=1; 
MATCH_DATE = datepart(begin_dt); run;

proc sort data=longitudinal_s;
by WL_ID_CODE MATCH_DATE descending LAS_match_exact begin_dt; run;

data longitudinal_s;
set longitudinal_s; by WL_ID_CODE MATCH_DATE;
include=first.MATCH_DATE; run;

data longitudinal_s1;
set longitudinal_s;
where include=1; run;

data longitudinal_s1; set longitudinal_s1;
drop include LAS_diff LAS_match_exact END_MATCH_LAS; run;

proc sort data=longitudinal_s1;
by MATCH_DATE; run;

data longitudinal_s1_total;
set longitudinal_s1; run;

data organ_dates_total; set organ_dates; run;

data public1_total; set public1; run;

/**********************************************************/
%let group_var="A";
%let group_var1=A;
%let max_simultaneous=8;

%include "code_output/experiment_generation.sas";

proc univariate data=total_final1;
where TRANSPLANTED=1; var MATCH_LAS; 
output out=test pctlpre=P_ pctlpts= 5,35,65,95;run;
proc print data=test; run;
%let match_las_k1=31.0338;
%let match_las_k2= 32.8567;
%let match_las_k3=34.7963;
%let match_las_k4=47.7804;
proc univariate data=total_final1;
where TRANSPLANTED=1; var mov_sum; 
output out=test pctlpre=P_ pctlpts= 5,35,65,95;run;
proc print data=test; run;
%let mov_sum_k1=16;
%let mov_sum_k2= 52;
%let mov_sum_k3= 92;
%let mov_sum_k4=231;
proc univariate data=total_final1;
where TRANSPLANTED=1; var INIT_AGE_USE; 
output out=test pctlpre=P_ pctlpts= 5,35,65,95;run;
proc print data=test; run;
%let INIT_AGE_USE_k1=43;
%let INIT_AGE_USE_k2= 56;
%let INIT_AGE_USE_k3= 61.34;
%let INIT_AGE_USE_k4=68;
proc univariate data=total_final1;
where TRANSPLANTED=1; var HGT_DIFF; 
output out=test pctlpre=P_ pctlpts= 5,35,65,95;run;
proc print data=test; run;
%let HGT_DIFF_k1=-20;
%let HGT_DIFF_k2= -8;
%let HGT_DIFF_k3= -2.5;
%let HGT_DIFF_k4=7.88;
proc univariate data=total_final1;
where TRANSPLANTED=1; var LIST_DAYS; 
output out=test pctlpre=P_ pctlpts= 5,35,65,95;run;
proc print data=test; run;
%let LIST_DAYS_k1=5;
%let LIST_DAYS_k2=54;
%let LIST_DAYS_k3=181;
%let LIST_DAYS_k4=791;

%include "code_output/prob_generation_r2.sas";

/**********************************************************/
%let group_var="B";
%let group_var1=B;
%let max_simultaneous=3;

%include "code_output/experiment_generation.sas";

proc univariate data=total_final1;
where TRANSPLANTED=1; var MATCH_LAS; 
output out=test pctlpre=P_ pctlpts= 5,35,65,95;run;
proc print data=test; run;
%let match_las_k1=31.8114;
%let match_las_k2= 36.3870;
%let match_las_k3=43.9240;
%let match_las_k4=68.4761 ;
proc univariate data=total_final1;
where TRANSPLANTED=1; var mov_sum; 
output out=test pctlpre=P_ pctlpts= 5,35,65,95;run;
proc print data=test; run;
%let mov_sum_k1=16;
%let mov_sum_k2= 62;
%let mov_sum_k3= 109;
%let mov_sum_k4=244;
proc univariate data=total_final1;
where TRANSPLANTED=1; var INIT_AGE_USE; 
output out=test pctlpre=P_ pctlpts= 5,35,65,95;run;
proc print data=test; run;
%let INIT_AGE_USE_k1=23;
%let INIT_AGE_USE_k2= 43;
%let INIT_AGE_USE_k3= 55.29;
%let INIT_AGE_USE_k4=65;
proc univariate data=total_final1;
where TRANSPLANTED=1; var HGT_DIFF; 
output out=test pctlpre=P_ pctlpts= 5,35,65,95;run;
proc print data=test; run;
%let HGT_DIFF_k1=-14.26;
%let HGT_DIFF_k2= -4.7;
%let HGT_DIFF_k3= 2.54;
%let HGT_DIFF_k4= 12.7;
proc univariate data=total_final1;
where TRANSPLANTED=1; var LIST_DAYS; 
output out=test pctlpre=P_ pctlpts= 5,35,65,95;run;
proc print data=test; run;
%let LIST_DAYS_k1=4;
%let LIST_DAYS_k2=50;
%let LIST_DAYS_k3=148;
%let LIST_DAYS_k4=538;

%include "code_output/prob_generation_r2.sas";

/**********************************************************/
%let group_var="C";
%let group_var1=C;
%let max_simultaneous=4;

%include "code_output/experiment_generation.sas";

proc univariate data=total_final1;
where TRANSPLANTED=1; var MATCH_LAS; 
output out=test pctlpre=P_ pctlpts= 5,35,65,95;run;
proc print data=test; run;
%let match_las_k1=33.8619;
%let match_las_k2= 38.2099;
%let match_las_k3=43.5650;
%let match_las_k4=86.9715;
proc univariate data=total_final1;
where TRANSPLANTED=1; var mov_sum; 
output out=test pctlpre=P_ pctlpts= 5,35,65,95;run;
proc print data=test; run;
%let mov_sum_k1=14;
%let mov_sum_k2= 51;
%let mov_sum_k3= 91;
%let mov_sum_k4=224;
proc univariate data=total_final1;
where TRANSPLANTED=1; var INIT_AGE_USE; 
output out=test pctlpre=P_ pctlpts= 5,35,65,95;run;
proc print data=test; run;
%let INIT_AGE_USE_k1=19;
%let INIT_AGE_USE_k2= 25;
%let INIT_AGE_USE_k3= 33;
%let INIT_AGE_USE_k4=50;
proc univariate data=total_final1;
where TRANSPLANTED=1; var HGT_DIFF; 
output out=test pctlpre=P_ pctlpts= 5,35,65,95;run;
proc print data=test; run;
%let HGT_DIFF_k1=-19.98;
%let HGT_DIFF_k2= -7.62;
%let HGT_DIFF_k3= 0;
%let HGT_DIFF_k4=10.16;
proc univariate data=total_final1;
where TRANSPLANTED=1; var LIST_DAYS; 
output out=test pctlpre=P_ pctlpts= 5,35,65,95;run;
proc print data=test; run;
%let LIST_DAYS_k1=3;
%let LIST_DAYS_k2=42;
%let LIST_DAYS_k3=153;
%let LIST_DAYS_k4=767;

%include "code_output/prob_generation_r2.sas";

/**********************************************************/
%let group_var="D";
%let group_var1=D;
%let max_simultaneous=9;

%include "code_output/experiment_generation.sas";

proc univariate data=total_final1;
where TRANSPLANTED=1; var MATCH_LAS; 
output out=test pctlpre=P_ pctlpts= 5,35,65,95;run;
proc print data=test; run;
%let match_las_k1=34.5972;
%let match_las_k2= 41.1277;
%let match_las_k3=50.0875;
%let match_las_k4=90.4015 ;
proc univariate data=total_final1;
where TRANSPLANTED=1; var mov_sum; 
output out=test pctlpre=P_ pctlpts= 5,35,65,95;run;
proc print data=test; run;
%let mov_sum_k1=17;
%let mov_sum_k2= 58;
%let mov_sum_k3= 96;
%let mov_sum_k4=231;
proc univariate data=total_final1;
where TRANSPLANTED=1; var INIT_AGE_USE; 
output out=test pctlpre=P_ pctlpts= 5,35,65,95;run;
proc print data=test; run;
%let INIT_AGE_USE_k1=39.80;
%let INIT_AGE_USE_k2= 56;
%let INIT_AGE_USE_k3= 63;
%let INIT_AGE_USE_k4=70;
proc univariate data=total_final1;
where TRANSPLANTED=1; var HGT_DIFF; 
output out=test pctlpre=P_ pctlpts= 5,35,65,95;run;
proc print data=test; run;
%let HGT_DIFF_k1=-12.7;
%let HGT_DIFF_k2= -2.54;
%let HGT_DIFF_k3= 5.08;
%let HGT_DIFF_k4=15.24;
proc univariate data=total_final1;
where TRANSPLANTED=1; var LIST_DAYS; 
output out=test pctlpre=P_ pctlpts= 5,35,65,95;run;
proc print data=test; run;
%let LIST_DAYS_k1=2;
%let LIST_DAYS_k2=22;
%let LIST_DAYS_k3=82;
%let LIST_DAYS_k4=409;

%include "code_output/prob_generation_r2.sas";
