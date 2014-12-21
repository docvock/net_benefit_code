/*
#######################################################
### SAS code to generate the tables needed for consort diagram
###
### PROGRAMS RELIED ON: none
### SAS DATA FILES USED: thoracic_public_use_data 
###                      thoracic_las_audit_data
###                      thoracic_las_history_data
### ADDITIONAL PROGRAMS USED: none
### SAS DATA FILES GENERATED: none
### TXT DATA FILES GENERATED: none
#######################################################
 */

options formdlim="_" nodate nonumber;

Libname newlib 'C:\\Users\\bstvock\\Documents\\UNOS data\\updated_data_12_31_2011';  

/*************************************************************/
/* read-in data */
/*************************************************************/
data public;
set newlib.thoracic_public_use_data (keep=
PT_CODE WL_ID_CODE TRR_ID_CODE /*patient identifiers*/
citizenship init_age ABO BMI_TCR PRI_PAYMENT_TCR ETHCAT RACE REGION GENDER HGT_CM_TCR TCR_DGN TCR_DGN_OSTXT /*baseline patient characteristics*/
init_date init_stat organ WL_ORG WLHL WLHR WLIN WLKI WLKP WLLI WLPA WLPI listing_ctr_code init_blu_flg init_llu_flg init_rlu_flg/*Listing information*/
end_match_las TX_DATE TXED TXHRT TXINT TXKID TXLIV TXPAN TX_PROCEDUR_TY ABO_DON ABO_MAT AGE AGE_DON ECD_DONOR PO2 COD_CAD_DON 
CONTIN_CIG_DON HIST_CIG_DON ISCHTIME SHARE_TY HGT_CM_DON_CALC GENDER_DON CMV_DON RACE_DON DIABETES_DON/*transplant and donor characteristics*/
PX_STAT PX_STAT_DATE  death_date end_date REM_CD RETXDATE SSDMF_DEATH_DATE /*status information*/    );run;

data LAS_audit; 
set newlib.thoracic_las_audit_data; run;

data LAS_hist;
set newlib.thoracic_las_history_data; run;

/************************************************************/
/* Combine LAS grouping covariate with audit data */
/************************************************************/
proc sort data=LAS_hist; by WL_ID_CODE; run;

data LAS_hist; set LAS_hist (keep=WL_ID_CODE GROUPING); 
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
first_patient = first.PT_CODE; run;

data public_death_date_s; set public_death_date;
where first_patient=1; run;

data public_death_date_s; set public_death_date_s;
DEATH_DATE_SOURCE = 0;
DEATH_DATE_COMBINED = death_date;
if death_date ne . then death_date_source=1;
if death_date= . then death_date_combined=PX_death_date;
if death_date= . & PX_death_date ne . then death_date_source=2;
if death_date= . & PX_death_date=. then death_date_combined=SSDMF_death_date;
if death_date= . & PX_death_date=. & SSDMF_death_date ne . then death_date_source=3;run;

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

/* Initial Inclusion Criteria */
data public1; set public;
where init_stat=7010 or init_stat=7999; 
LIST_TRANS = 0;
if TRR_ID_CODE ne "" & TX_DATE le '30SEP2011'd then LIST_TRANS=1; run;

data public1; set public1;
END_DATE_ALT =END_DATE;
if init_date le tx_date1 then END_DATE_ALT=TX_DATE1; run;
 
data public1; set public1;
where end_date_alt ge '4MAY2005'd & init_date le '30SEP2011'd ; run;


data public1; set public1;
TRANSPLANTED = 0;
if tx_date1 ne .  & tx_date1 ge '4MAY2005'd then TRANSPLANTED=1; 
TRANSPLANTED1 = 1;
if tx_date1 ne .  then TRANSPLANTED1=1;
run;

data public_nodup;
set public1;  run;
/*
proc freq data=public_nodup;
tables LIST_TRANS TRANSPLANTED; run; */

proc sort data=public_nodup;
by PT_CODE INIT_DATE; run;

data public_nodup; set public_nodup;
end_date_use = min(end_date,tx_date1,death_date_combined,'30SEP2011'd); 
LAS_AVAIL = 0;
if INIT_DATE_LAS ne . then LAS_AVAIL=1;
AFTER_TRANS=0;
AFTER_TRANS1=0;
if TX_DATE1 ne . & TX_DATE1 ge '4MAY2005'd & INIT_DATE > TX_DATE1 then AFTER_TRANS=1;
if TX_DATE1 ne . & TX_DATE1 ge '4MAY2005'd & INIT_DATE = TX_DATE1 & TX_DATE > TX_DATE1 then AFTER_TRANS=1;
if TX_DATE1 ne . & TX_DATE1 ge '4MAY2005'd & INIT_DATE = TX_DATE1 & TX_DATE =. then AFTER_TRANS=1;
if TX_DATE1 ne . & INIT_DATE > TX_DATE1 then AFTER_TRANS1=1;
if TX_DATE1 ne . & INIT_DATE = TX_DATE1 & TX_DATE > TX_DATE1 then AFTER_TRANS1=1;
if TX_DATE1 ne . & INIT_DATE = TX_DATE1 & TX_DATE =. then AFTER_TRANS1=1;
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
end;
run;

proc freq data=public_nodup;
tables AFTER_TRANS*AFTER_TRANS1; run;
/*
proc print data=public_nodup;
where AFTER_TRANS=0 & AFTER_TRANS1=1; run;
proc print data=public_nodup; where PT_CODE=7521; run; */

proc sort data=public_nodup; 
by PT_CODE descending LAS_AVAIL AFTER_TRANS descending end_date_use descending REM_TRANS descending REM_DEATH descending REM_DETER descending REM_HEALTHY descending REM_OTHER descending END_MATCH_LAS descending INIT_DATE; run;

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

data public_nodup1; set public_nodup1 (keep=PT_CODE INIT_DATE_LAS INIT_MATCH_LAS INIT_DATE INIT_AGE);
rename INIT_DATE = INIT_DATE1;
rename INIT_DATE_LAS=INIT_DATE_LAS1;
rename INIT_MATCH_LAS=INIT_MATCH_LAS1; 
INIT_AGE_LAS = INIT_AGE+(INIT_DATE_LAS-INIT_DATE)/365;
run;

data public_nodup1; set public_nodup1;
drop INIT_AGE;

data public_nodup_final; merge public_nodup_final public_nodup1; 
by PT_CODE; run;

data public_nodup_final; set public_nodup_final;
FAILURE = 0;
if DEATH_DATE_COMBINED ne . & DEATH_DATE_COMBINED le '30SEP2011'd & DEATH_DATE_COMBINED ge '01JAN2001'd then FAILURE=1; run;

data public1; set public_nodup_final;
run;

/*************************************************************/
/* restrict analysis using the inclusion/exclusion criteria described */
/*************************************************************/
title "Initial Inclusion Criteria"; run;
proc freq data=public1; tables TRANSPLANTED*FAILURE; run;

title "Exclude Prior Transplants"; run;
data public1; set public1;
PRIOR_TX = 0;
if TX_DATE1 ne . & INIT_DATE > TX_DATE1 then PRIOR_TX=1; run;

data public1; set public1;
where PRIOR_TX=0; run;

proc freq data=public1; tables TRANSPLANTED*FAILURE; run;

title "Exclude Listed for Multi-Organ Transplants"; run;
data public1; set public1;
where WLHR="" & WLHL="" & WLIN="" & WLKI="" & WLKP="" & WLLI="" & WLPA="" & WLPI="" &
PT_CODE ne 125779 & PT_CODE ne 619743 & PT_CODE ne 632410 & PT_CODE ne 638611 & PT_CODE ne 644831 &
PT_CODE ne 663088 & PT_CODE ne 665398 & PT_CODE ne 734998 & PT_CODE ne 750685 & PT_CODE ne 778780 & 
PT_CODE ne 792102 & PT_CODE ne 871377; run;
/* PT_CODE 125779 listed initially for lung only and then lung+liver subsequently */

proc freq data=public1; tables TRANSPLANTED*FAILURE; run;

title "Exclude Only Zero LAS Values"; run;

data public1; set public1;
where INIT_DATE_LAS ne . & INIT_DATE_LAS le '30SEP2011'd;; run;

data public1; set public1;
where TRR_ID_CODE = "" or END_MATCH_LAS > 0; run;

proc freq data=public1; tables TRANSPLANTED*FAILURE; run;

title "Known to be Over 18 years of Age"; run;

data public1; set public1;
where init_age ge 18 or ((init_date_las-init_date)/365+init_age ge 18); run;

proc freq data=public1; tables TRANSPLANTED*FAILURE; run;

title "Exclude Living Donors and Single LTx"; run;

data public1; set public1;
nonliving = 1;
if TRR_ID_CODE ne "" & TX_DATE le '30SEP2011'd & TXED=0 then nonliving=0; 
if REM_CD=15 & END_DATE le '30SEP2011'd then nonliving=0; run;

data public1; set public1;
where nonliving=1; run;

proc freq data=public1; tables TRANSPLANTED*FAILURE; run;

title "Exclude Non-citizens"; run;

data public1; set public1;
where citizenship=1 & PT_CODE ne 635399 & PT_CODE ne 778539 & PT_CODE ne 862252; run;
/* PT_CODE=635399 was listed once as citizen and then once as permanent resident*/

data public1;
set public1; 
if GROUPING="" then GROUPING="D"; run; 
/* PT_CODE=573640 has missing grouping but TCR_DGN=422 which is always assigned group D */

proc freq data=public1; tables TRANSPLANTED*FAILURE; run;

proc freq data=public1; tables TRANSPLANTED*FAILURE*GROUPING; run;

proc freq data=public1; tables FAILURE*GROUPING; run;

