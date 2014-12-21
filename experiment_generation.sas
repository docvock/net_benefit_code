/*
#########################################################
### function to generate the patients in each experiment 
### this will be included in the main function
########################################################## */


/**************************************************************/
data longitudinal_s1;
set  longitudinal_s1_total; where GROUPING=&group_var; run;

data organ_dates; set organ_dates_total;
where GROUPING=&group_var; run;

data public1; set public1_total;
where GROUPING=&group_var; run;
/**************************************************************/

data organ_dates_a;
set organ_dates;
rename approx_match_date1=MATCH_DATE; 
drop LIST_TRANS; run;

proc sort data=organ_dates_a;
by match_date; run;

data organ_dates_loop;
set organ_dates_a (keep=MATCH_DATE T_TRR_ID_CODE T_TX_PROCEDUR_TY T_AGE_DON T_ABO_DON T_HGT_CM_DON_CALC T_PO2 T_ECD_DONOR
T_GENDER_DON T_COD_CAD_DON T_HIST_CIG_DON T_CMV_DON T_RACE_DON T_DIABETES_DON); run;

data organ_dates_loop;
set organ_dates_loop;
by match_date;
organ1 = first.match_date; run;


/**********************/
/* Round 1 */

data organ_dates1;
set organ_dates_loop; where organ1=1;
avail=1; run;

data test; merge organ_dates1 longitudinal_s1;
by MATCH_DATE; run;

data test; set test;
where avail=1; run;

data total; set test; run;

/**********************/
/* Round 2-&max_simultaneous */


%macro do_round;
%do i=1 %to &max_simultaneous;

data organ_dates_loop;
set organ_dates_loop; where organ1=0; run;

data organ_dates_loop;
set organ_dates_loop;
by MATCH_DATE; 
organ1 = first.match_date; run;

data organ_dates1;
set organ_dates_loop; where organ1=1;
avail=1; run;

data test; merge organ_dates1 longitudinal_s1;
by MATCH_DATE; run;

data test; set test;
where avail=1; run;

data total; set total test; run;
%end;
%mend;

%do_round;

/*******************************/

data total; set total;
drop organ1 avail FVC_PRE_LBV; 
TRANSPLANTED = 0;
if T_TRR_ID_CODE = TRR_ID_CODE then TRANSPLANTED=1; run;

proc sort data=total; by WL_ID_CODE; run;

proc sort data=public1;
 by WL_ID_CODE; run;

data public1; set public1;
rename TRANSPLANTED=EVER_TRANSPLANTED; run;

data total_final; merge total public1;
by WL_ID_CODE; run;

proc sort data=total_final;
by T_TRR_ID_CODE; run;

data total_final; set total_final;
where T_TRR_ID_CODE ne ""; run; 

data total_final1; set total_final;
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
if TX_DATE1=DEATH_DATE_COMBINED_USE & TX_DATE1 ne . then EVENT_DATE_USE=EVENT_DATE_USE+0.5;
RESID_LIFE = EVENT_DATE_USE-MATCH_DATE;
POT_CENSOR = '30SEP2011'd-match_date;
if TX_DATE1 ne . then do;
RESID_AFTER_TRANS=EVENT_DATE_USE-TX_DATE1;
RESID_TO_TRANS=TX_DATE1-MATCH_DATE;
end;
if TX_DATE1 = . then do;
RESID_AFTER_TRANS=0;
RESID_TO_TRANS=EVENT_DATE_USE-MATCH_DATE;
end;
DUMMY_TIME = 2-TRANSPLANTED;
run;

data total_final1; set total_final1;
cohort =1; run;
proc sort data=total_final1; by listing_ctr_code match_date; run;
data public_trans_tot_s; set public_trans_tot_s;
rename tx_date=match_date; run;

data total_final1; merge total_final1 public_trans_tot_s; 
by LISTING_CTR_CODE MATCH_DATE; run; 

data total_final1; set total_final1; where cohort=1; run;

proc print data=total_final1 (obs=10); run;

data total_final_add; set public_nodup_final (keep=PT_CODE INIT_AGE HGT_CM_TCR); 
rename INIT_AGE=LAST_AGE;
rename HGT_CM_TCR=FINAL_HGT; run;

proc sort data=total_final1; by PT_CODE; run;
proc sort data=total_final_add; by PT_CODE; run;
data total_final1; merge total_final1 total_final_add; by PT_CODE; run;
data total_final1; set total_final1; where cohort=1; run;
/*
proc export data=total_final1 outfile='C:\\Users\\bstvock\\Documents\\research\\net_benefit_CF\\combined_longitudinal_listing_GROUP_A.txt' dbms=tab replace; run;*/   

data total_final1;
set total_final1;
INIT_AGE_USE = ((init_date_las-init_date)/365+init_age); 
if HGT_CM_TCR <80 then HGT_CM_TCR=HGT_CM_TCR*2.54;
if INIT_AGE < 18 then HGT_CM_TCR=FINAL_HGT;
if LAST_AGE <15 & GENDER="M" then HGT_CM_TCR=172;
if LAST_AGE <15 & GENDER="F" then HGT_CM_TCR=162;
if HGT_CM_TCR=. & GENDER="M" then HGT_CM_TCR=172;
if HGT_CM_TCR=. & GENDER="F" then HGT_CM_TCR=162;
HGT_DIFF = HGT_CM_TCR-T_HGT_CM_DON_CALC;
LIST_DAYS = MATCH_DATE-INIT_DATE_LAS;
if T_PO2=. then T_PO2=425.45; 
if T_HIST_CIG_DON="Y" then T_HIST_CIG_DON_Y=1;
else T_HIST_CIG_DON_Y=0;
if T_CMV_DON="P" then T_CMV_DON_P=1;
else T_CMV_DON_P=0;
if T_DIABETES_DON="Y" then T_DIABETES_DON_Y=1;
else T_DIABETES_DON_Y=0;
if T_COD_CAD_DON=3 then T_COD_CAD_DON_HT=1;
else T_COD_CAD_DON_HT=0;
run;




