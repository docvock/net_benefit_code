## Overall code which runs all necessary R files to give outout for net benefit analysis  
## Assumes that all code is in current working directory


## code to generate tables for characteristics of recipients, donors, surgery, and transplant type
source("baseline_demographics.R")
rm(list = ls())


## code to generate pdfs for waitlist dynamics and post-LTx survival
source("WL_post_LTx_survival_graphics.R")
rm(list = ls())


## code to estimate structural nested accelerated failure time models
start_overall <- proc.time()

## read-in functions and libraries used
source("library_func.R")

## read-in data
delta.use <- 5
c.adj <- 90
source("data_read_in.R")


## obtain starting value estimates of the baseline residual survival and the parameters in the scale 
## accelerated failure time model
source("gamma_generate.R")


## obtain estimates of the parameters in the structural nested Scale Accelerated Failure Time Model
## (SAFTM)
eps2 <- 5
source("code_output/coef_est_SAFTM_streamlined_r2.R")

###################
### obtain estimates of the parameters in the structural nested non-scale accelerated failure time model 
###################
eps <- 5
eps2 <- 5

source("code_output/coef_est_NSAFTM_streamlined_r2_final.R")

proc.time()-start_overall
###################
### output tables and graphics of model 
###################

source("code_output/model_output_graphics_streamlined_r2.R")

###################
### model diagnostics 
###################

source("code_output/internal_validation_streamlined_r2.R")


proc.time()-start_overall
