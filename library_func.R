##########################################################
### libraries used in the analysis
##########################################################

library(xtable)
library(survival)
library(data.table)

######################################
### functions used during the optimization
######################################

smooth.indicator <- function(x,c,eps) {
smooth.indicator <- ifelse(x<(c-eps),1,
	ifelse(x>=(c-eps)&x<(c-eps/2),-2/eps^2*(x-(c-eps))^2+1,
	ifelse(x>=(c-eps/2)&x<c, 2/eps^2*(x-c)^2,0)))
return(smooth.indicator) }
###
smooth.indicator.deriv.x <- function(x,c,eps) {
smooth.indicator.deriv.x <- ifelse(x<(c-eps),0,
	ifelse(x>=(c-eps)&x<(c-eps/2),-4/eps^2*(x-(c-eps)),
	ifelse(x>=(c-eps/2)&x<c, 4/eps^2*(x-c),0)))
return(smooth.indicator.deriv.x) }
###
smooth.indicator.deriv.c <- function(x,c,eps) {
smooth.indicator.deriv.c <- ifelse(x<(c-eps),0,
	ifelse(x>=(c-eps)&x<(c-eps/2),4/eps^2*(x-(c-eps)),
	ifelse(x>=(c-eps/2)&x<c, -4/eps^2*(x-c),0)))
return(smooth.indicator.deriv.c) }
###
smooth.min <- function(x,c,eps) {
smooth.min <- ifelse(x-c<0,x-c,0)*smooth.indicator(x,c,eps)+c
return(smooth.min) }
###
smooth.min.deriv.x <- function(x,c,eps) {
smooth.min.deriv.x <- ifelse(x<c,smooth.indicator(x,c,eps)+smooth.indicator.deriv.x(x,c,eps)*(x-c),0)
return(smooth.min.deriv.x) }
###
smooth.min.deriv.c <- function(x,c,eps) {
smooth.min.deriv.c <- ifelse(x<c,-smooth.indicator(x,c,eps)+smooth.indicator.deriv.c(x,c,eps)*(x-c)+1,1)
return(smooth.min.deriv.c) }