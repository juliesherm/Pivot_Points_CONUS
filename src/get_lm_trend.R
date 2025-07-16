## ---------------------------
##
## Script name: get_lm_trend.R
##
## Purpose of script: gets y-intercept, slope, r-squared for linear model fit to 
## vector of yearly dates
##
## Author: Carolyn Lober
##
## Date Created:  2/14/2024
## Last edited: 2/14/2024
##
## ---------------------------
##
## Notes:
## - data should be yearly
## - could probably be combined with getlm? 
##
## ---------------------------

get_lm_trend <- function(vec,yrs) {
  if (any(is.na(vec))) {
    as.numeric(c(NA, NA, NA))
  } else {
    my_lm <- lm(vec ~ yrs)
    # return y-intercept, slope, r squared
    return(as.numeric(c(
      my_lm$coefficients, summary(my_lm)$r.squared
    )))
  }
}
