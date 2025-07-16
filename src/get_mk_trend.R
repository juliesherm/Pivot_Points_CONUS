## ---------------------------
##
## Script name: get_mk_trend.R
##
## Purpose of script: gets Mann-Kendall tau and p-value for years with complete data 
## from vector of yearly dates
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

get_mk_trend <- function(vec) {
  
  if (any(is.na(vec))) {
    
    # fill NAs with linear interpolation (mean of surrounding two data points)
    which_missing <- which(is.na(vec))
    vec[which_missing] <- sapply(which_missing, function(i) {
      (vec[i-1] + vec[i+1])/2
    })
    
    # if there is still a missing number 
    if (any(is.na(vec))) {
      return(c(NA,NA))
    } 
  } 

  my_mk <- MannKendall(unlist(vec))
    
  return(c(my_mk$tau, my_mk$sl))

}
