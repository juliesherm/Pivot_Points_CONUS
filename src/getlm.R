## ---------------------------
##
## Script name: getlm.R
##
## Purpose of script: From a vector, fits a linear model using formula 1st half ~ 2nd half
## and returns y.intercept, slope, r.squared, and x.intercept
##
## Author: Carolyn Lober, Julie Sherman
##
## Date Created: 2/9/2024
## Last edited: 3/14/2024
##
## ---------------------------
##
## Notes:
##  - vector MUST be divisible by 2
##  - first half of vector should be vegetation index (dependent variable),
##    and second half should be water balance (independent variable)
##
## ---------------------------

getlm <- function(vi_wb_vec) {
  len <- length(vi_wb_vec) / 2

  if (sum(is.na(vi_wb_vec)) > 2) {
    out <- as.numeric(rep(NA, times=8))
  } else {
    my_lm <- lm(vi_wb_vec[1:len] ~ vi_wb_vec[(1 + len):(len * 2)])
    
    # return y-intercept, slope, r squared, x-intercept, also return range(x)
    out <- as.numeric(c(
      my_lm$coefficients, # intercept and slope
      summary(my_lm)$r.squared, # r-squared
      (0 - my_lm$coefficients[1]) / my_lm$coefficients[2], # x-intercept
      summary(my_lm)$coefficients[,4], # p-value
      range(vi_wb_vec[(1 + len):(len * 2)], na.rm = TRUE)
    ))
    
    
  }
  
  #if (typeof(out) != "double" | length(out) != 8) {
  #  print(vi_wb_vec)
  #  stop("Something didn't calculate correctly")
  #}
  
  return(out)
}
