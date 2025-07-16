## ---------------------------
##
## Script name: stk_trend.R
##
## Purpose of script: gets y-intercept, slope, r-squared for linear model, 
## and Mann-Kendall tau and p-value for years with complete data from raster stack 
## and dates index vector
##
## Author: Carolyn Lober
##
## Date Created:  2/14/2024
## Last edited: 2/14/2024
##
## ---------------------------
##
## Notes:
## - if summary_choice is "gs", stk should be pre-summarized growing season integrated VI 
## and in this case dates is not used
##
## ---------------------------

stk_trend <- function(stk, dates, summary_choice) {
  
  if (class(stk) != "SpatRaster") {
    stop("argument 'stk' is not valid. must be a SpatRaster")
  }
  
  if (class(dates) != "Date" & summary_choice != "gs") {
    stop("argument 'dates' is not valid. must be a Date")
  }
  
  if (dim(stk)[3] != length(dates) & summary_choice != "gs") {
    stop("length of dates does not match number of layers.")
  }
  
  if (summary_choice %in% c("mean", "max")) {
    full_yrs <- as.data.frame(dates) %>% count(year(dates)) %>% dplyr::filter(n > 10) %>% .$`year(dates)`
    full_yrs_stk <- subset(stk, year(dates) %in% full_yrs)
    full_yrs_dates <- subset(dates, year(dates) %in% full_yrs)
  }
    
  if (summary_choice == "mean") {
    annual_summary <- tapp(full_yrs_stk, index=year(full_yrs_dates), fun = mean, na.rm = TRUE)
  } else if (summary_choice == "max") {
    annual_summary <- tapp(full_yrs_stk, index=year(full_yrs_dates), fun = max, na.rm = TRUE)
  } else if (summary_choice == "gs") {
    annual_summary <- stk
    full_yrs <- as.numeric(names(stk))
  } else {
    stop("only functions 'mean' and 'max' are defined, please choose one of those")
  }
    
  out_lm <- terra::app(annual_summary, get_lm_trend, yrs=full_yrs)
  names(out_lm) <- c("y.intercept","slope","r.squared")
  out_mk <- terra::app(annual_summary, get_mk_trend)
  names(out_mk) <- c("tau","p.value")

  return(c(out_lm,out_mk))
}
