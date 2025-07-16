## ---------------------------
##
## Script name: smooth_ts.R
##
## Purpose of script: From a vector (i.e. raster stack values for one pixel), 
## set winter values to 60% of the winter (Jan) max and apply a Savitzky-Golay 
## smoothing filter
##
## Author: Carolyn Lober
##
## Date Created: 2/15/2024
## Last edited: 3/15/2024
##
## ---------------------------
##
## Notes:
##  - requires packages 'signal' and 'xts'
##
## ---------------------------

smooth_ts <- function(vec, dates) {
  
  if (length(vec) != length(dates)) {
    stop("arguments 'vec' and 'dates' must be the same length")
  }
  
  if (!is.Date(dates)) {
    stop("'dates' must be a date object")
  }
  
  # convert to time series object
  vec_ts <- xts(vec, order.by = dates)
  
  # set anything below 60% of winter max to that value (minimum realistic value per Wang et al.)
  min_winter_val <- max(vec_ts[month(vec_ts) == c(12,1,2)], na.rm = TRUE) * 0.6
  #max_winter_val <- # to filter out winter spikes
  if (is.na(min_winter_val)) {
    stop("no min winter value")
  }
  
  vec_ts[which(vec_ts < min_winter_val)] <- min_winter_val
  
  # pad beginning with numbers before na.approx with na.approx
  vec_ts <- na.fill(vec_ts, fill = c("extend","extend","extend"))
  
  # fill all gaps with linear interpolation, longer gaps with min winter value 
  #vec_ts <- na.approx(vec_ts, na.rm = FALSE) 
  
  smoothed_vec <- sgolayfilt(vec_ts, p = 1 ,n = 5, m = 0, ts = 3) %>%
    xts(.,order.by = index(vec_ts)) 
  
  #smoothed_vec_2 <-  sgolayfilt(vec_ts, p = 1 ,n = 5, m = 0, ts = 3) %>%
  #  xts(.,order.by = index(vec_ts)) 
  #smoothed_vec_2 <- approx(y = as.numeric(smoothed_vec_2), x = index(smoothed_vec_2),
  #                       xout = seq.Date(from=min(index(smoothed_vec_2)), to=max(index(smoothed_vec_2)), by="day"))
  #smoothed_vec_2 <- xts(smoothed_vec_2$y, order.by = smoothed_vec_2$x)
  
  #smoothed_vec_2001 <- sgolayfilt(vec_ts["2020"], p = 1 ,n = 5, m = 0, ts = 3) %>%
  #  xts(.,order.by = index(vec_ts["2020"]))
  
  #plot(vec_ts)
  #lines(smoothed_vec,col="maroon")
  
  # return yearly summaries instead?
  return(as.numeric(smoothed_vec))
}
