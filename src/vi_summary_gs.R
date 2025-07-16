vi_summary_gs <- function(vi_stk, vi_dates, gs_dates, var) {

  print(paste0("Calculating growing season summaries for ", var))
  
  #for each year, apply provided summary function to the day of years in gs_dates
  #save the summary layer
  yrs <- year(vi_dates)
  
  
  vi_gs <- lapply(unique(yrs), function(yr) {
    
    start <- Sys.time()
    
    if(paste0("first",yr) %in% names(gs_dates)) {
      gs_rast <- c(gs_dates[paste0("first",yr)], gs_dates[paste0("last",yr)], vi_stk[[yr == yrs]])
      
    } else {
      gs_rast <- c(gs_dates["medianfirst"], gs_dates["medianlast"], vi_stk[[yr == yrs]])
    }
    
    vi_yday <- yday(vi_dates[yr == yrs])
    
    out <- terra::app(gs_rast, 
                      fun = integrate_gs, 
                      ydays = vi_yday)
    end <- Sys.time()
    print(paste0("calculating ", yr, " took ", difftime(end, start, units="secs"), " seconds"))
    
    return(out)
  })
  vi_gs <- rast(vi_gs)
  names(vi_gs) <- unique(yrs)
  
  return(vi_gs)
}
