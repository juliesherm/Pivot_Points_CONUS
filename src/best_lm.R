## ---------------------------
##
## Script name: best_lm.R
##
## Purpose of script: returns the most favorable lm for each pixel by 
## comparing R squared over all VI and WB combinations
##
## Author: Julie Sherman
##
## Date Created: 2/13/2024
## Last edited: 2/14/2024
##
## ---------------------------
##
## Notes:
## - input argument is a list of rasters with different VI/WB combinations
## - each raster must have the layers "r.squared", "x.intercept", "y.intercept" and "slope"
##   such as returned by getlm
## - returns a raster with the same dimensions as one of the input rasters
##   and an additional layer indicating which combination had the lowest R2 and
##   thus provided the data for the returned fit
##
## ---------------------------

best_lm<-function(rasts){
  rsqs<-rast(lapply(rasts, subset, "r.squared"))
  best_combo<-which.max(rsqs) #returns raster of indices of highest rsquared for each pixel
  #create blank raster of the same resolution
  best_rast<-rasts[[1]]  
  values(best_rast)<-NA
  #fill in best_raster with values from the model with highest r squared
  for(i in 1:length(rasts)){
    best_rast$slope[best_combo==i]<-rasts[[i]]$slope[best_combo == i]
    best_rast$y.intercept[best_combo==i]<-rasts[[i]]$y.intercept[best_combo == i]
    best_rast$x.intercept[best_combo==i]<-rasts[[i]]$x.intercept[best_combo == i]
    best_rast$r.squared[best_combo==i]<-rasts[[i]]$r.squared[best_combo == i]
  }
  #append a layer indicating which variable the data at each cell derives from
  best_rast<-c(best_rast,best_combo)
  return(best_rast)
}