###############################################################
# pivot var
library(here)
library(terra)
library(dplyr)
piv_var = "rsq"

VIS = c("EVI",
  "SAVI")
WBS = c("accumswe","aet","deficit","pet",
        "rain","runoff","soilwater",
        "precip","tmin","tmax","tavg")


pdf(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear",
         paste0(piv_var,"_all_pivot_plots.pdf")))
par(mfrow = c(6,4))
for(vi in VIS){
  for(wb in WBS){
    print(wb)
    gc()
    pivs<-terra::rast(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear",
                           paste0(vi,"_",wb,"_",piv_var,".tif")))
    print("Calculating bounds")
    ub = global(pivs, fun=quantile, probs = 0.99, na.rm = TRUE)$X99
    lb = global(pivs, fun=quantile, probs = 0.01, na.rm = TRUE)$X1
    if(piv_var=="rsq" || piv_var=="pval"){
      ub=1
      lb=0
    }
    print("Plotting")
    plot(pivs,1,range = c(lb,ub), main = paste(vi,wb))
    remove(pivs)
  }
  # blank plot so that EVI and SAVI show up "evenly"
  pivs<-terra::rast(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear",
                         paste0(vi,"_",wb,"_",piv_var,".tif")))
  plot(pivs,1, range = c(0,0.01))
  remove(pivs)
}
dev.off()

###############################################################
# plot where relationships are significant
library(here)
library(terra)
library(dplyr)
piv_var = "pval"

VIS = c("EVI",
        "SAVI")
WBS = c("accumswe","aet","deficit","pet",
        "rain","runoff","soilwater",
        "precip","tmin","tmax","tavg")


pdf(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear",
         paste0(piv_var,"_sig_pivot_plots.pdf")))
par(mfrow = c(6,2))
for(wb in WBS){
  for(vi in VIS){
    print(wb)
    gc()
    pivs<-terra::rast(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear",
                           paste0("SAVI","_",wb,"_",piv_var,".tif")))
    pivs2<-terra::rast(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear",
                           paste0("EVI","_",wb,"_",piv_var,".tif")))
    print("Calculating bounds")
    pivs<-pivs<0.1
    print("Plotting")
    plot(pivs,1, main = paste("p<0.1",vi,wb))
    remove(pivs)
  }
}
dev.off()
