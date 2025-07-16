###############################################################
# plot one of the outputs of pivot point analysis for all of CONUS
# save to pdf
# pivot var
library(here)
library(terra)
library(dplyr)
piv_var = "slope"

VIS = c("SAVI")
WBS = c("accumswe","aet","deficit","pet",
        "rain","runoff","soilwater")#,"precip","tmin","tmax","tavg")
sig = TRUE

pdf(here("results","CONUS","figures",
         paste0(piv_var,"_sig_SAVI_pivot_plots.pdf")))
par(mfrow = c(3,3))
for(vi in VIS){
  for(wb in WBS){
    print(wb)
    gc()
    pivs<-terra::rast(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear",
                           paste0(vi,"_",wb,"_",piv_var,".tif")))
    if(sig){
      print("Masking non-significant p-values")
      pvals<-terra::rast(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear",
                              paste0(vi,"_",wb,"_pval.tif")))
      pivs<-mask(pivs,pvals>0.1, maskvalues = TRUE, updatevalue = NA)      
      remove(pvals)
    }
    
    print("Calculating bounds")
    if(piv_var=="rsq" || piv_var=="pval"){
      ub=1
      lb=0
    }else{
      ub = global(pivs, fun=quantile, probs = 0.99, na.rm = TRUE)$X99
      lb = global(pivs, fun=quantile, probs = 0.01, na.rm = TRUE)$X1
      if(is.infinite(ub)){ #soilwater slope has inf for 95 percentile as well
        ub = global(pivs, fun=quantile, probs = 0.90, na.rm = TRUE)$X90
      }
      if(is.infinite(lb)){
        lb = global(pivs, fun=quantile, probs = 0.10, na.rm = TRUE)$X10
      }
      if(piv_var=="slope"){
        ub = max(abs(ub),abs(lb))
        lb = -ub
        if(is.infinite(lb)){
          lb = global(pivs, fun=quantile, probs = 0.81, na.rm = TRUE)$X81
          ub = min(abs(ub),abs(lb))
          lb = -ub
        }
        pal<-leaflet::colorNumeric(palette="RdBu", domain = c(lb,ub), reverse = T)
      }
    }
    print("Plotting")
    plot(pivs,1,range = c(lb,ub), main = paste(vi,wb), 
         col = pal(seq(lb,ub,length.out = 100)),
         colNA = "lightgrey")
    remove(pivs)
  }
  # blank plot so that EVI and SAVI show up "evenly"
  #pivs<-terra::rast(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear",
  #                       paste0(vi,"_",wb,"_",piv_var,".tif")))
  #plot(pivs,1, range = c(0,0.01))
  #remove(pivs)
}
dev.off()

###############################################################
# plot where relationships are significant side-by-side
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


###############################################################
# plot where relationships are significant in 1 plot both vis
library(here)
library(terra)
library(dplyr)
library(tidyterra)
piv_var = "pval"

WBS = c("accumswe","aet","deficit","pet",
        "rain","runoff","soilwater",
        "precip","tmin","tmax","tavg")
pas<-vect(here("raw_data", "PADUS4_0VectorAnalysis_State_ClipCensus_CONUS.gdb",
               "a00000028.gdbtable"))
#pas<-vect(here("raw_data", "PADUS4_0Geodatabase","PADUS4_0_Geodatabase.gdb","a00000018.gdbtable")) #contains all US (not just CONUS)

#Unit "Non-PAD-US Area"
#Area of Critical Environmental Concern"
#Oklahoma largely tribal lands?

pas<-pas %>% filter(FeatClass != "Marine" & Unit_Nm != "Non-PAD-US Area" )
pas<-pas %>% filter(GAP_Sts < 3 )

pas<-project(pas,"epsg:4326")



pdf(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear",
         paste0(piv_var,"_sig_pivot_plots_PAD.pdf")))
par(mfrow = c(4,3))
for(wb in WBS){
  print(wb)
  gc()
  pivs<-terra::rast(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear",
                         paste0("SAVI","_",wb,"_",piv_var,".tif")))
  pivs2<-terra::rast(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear",
                          paste0("EVI","_",wb,"_",piv_var,".tif")))
  pivs<-mask(pivs, pas)
  pivs2<-mask(pivs2, pas)
  
  print("Calculating bounds")
  pivs<-pivs<0.1
  pivs2<-pivs2<0.1
  comb <-1*(pivs & !pivs2)+2*(!pivs & pivs2)+3*(pivs & pivs2)
  sav = global(pivs & !pivs2, "mean", na.rm=TRUE)
  ev = global(!pivs & pivs2, "mean", na.rm=TRUE)  
  both = global(pivs & pivs2, "mean", na.rm=TRUE)
  non = 1-ev-sav-both
  remove(pivs)
  remove(pivs2)
  print("Plotting")
  plot(comb,1, main = paste("p<0.1",wb),
       #plg=list(legend=c("None", "SAVI","EVI","Both")),
       plg=list(legend=c(paste("None\n", round(non,3)), paste("SAVI\n",round(sav,3)),paste("EVI\n",round(ev,3)),paste("Both\n",round(both,3))),
                cex = 0.8),
       col = c("grey81","orange2","cornflowerblue","forestgreen"))
  remove(comb)
}
dev.off()


###############################################################
# just  protected areas
library(here)
library(terra)
library(dplyr)
library(tidyterra)
piv_var = "pval"

pivs<-rast(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear",
                "SAVI_precip_xint.tif"))
pas<-vect(here("raw_data", "PADUS4_0VectorAnalysis_State_ClipCensus_CONUS.gdb",
               "a00000028.gdbtable"))
#pas<-vect(here("raw_data", "PADUS4_0Geodatabase","PADUS4_0_Geodatabase.gdb","a00000018.gdbtable")) #contains all US (not just CONUS)

names(pas)
#Unit "Non-PAD-US Area"
#Area of Critical Environmental Concern"
#Oklahoma largely tribal lands?

pas<-pas %>% filter(FeatClass != "Marine" & Unit_Nm != "Non-PAD-US Area" )
paslow<-pas %>% filter(GAP_Sts < 3 )

pas<-project(pas,"epsg:4326")

pivs<-mask(pivs, pas)
pivslow<-mask(pivs, paslow)
plot(pas)

plot(pivs, range = c(-10,500))
plot(pivslow, range = c(-10,500))

polys(pas)
plot(pas, add = TRUE)

###########################

WBS = c("accumswe","aet","deficit","pet",
        "rain","runoff","soilwater",
        "precip","tmin","tmax","tavg")


pdf(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear",
         paste0(piv_var,"_sig_pivot_plots_vis.pdf")))
par(mfrow = c(4,3))
for(wb in WBS){
  print(wb)
  gc()
  pivs<-terra::rast(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear",
                         paste0("SAVI","_",wb,"_",piv_var,".tif")))
  pivs2<-terra::rast(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear",
                          paste0("EVI","_",wb,"_",piv_var,".tif")))
  print("Calculating bounds")
  pivs<-pivs<0.1
  pivs2<-pivs2<0.1
  comb <-1*(pivs & !pivs2)+2*(!pivs & pivs2)+3*(pivs & pivs2)
  sav = global(pivs & !pivs2, "mean", na.rm=TRUE)
  ev = global(!pivs & pivs2, "mean", na.rm=TRUE)  
  both = global(pivs & pivs2, "mean", na.rm=TRUE)
  non = 1-ev-sav-both
  remove(pivs)
  remove(pivs2)
  print("Plotting")
  plot(comb,1, main = paste("p<0.1",wb),
       #plg=list(legend=c("None", "SAVI","EVI","Both")),
       plg=list(legend=c(paste("None\n", round(non,3)), paste("SAVI\n",round(sav,3)),paste("EVI\n",round(ev,3)),paste("Both\n",round(both,3))),
                cex = 0.8),
       col = c("grey81","orange2","cornflowerblue","honeydew2"))
  remove(comb)
}
dev.off()

