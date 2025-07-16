###############################################################
# plot one of the outputs of pivot point analysis for all of CONUS
# save to pdf
# pivot var
setwd("C:/Users/Julie/Documents/MATH/NPS SIP/pivot_pts")
setwd("/home/1024/ma/sherman/NPS_SIP/pivot_pts")

library(here)
library(terra)
library(dplyr)
library(ggplot2)
library(tidyterra)
library(scales)
source(here("src","get_ext.R"))



piv_out_plot<-function(wb = "deficit",
                      vi = "EVI",
                      outvar = "rsq", 
                      park_extent = NULL,
                      filter_nlcd = TRUE){
  print(outvar)
  gc()
  if(is.null(park_extent)){
      varrast<-terra::rast(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear",
                               paste0(vi,"_",wb,"_",outvar,".tif")),
                           win = park_extent)
      states<-terra::vect(here("raw_data","states","s_18mr25.shp"))
      states<-crop(states,  varrast)
      
  }else{
    varrast<-terra::rast(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear","tiled",vi,
                              paste0("pivs_",vi,"_",wb,"_0_5222.tif")),
                         win = park_extent)
    lyrnm<-switch(outvar,
      pval="pvalue",
      xint="x-intercept",
      slope="slope",
      rsq="rsqaured"
    )
    varrast<-subset(varrast, lyrnm)
  }
  
  # mask to significant p-vals
  # if filtering out developed...
  fname<-ifelse(filter_nlcd,"sig_pvals_masked.tif","sig_pvals.tif")
  sig_mask<-terra::rast(here("results","CONUS","pivot_pts",fname),
                        win = park_extent) 
  sig_mask<-subset(sig_mask,paste0(vi,"_",wb,"_pval"))
  
  tit<-switch(outvar,
              pval = "p-value",
              xint = "x-intercept",
              rsq = bquote(R^2),
              rsqaured = bquote(R^2),
              slope = "slope")
  if(outvar=="pval" |outvar=="pvalue"){
    sig_f<-as.factor(sig_mask)
    plt<-ggplot()+
      geom_spatraster(data =sig_f)+  
      scale_fill_manual(values = c("darkgray","black"), name = NULL,
                        labels = c("p>0.1","p<0.1"),
                        na.value = "white")
  }else{
    print("Masking to CONUS")
    varrast<-mask(varrast,sig_mask, maskvalue=NA,updatevalue = NA)
    print("Masking to significant")
    varrast_sig<-mask(varrast,sig_mask, maskvalue=0,updatevalue = NaN)
  }
  if(outvar=="xint"|outvar=="x-intercept"){
    plt<-ggplot()+geom_spatraster(data =varrast_sig)+
      scale_fill_continuous(na.value = "white",
                           name = NULL)
  }
  if(outvar=="rsq"|  outvar =="rsqaured"){
    plt<-ggplot()+geom_spatraster(data =varrast_sig)+
      scale_fill_viridis_c(limits = c(0,1),na.value = "white",
                           name = NULL)
  }
  if(outvar=="slope"){
      print("Calculating bounds")
      lims = global(varrast, fun=quantile, probs = c(0.01,0.99), 
                      na.rm = TRUE)
      plt<-ggplot()+geom_spatraster(data =varrast_sig)+
        scale_fill_gradient2(low="darkorange3", mid="white", high="orchid4",
                             midpoint = 0,
                             limits = unlist(lims),
                             oob=squish,na.value = "white",
                             name = NULL)
  }
  plt<-plt+labs(subtitle = tit)+theme_void()+
    theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"),
          legend.key.size = unit(0.35, 'cm'))
  # add park boundary or state outlines
  if(is.null(park_extent)){
    plt<-plt+geom_spatvector(data = states, col = "gray",
                             show.legend = FALSE,fill = NA, lwd = 0.25)
  }
  
  return(plt)
}

# *_rsq.tif = rsquared
# *_pval.tif = p-value 
# *_slope.tif = slope 


#### Read in park boundary
park <- "BLCA"
park_extent<-get_ext(park)
park_boundary<-get_ext(park, boundary=TRUE)



WBS = c("accumswe","aet","deficit","pet","runoff","soilwater")
WBS = c("soilwater")
par(mfrow = c(3,2))

for(wb in WBS){
  print(wb)
  
  plt1<-piv_out_plot(wb = wb,vi = "EVI",outvar = "pval", park_extent = park_extent)
  plt1<-plt1+labs(title=wb)
  print(plt1)
#  plt2<-piv_out_plot(wb = wb,vi = "EVI",outvar = "rsq", park_extent = park_extent)
#  plt3<-piv_out_plot(wb = wb,vi = "EVI",outvar = "slope", park_extent = park_extent)
#  plt4<-piv_out_plot(wb = wb,vi = "EVI",outvar = "xint", park_extent = park_extent)
  
}
pdf(here("results",park,paste0(park," figures"),
         paste0("piv_out_",wb,"_",park,".pdf")),
    width =  10,height = 5)

ggpubr::ggarrange(plt1,plt2,plt3,plt4)

dev.off()


# for BLCA
plist<-lapply(WBS,function(wb) piv_out_plot(wb = wb,vi = "EVI",outvar = "pval", park_extent = park_extent)+
                labs(subtitle=wb))
plist<-lapply(plist, function(x) x+geom_spatvector(data=park_boundary,
                                                   col="red", fill=NA,  
                                                   lwd=1))
ggarrange(plotlist=plist, common.legend = TRUE,  legend="right")
