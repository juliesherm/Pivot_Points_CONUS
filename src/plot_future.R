###############################################################
# plot one of the outputs of pivot point analysis for all of CONUS
# save to pdf
# pivot var
setwd("C:/Users/Julie/Documents/MATH/NPS SIP/pivot_pts")

#setwd("/home/1024/ma/sherman/NPS_SIP/pivot_pts")

library(here)
library(terra)
library(dplyr)
library(ggpubr)
library(tidyterra)
library(scales)
library(ggnewscale)
terraOptions(memfrac = 0.7)
WBS = c("accumswe","aet","deficit","pet","soilwater")
wb = "deficit"

rcp = "rcp45"
era = "midcent"
measure = "mag"
win = NULL
filter_nlcd = TRUE
# *_total.tif = total years under/over piv
# *_mag.tif = integrated (delta) vi 
# *_longest.tif = longest consecutive years under/over piv 


future_plot<-function(wb = "deficit",
                      rcp = "rcp45",
                      era = "midcent",
                      measure = "total",
                      win = NULL,
                      filter_nlcd = TRUE,
                      q = 0.5){
  print(wb)
  gc()
  if(measure=="last"){
    future<-terra::rast(here("results","CONUS","future_comparison",
                             paste0(rcp,"_",wb,"_EVI_",measure,"_CONUS.tif")),
                        win = win)
    era = "none"
  }else{
    if(!is.null(win)){
      #  FIX THIS..  HOW MANY TILES?
      tile<-get_ext(ext(win))
      xoff<-tile[[1]][1]
      yoff<-tile[[1]][2]
      xoff<-tile[[2]][1]
      yoff<-tile[[2]][2]
      future<-terra::rast(here("results","CONUS","future_comparison","tiled",
                               paste0(rcp,"_",wb,"_EVI_",era,"_",measure,"_",xoff,"_",yoff,".tif")),
                          win = win)
    }else{
      future<-terra::rast(here("results","CONUS","future_comparison",
                               paste0(rcp,"_",wb,"_EVI_",era,"_",measure,
                                      ".tif")))
    }
    #each layer is a model (in alphabetical order)
  }
  print(measure)
  
  # mask to significant p-vals
  fname<-ifelse(filter_nlcd,"sig_pvals_masked.tif","sig_pvals.tif")
  sig_mask<-terra::rast(here("results","CONUS","pivot_pts",
                             fname),
                        win = win) 
  sig_mask<-subset(sig_mask,paste0("EVI_",wb,"_pval"))
  
  print("Masking to CONUS")
  future<-mask(future,sig_mask, maskvalue=NA,updatevalue = NA,
               verbose = TRUE, steps = 2611)
  if(measure =="mag"){
    print("Masking to significant")
    future_sig<-mask(future,sig_mask, maskvalue=0,updatevalue = NA)
    rm(sig_mask)
    rm(future)
    ##Convert to per-year amount
    nyears = 35
    future_sig<-future_sig/nyears
    vi_mean<-terra::rast(here("results","CONUS","integrated_vi_gs","mean_vi",
                              "meanEVI.tif"), win = win)
    
    print("Scaling by historic integrated VI mean")
    vi_mean<-vi_mean/10000
    future_sig_rel<-future_sig/vi_mean+1
    rm(vi_mean)
    rm(future_sig)
  }else if(measure=="last"){
    future_sig_rel<-future
  }else{
    ##Convert to a percentage of years
    nyears = 35
    future_sig_rel<-future/nyears
  }
  
  print("Computing model quantiles")
  #breaks here for CONUS
  bounds = quantile(future_sig_rel, q ,na.rm = TRUE)
  rm(future_sig_rel)
  names(bounds)<-"quant"
  #lims = global(bounds,quantile, probs = c(0.05,0.95),na.rm = TRUE)
  #defining plot titles
  eranm<-switch(era,
              midcent="2030-2064",
              endcent="2065-2099",
              none = "")
  plt<-ggplot()+theme_void()+
    geom_spatraster(aes(fill = quant),data = bounds,
                    maxcell=1e7)+
    theme(plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm"),
          legend.key.size = unit(0.5, 'cm'))
  keyname<-switch(measure,
                  total = "",#years above pivot",
                  duration = "years above pivot",
                  last = "",#"last `good' year",
                  mag = "Rel. Veg.")
  lims<-switch(measure,
                  total = c(0,1),
                  duration = c(0,1),
                  last = c(2030,2100),
                  mag = c(0.5,1.25))
  colfil<-switch(measure,
                 total = scale_fill_gradient2(low="darkblue", mid="white", high="darkred", midpoint = 0.5, 
                                              na.value="white",limits = lims,
                                              name = keyname),
                 duration = scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0.5, 
                                                 na.value="white",limits = lims, 
                                                 name = keyname),
                 last = scale_fill_viridis_c(direction = -1, name = keyname,limits = lims,na.value="white"),
#                 mag = scale_fill_stepsn(low="brown", mid="cornsilk3", high="forestgreen", midpoint = 1, 
#                                            name = keyname,na.value="white",
#                                            oob = squish,limit = lims,
#                                            breaks = seq(0.5, 1.25, length.out = 11), 
#                                            labels = c("<0.5","0.75", "1", ">1.25")))
                mag = scale_fill_gradient2(low="brown", mid="cornsilk3", high="forestgreen", midpoint = 1, 
                                            name = keyname,na.value="white",
                                            oob = squish,limit = lims,
                                            breaks = seq(0.5, 1.25, length.out = 4), 
                                            labels = c("<0.5","0.75", "1", ">1.25")))
  plt<-plt+colfil
  if(!is.null(win)){
    plt<-plt+ geom_spatvector(data = park_boundary, aes(col = GEOMETRYID), 
                              fill = NA,lwd = 1, show.legend = FALSE,
                              col = "black")+
      theme(axis.text.x=element_blank(),axis.text.y=element_blank(),
            axis.ticks=element_blank())
  }
  rm(bounds)
  return(plt)
}
################3
## labeling plots
midlab<-ggplot()+geom_text(aes(x = 0, y = 0, label = "2030-2064"), 
                           parse  = FALSE, vjust = 1)+theme_void()
endlab<-ggplot()+geom_text(aes(x = 0, y = 0, label = "2065-2099"), 
                           parse  = FALSE, vjust = 1)+theme_void()

r45lab<-ggplot()+geom_text(aes(x = 0, y = 0, label = "RCP 4.5"), 
                           parse  =  FALSE,vjust = 0.5, hjust = 1, angle = 90)+theme_void()
r85lab<-ggplot()+geom_text(aes(x = 0, y = 0, label = "RCP 8.5"), 
                           parse  =  FALSE,vjust = 0.5, hjust = 1, angle = 90)+theme_void()

medlab<-ggplot()+geom_text(aes(x = 0, y = 0, label = "median"), 
                           parse  = FALSE, vjust = 1)+theme_void()
minlab<-ggplot()+geom_text(aes(x = 0, y = 0, label = "minimum"), 
                           parse  = FALSE, vjust = 1)+theme_void()
maxlab<-ggplot()+geom_text(aes(x = 0, y = 0, label = "maximum"), 
                           parse  = FALSE, vjust = 1)+theme_void()
blankplt<-ggplot()+theme_void()
#########################################

source(here("src","get_ext.R"))
park <- "CONUS"
park_ext<-get_ext(park)
win<-park_ext
park_boundary<-get_ext(park, boundary = TRUE)
### read in some extra plotting things
cultivated<-terra::rast(here("results","CONUS","nlcd_masked.tif"))
cultivated<-as.factor(cultivated)
states<-terra::vect(here("raw_data","states","s_18mr25.shp"))
states<-crop(states,  cultivated)

### set up the loops
eras = c("midcent","endcent")
rcps = c("rcp45","rcp85")
measure<-"total"
quant=0.5
plist<-list()
i=1
wb= "soilwater"
for(rcp in rcps){
  for(era in eras){
    plist[[i]]<-future_plot(wb = wb,rcp = rcp,era = era,
                     measure = measure,win = park_ext,
                     q= quant)
    # q defines the quantile to look at
    i  = i+1
  }
}
# add gray over the cultivated areas
mask_cult<-TRUE
if(mask_cult){
  plist<-lapply(plist, function(x) x+new_scale_fill()+
                  geom_spatraster(data=cultivated,aes(fill=`NLCD Land Cover Class`),
                                  show.legend=FALSE, maxcell=1e7)+
                  scale_fill_manual(limits=c("1","0"),values=c("gray30","transparent")))
}
# add state lines
if(park=="CONUS"){
  plist<-lapply(plist, function(x) x+geom_spatvector(data = states, col = "gray",
                                                     show.legend = FALSE,fill = NA, lwd = 0.5))
}

metnm<-switch(measure,
  total="Frequency",
  mag="Intensity",
  last="Last",
  duration="Duration")

qnm<-switch(as.character(quant),
            '0.5'="(median)",
            '0'="(minimum)",
            '1'="(maximum)")

pdf(here("results","CONUS","CONUS figures",
         paste0("all_",wb,"_",measure,"_",qnm,"_",park,".pdf")),
    width =  10,height = 5)

if(length(plist)==4){
  plistnew<-list(blankplt, midlab, endlab, r45lab, plist[[1]],plist[[2]],
                 r85lab, plist[[3]],plist[[4]])
  
  allplots<-ggarrange(plotlist = plistnew, ncol=3, nrow=3, 
            common.legend = TRUE, legend="right", 
            widths = c(2,10,10), heights = c(1,7,7))
  annotate_figure(allplots, fig.lab = paste0(metnm,' \n ', qnm),
                  fig.lab.size = 16, fig.lab.face="bold")
}
dev.off()


#########################
# look at several quantiles...
measure<-"last"
era <- "endcent"
qs<- c(0,0.5,1)
plist<-list()
i=1
wb= "deficit"
for(rcp in rcps){
  for(q in qs){
    plist[[i]]<-future_plot(wb = wb,rcp = rcp,era = era,
                            measure = measure,win = park_ext,
                            q= q)
    # q defines the quantile to look at
    i  = i+1
  }
}

pdf(here("results","CONUS","CONUS figures",
         paste0("all_",wb,"_",measure,"_",park,".pdf")),
    width =  12,height = 5)
if(length(plist)==6){
  plistnew<-list(blankplt, minlab, medlab,maxlab, r45lab, plist[[1]],plist[[2]],plist[[3]],
                 r85lab, plist[[4]],plist[[5]],plist[[6]])
  
  allplots<-ggarrange(plotlist = plistnew, ncol=4, nrow=3, 
                      common.legend = TRUE, legend="right", 
                      widths = c(2,10,10,10), heights = c(1,7,7))
  annotate_figure(allplots, fig.lab = 'Last',
                  fig.lab.size = 16, fig.lab.face="bold")
  
}
dev.off()
