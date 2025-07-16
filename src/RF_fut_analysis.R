#
setwd("C:/Users/Julie/Documents/MATH/NPS SIP/pivot_pts")
library(here)
library(terra)
library(stringr)
library(ggplot2)
library(tidyr)
library(tidyterra)
library(dplyr)
library(scales)
library(ggpubr)
library(ggdist)

source(here("src","get_ext.R"))
source(here("src","RF_fut_plot.R"))


park <- "BLCA"
RF<-"withwbxy"

#### Read in park boundary
park_extent<-get_ext(park)
park_boundary<-get_ext(park, boundary = TRUE)
#determine  which tile  the  park is in
#tile<-get_ext(park_extent)[[1]]
#xoff = tile[1]
#yoff = tile[2]
# For BLCA
xoff = 0
yoff = 5222

#read in value polygons (from park managers)
value_polys <- terra::vect(here("raw_data", "NCPN_PJ_Values_2024.gdb"))
value_polys<- subset(value_polys, value_polys$Park==park)
value_polys$poly_num<-1:9

#project to MODIS crs
ref_crs<-lapply(list.files(here("raw_data","MODIS","CONUS","VI_16Days_250m_v61","EVI","tiled"),
                           pattern = paste0("EVI","_2000_*"), full.names = T),
                terra::rast)
if(length(ref_crs)==0){
  ref_crs <- terra::rast(list.files(here("results","CONUS","wb_proj","wateryear","bilinear","tiled"),
                                    pattern = paste0("*_resampled_*"), full.names = T)[1])
}else{
  ref_crs = ref_crs[[1]]
}
value_polys <- terra::project(value_polys,crs(ref_crs))
rm(ref_crs)
#read in historical mean
vi_mean<-terra::rast(here("results","CONUS","integrated_vi_gs","mean_vi",
                          paste0("mean_EVI_",xoff,"_",yoff,".tif")),
                     win = park_extent)
vi_mean<-vi_mean/10000

#read in future predicted vi from the RF model
fut_files<-list.files(here::here("results",park,"RF_future_fit"), pattern = RF, full.names = TRUE)
fut_vi<-lapply(fut_files,terra::rast, win = park_extent)

# determine which  pixels to plot
#read in historic vi
historic_vi<-terra::rast(here::here("results","CONUS","integrated_vi_gs","window_14days","cutoff_2.5","bilinear","EVI","tiled",paste0("EVI_int_",xoff,"_",yoff,".tif")), win = park_extent)

##pick the pixels by clicking on a map
dev.new(noRStudioGD = TRUE)
plot(historic_vi,1)
lines(value_polys)
click(historic_vi, n = 2, cell = TRUE)

##alternatively, hardcode the pixel numbers
cell1 = 354
cell2 = 6270
cell1  = 3019
cell2 = 483 # pretty good
rw1 = rowFromCell(historic_vi,cell1)
col1 = colFromCell(historic_vi,cell1)
rw2 = rowFromCell(historic_vi,cell2)
col2 = colFromCell(historic_vi,cell2)

#data frame to plot the point
x1<-xFromCell(historic_vi, cell1)
y1<-yFromCell(historic_vi, cell1)
x2<-xFromCell(historic_vi, cell2)
y2<-yFromCell(historic_vi, cell2)
pts<-data.frame(lon = c(x1,x2),lat = c(y1,y2), pixid = c(1,2))


############################################
# Spatial Future Plotting
############################################

# plot of min, median, and  max for ONLY midcent, split by RCP

source(here("src","RF_fut_plot.R"))
quants<-c(1,0.5,0)
RF = "withwbxy"
rcps<-c("rcp45","rcp85")

i=1
plist<-list()
for(rcp in rcps){
  for(quant in quants){
    plist[[i]]<-RF_future_plot(RF = RF,rcp = rcp,
                               era = "midcent",  quant = quant,
                               pj_shade = TRUE)
    i  = i+1
  }
}
#add value poly numbers to the last panel
plist[[6]]<-plist[[6]]+
  geom_spatvector_text(aes(label = poly_num),data = value_polys,
                       fontface = "bold",
                       color = "black")
#add pixel locations for the time series plots to the second-to-last panel
plist[[5]]<-plist[[5]]+
  geom_point(aes(x=lon,y=lat, col = pixid),data = pts,
                       pch = 8,size = 1.5,
                       color = "black") 

#plot all together
print(ggarrange(plotlist = plist, 
                ncol=3, nrow=2, 
                common.legend = TRUE, legend="right"))

#################################################################
#  Pixel time series plots
#################################################################

## create dataframe for ggplotting
vidf<-data.frame(year=2001:2024, EVI = unlist(historic_vi[rw1,col1,]/10000),  
                 pix = 1,rcp = 0,  mod = "obs")
vidf<-rbind(vidf,data.frame(year=2001:2024, EVI = unlist(historic_vi[rw2,col2,]/10000),
                            pix = 2,rcp = 0,  mod = "obs"))
#add in future values from each model
for(i in 1:length(fut_vi)){
  rcp<-str_split_i(str_split_i(sources(fut_vi[[i]]),"_",-1),".tif",1)
  mod<-str_split_i(sources(fut_vi[[i]]),"_",-2)  
  vidf<-rbind(vidf, data.frame(year = as.numeric(names(fut_vi[[i]])), 
                               EVI = unlist(fut_vi[[i]][rw1,col1,]),
                               pix = 1,rcp = rcp,mod = mod))
  vidf<-rbind(vidf, data.frame(year = as.numeric(names(fut_vi[[i]])), 
                               EVI = unlist(fut_vi[[i]][rw2,col2,]),
                               pix = 2,rcp = rcp,mod = mod))
}
vidf$decade<-paste0(str_sub(as.character(vidf$year), start = 1, end =  3),  "0's")

#add in RF  fit for historic EVI
preds<-readRDS(here("results","BLCA","RF_preds_all.Rdata"))
histpred<-subset(preds,subset = x==x1 & y==y1, select = c(year,EVI,pred2))
histpred$pix<-1
histpred<-rbind(histpred, cbind(subset(preds,subset = x==x2 & y==y2, select = c(year,EVI,pred2)),
                                data.frame(pix = 2)))

##  decade-wise  density plot
pixid=1
library(ggridges)
denp<-ggplot()+geom_density(data = subset(vidf,pix==pixid),
                      aes(y = EVI, fill = rcp),  alpha = 0.4)+theme_bw()
denp+facet_grid(~decade)

denp<-ggplot(data = subset(vidf,pix==pixid),
             aes(x = EVI,y = decade, group = interaction(rcp,decade ), 
                 fill = rcp, alpha = 0.5))+
  ggridges::geom_density_ridges(scale = 1, 
    linewidth = 0.3, rel_min_height = 0.01) +
  scale_fill_manual( values = c("gray","blue","red"), guide=NULL)+
  labs(title = element_blank(),
       subtitle = paste0('P',pixid))+ylab(element_blank())+
  xlab('iEVI')+theme_bw()
denp
# time series plot
## make the plot
pixid = 1
if(pixid==2){
  plotid=9
}else{
  plotid=1
}
pdf(here('results','BLCA','BLCA figures',
         paste0('p',pixid,'_timeseries_EVI_RF_withxy.pdf')),
    width = 5,height =5*2/3)
ggplot() +
  stat_lineribbon(data = subset(vidf,pix==pixid & rcp=="rcp85"),
                  aes(x = year, y = EVI, fill_ramp = after_stat(.width)),
                  .width = ppoints(100), fill ="brown3",
                  alpha = 0.35, col = "brown3", lwd = 0.5) +
  stat_lineribbon(data = subset(vidf,pix==pixid & rcp=="rcp45"),
                  aes(x = year, y = EVI, fill_ramp = after_stat(.width)),
                  .width = ppoints(100), fill = "#2171b5",
                  alpha = 0.25,col = "#2171b5", lwd = 0.5) +
  scale_fill_ramp_continuous(#range = c(1, 0), 
                             guide = NULL)+
  geom_line(data = subset(vidf,pix==pixid & rcp=="0"),
            aes(x = year, y = EVI), lwd = 1.25)+
  geom_line(data = subset(histpred, pix==pixid),
            aes(x = year, y = pred2),lty = 2, lwd = 1, col= "gray30")+
  xlim(c(2000,2100))+
  ylim(c(20,54))+
  labs(title = element_blank(),
       subtitle = paste0('P',plotid))+xlab(element_blank())+
  ylab("iEVI")+
  theme_bw()
dev.off()

#################################################
## Under construction
#################################################

#######  pixel-wise WB time series
# is this interesting? Not obvious which wb to choose, etc.
# most important??
#read in historic WB
historic_wb<-terra::rast(here::here("results","CONUS","wb_proj","wateryear","bilinear","tiled",
                                    paste0("V_1_5_wateryears_gridmet_historical_runoff_resampled_",xoff,"_",yoff,".tif")),
                         win = park_extent)

#read in future predicted wb from the GCMs
fut_wb_files<-list.files(here::here("results","CONUS","wb_proj","futures","tiled"),
                         pattern = "runoff", full.names = TRUE)
fut_wb<-lapply(fut_wb_files,terra::rast, win = park_extent)

hist_wb1<-unlist(historic_wb[rw1,col1,]/10)
hist_wb2<-unlist(historic_wb[rw2,col2,]/10)


wbdf<-data.frame(year=2001:2024, runoff = hist_wb1,  pix = 1,rcp = 0,  mod = "obs")
wbdf<-rbind(wbdf,data.frame(year=2001:2024, runoff = hist_wb2,  pix = 2,rcp = 0,  mod = "obs"))

for(i in 1:length(fut_wb)){
  rcp<-str_split_i(sources(fut_wb[[i]]),"_",-4)
  mod<-str_split_i(sources(fut_wb[[i]]),"_",-5)  
  wbdf<-rbind(wbdf, data.frame(year = 2030:2099, 
                               runoff = unlist(fut_wb[[i]][rw1,col1,]/10),
                               pix = 1,rcp = rcp,mod = mod))
  wbdf<-rbind(wbdf, data.frame(year = 2030:2099, 
                               runoff = unlist(fut_wb[[i]][rw2,col2,]/10),
                               pix = 2,rcp = rcp,mod = mod))
}

pixid = 2
ggplot() +
  stat_lineribbon(data = subset(wbdf,pix==pixid & rcp=="rcp85"),
                  aes(x = year, y = runoff, fill_ramp = after_stat(.width)),
                  .width = ppoints(10), fill ="brown3",
                  alpha = 0.25, col = "brown3", lwd = 1) +
  stat_lineribbon(data = subset(wbdf,pix==pixid & rcp=="rcp45"),
                  aes(x = year, y = runoff, fill_ramp = after_stat(.width)),
                  .width = ppoints(10), fill = "#2171b5",
                  alpha = 0.2,col = "#2171b5", lwd = 1) +
  scale_fill_ramp_continuous(range = c(1, 0), guide = NULL)+
  geom_line(data = subset(wbdf,pix==pixid & rcp=="0"),
            aes(x = year, y = runoff), lwd = 1.25)+
  xlim(c(2000,2100))+#ylim(c(0,300))+
  labs(title = element_blank(),subtitle = paste0('P',pixid))+xlab(element_blank())+
  theme_bw()



#################################################################
#  Project to value/vulnerability dimensions
#################################################################
RF = "withxy"
rcp = "rcp85"
fut_evis<-list.files(here("results",park,"RF_future_fit"),
                     pattern = regex(paste0(RF,".+",rcp,".tif")),full.names = TRUE)
fut<-lapply(fut_evis,rast)
halfcent_mag<-lapply(fut, function(x) tapp(x,c(rep(TRUE, 35),rep(FALSE, 35)),sum))
midcent_mag<-rast(lapply(halfcent_mag, function(x) subset(x, 1)/35))
mid_diff<-median(midcent_mag)/vi_mean

#summarize over the polygons
value_polys<-zonal(mid_diff, value_polys,fun = "median", as.polygons = TRUE)
names(value_polys)[names(value_polys)=='median']<-paste0("median",rcp)

#### Assign RAD decision...
value_polys$RAD<-NA
for(i in 1:length(value_polys)){
  if(value_polys$median[i]>1){ #veg expected to increase
    value_polys$RAD[i] = "accept"
    next
  }else if(value_polys$median[i]<0.9 & #veg decreasing a bit, and somewhat important
           str_sub(value_polys$ValueRate[i],1,1)=="2"){
    value_polys$RAD[i] = "direct"
  }else if(value_polys$median[i]>0.9 & #veg decreasing only a little, and most important
           str_sub(value_polys$ValueRate[i],1,1)=="3"){
    value_polys$RAD[i] = "resist"
  }else if (value_polys$median[i]>0.8 & #veg decreasing a bit, and most important
            str_sub(value_polys$ValueRate[i],1,1)=="3"){
    value_polys$RAD[i] = "direct"
  }
}

polydf<-data.frame(list(vuln = value_polys$median, val = as.numeric(substr(value_polys$ValueRate, 1, 1)),
                        poly_num = value_polys$poly_num))
polydf$val
polydf$val[is.na(polydf$val)]<-3

#### Plot in value-sensitivity plane
vulnerability <-rep(seq(0.5,1.5, length.out = 1000),3)
value<-rep(1:3, each = 1000)
man_priority<-(value)^1*(1-vulnerability/(value/0.2)^(1/10))*0.5
ggplot(mapping = aes(vulnerability, value, fill = man_priority))+
  geom_tile()+theme_bw()+ylab("value score")+
  scale_x_continuous("vulnerability", breaks = c(0.5,0.5,1.5),labels =  c("least","","most"))+
  scale_fill_gradient(low = "white", high = "red",breaks = c(min(man_priority),median(man_priority),max(man_priority)),
                      labels = c("accept", "direct", "resist"))+
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 10,
                                title = ""))+
  annotate(geom = "text", x = polydf$vuln, y = polydf$val, label = as.character(polydf$poly_num))
  
