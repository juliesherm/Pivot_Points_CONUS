## fit random forest to black canyon pixels
## provide better fits, which can be forecast to future
setwd("C:/Users/Julie/Documents/MATH/NPS SIP/pivot_pts")
library(here)
library(terra)
library(stringr)
library(ggplot2)
library(tidyr)
library(tidyterra)
library(dplyr)
source(here("src","get_ext.R"))

park <- "BLCA"

#### Read in park boundary
park_extent<-get_ext(park)
park_boundary<-get_ext(park, boundary = TRUE)
#determine  which tile  the  park is in
tile<-get_ext(park_extent)[[1]]
xoff = tile[1]
yoff = tile[2]


#read in wb
wb_files<-list.files(here("results","CONUS","wb_proj","wateryear","bilinear","tiled"),pattern = paste0("*_",xoff,"_",yoff,".tif"), full.names = TRUE)
wbs<-str_split_i(wb_files,"_",-4)
gm_files<-list.files(here("results","CONUS","gm_proj","bilinear","tiled"),pattern = paste0("*_",xoff,"_",yoff,".tif"), full.names = TRUE)
gms<-str_split_i(str_split_i(gm_files,  "_",-3),"/",-1)
wb_files<-c(wb_files, gm_files)
rast_wb<-lapply(wb_files,terra::rast,  win = park_extent)
rast_wb<-lapply(rast_wb,function(x) x/10)
wbs<-c(wbs,gms)
names(rast_wb)<-wbs

#read in vi
vi = "EVI"
rast_vi<-terra::rast(list.files(here("results","CONUS","integrated_vi_gs","window_14days","cutoff_2.5","bilinear",vi,"tiled"),
pattern = paste0("*_",xoff,"_",yoff,".tif"), full.names = TRUE),  win = park_extent)
rast_vi<-rast_vi/10000

## mask out  cultivate land (82), water (11), and barren (31)
land_class_proj<-terra::rast(here("raw_data","nlcd_2021","nlcd_2021_proj","land_cover_2021_proj.tif"),
                             win = park_extent)
veg_mask<-land_class_proj$`NLCD Land Cover Class`>31 & (land_class_proj$`NLCD Land Cover Class`<82) 

rast_vi<-mask(rast_vi, veg_mask,maskvalues=FALSE,updatevalue=NA)
rast_wb<-lapply(rast_wb,function (x) mask(x, veg_mask,maskvalues=FALSE,updatevalue=NA))

#### identify just PJ areas
# national park service vegetation map
nps_class<-terra::vect(here("raw_data","vegmaps","blcageodata",
                            "blcageodata.gdb"),
                       layer = "fcl_Veg_Polys")

nps_class <- terra::project(nps_class,crs(rast_vi))



# turn into dataframe with x and y coordinates
vi_all<-as.data.frame(rast_vi, xy = TRUE)

# extract the wb values for the locations
wb_all <- lapply(rast_wb,terra::extract, vi_all[,c("x","y")], ID=FALSE, xy = TRUE)


# put into one long dataframe
wb = wbs[1]
names(wb_all[[wb]])<-paste0(str_split_i(names(wb_all[[wb]]),"_",-5),str_split_i(names(wb_all[[wb]]),"_",-1))
names(wb_all[[wb]])[(length(names(wb_all[[wb]]))-1):length(names(wb_all[[wb]]))]<-c("x","y")

alldat<-wb_all[[wb]] %>%
  pivot_longer(
    cols = starts_with(wb),
    names_to = "year",
    names_prefix = wb,
    values_to = wb,
    values_drop_na = FALSE
  )

for(wb in wbs[-1]){
  print(wb)
  if(wb=="precip"| wb=="tavg"| wb=="tmin"| wb=="tmax"){
    names(wb_all[[wb]])<-paste0(str_split_i(names(wb_all[[wb]]),"_",-4),str_split_i(names(wb_all[[wb]]),"_",-1))
    names(wb_all[[wb]])[(length(names(wb_all[[wb]]))-1):length(names(wb_all[[wb]]))]<-c("x","y")
  }else{
    names(wb_all[[wb]])<-paste0(str_split_i(names(wb_all[[wb]]),"_",-5),str_split_i(names(wb_all[[wb]]),"_",-1))
    names(wb_all[[wb]])[(length(names(wb_all[[wb]]))-1):length(names(wb_all[[wb]]))]<-c("x","y")
  }

  alldat<-merge(alldat,wb_all[[wb]] %>%
  pivot_longer(
    cols = starts_with(wb),
    names_to = "year",
    names_prefix = wb,
    values_to = wb,
    values_drop_na = FALSE
  ))
}

#add in vi data
names(vi_all)<-paste0(str_split_i(names(vi_all),"_",-5),str_split_i(names(vi_all),"_",-1))
names(vi_all)[1:2]<-c("x","y")

alldat<-merge(alldat,vi_all %>%
  pivot_longer(
    cols = starts_with(vi),
    names_to = "year",
    names_prefix = vi,
    values_to = vi,
    values_drop_na = FALSE
  ))
alldat$year = as.numeric(alldat$year)+1999


##################################################################################
## build the  model
##################################################################################

#random forest
library(randomForest)
#subset to variables of interest
wbs = c("accumswe","aet","deficit","pet","runoff","soilwater","tavg","precip","tmax","tmin")
subdat<-subset(alldat,select = c(wbs,"x","y","EVI","year"))
#split data into training and testing
set.seed(4242)
train <- subdat %>% sample_frac(0.7, replace = FALSE)
test <- setdiff(subdat, train)
# Train the models


modlist<-list()
###
# all BLCA
print("Fitting model with only water balance")
modlist[[1]] <- randomForest(EVI~accumswe+deficit+aet+pet+runoff+soilwater, 
                             data = train, importance=TRUE, 
                             ntree=500,mtry = 3,
                             na.action = na.omit)
print("Fitting model with only water balance and x/y coordinates")
modlist[[2]] <- randomForest(EVI~accumswe+deficit+aet+pet+runoff+soilwater+x+y, 
                             data = train, importance=TRUE, 
                             ntree=500,mtry = 3,
                             na.action = na.omit)

print("Fitting model with all vars  but no x/y") #don't use x, y and year? given the idea that wb <-> vi
modlist[[3]] <- randomForest(EVI~accumswe+deficit+aet+pet+runoff+soilwater+tmin+tmax+precip, 
                             data = train, importance=TRUE, 
                             ntree=500,mtry = 3,
                             na.action = na.omit)
print("Fitting model with all vars and x/y")
modlist[[4]]<- randomForest(EVI~accumswe+deficit+aet+pet+runoff+soilwater+tmin+tmax+precip+x+y, 
                            data = train, importance=TRUE,
                            ntree=500,mtry = 3,
                            na.action = na.omit)

saveRDS(modlist, file=here("results","BLCA","RF_mods.RData"))

#print("Fitting model with bio vars and x/y")
#modlist[[5]]<- randomForest(EVI~deficit+aet+runoff+soilwater+tmin+tmax+precip+x+y, 
#                            data = train, importance=TRUE,
#                            ntree=500,mtry = 3,
#                            na.action = na.omit)

# just pj  in BLCA
##################################################################################
# Results
##################################################################################

# Plot the variable importance
pdf(here("results",park,"futwb_varImportance.pdf"))
varImpPlot(modlist[[1]])
varImpPlot(modlist[[2]])
varImpPlot(modlist[[3]])
varImpPlot(modlist[[4]])
dev.off()


# predict over all of BLCA (30% of which was withheld during training)
alldat$pred1 <- predict(modlist[[1]], alldat)
alldat$pred2 <- predict(modlist[[2]], alldat)
alldat$pred3 <- predict(modlist[[3]], alldat)
alldat$pred4 <- predict(modlist[[4]], alldat)
test$pred1 <- predict(modlist[[1]], test)
test$pred2 <- predict(modlist[[2]], test)
test$pred3 <- predict(modlist[[3]], test)
test$pred4 <- predict(modlist[[4]], test)

saveRDS(alldat, file=here("results","BLCA","RF_preds_all.RData"))
saveRDS(test, file=here("results","BLCA","RF_preds_test.RData"))
##################################################################
#  or just read in the saved data
##################################################################
alldat<-readRDS(here("results","BLCA","RF_preds_all.RData"))
test<-readRDS(here("results","BLCA","RF_preds_test.RData"))

### rmse over test data
print(sqrt(mean((test$pred1-test$EVI)^2)))
#[1] 4.998562
print(sqrt(mean((test$pred2-test$EVI)^2)))
#[1] 3.645936
print(sqrt(mean((test$pred3-test$EVI)^2)))
#[1] 7.367839
print(sqrt(mean((test$pred4-test$EVI)^2)))
#[1] 3.467265

### R2 over test data
summary(lm(test$EVI~test$pred1))
#0.892
summary(lm(test$EVI~test$pred2))
#0.943
summary(lm(test$EVI~test$pred3))
#0.760
summary(lm(test$EVI~test$pred4))
#0.948

#convert to raster
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

#each layer is a year
yrs <-unique(alldat$year)
yr_pred_df1 <- split(alldat[,c("x","y","pred1")], alldat$year)
yr_pred_df2 <- split(alldat[,c("x","y","pred2")], alldat$year)
yr_pred_df3 <- split(alldat[,c("x","y","pred3")], alldat$year)
yr_pred_df4 <- split(alldat[,c("x","y","pred4")], alldat$year)
yr_pred_rast1 <- lapply(yr_pred_df1, \(i) rast(i, type="xyz", crs = crs(ref_crs))) 
yr_pred_rast2 <- lapply(yr_pred_df2, \(i) rast(i, type="xyz", crs = crs(ref_crs))) 
yr_pred_rast3 <- lapply(yr_pred_df3, \(i) rast(i, type="xyz", crs = crs(ref_crs))) 
yr_pred_rast4 <- lapply(yr_pred_df4, \(i) rast(i, type="xyz", crs = crs(ref_crs))) 
pred_rast1 <- rast(yr_pred_rast1)
pred_rast2 <- rast(yr_pred_rast2)
pred_rast3 <- rast(yr_pred_rast3)
pred_rast4 <- rast(yr_pred_rast4)

# compute rmse
er1<-pred_rast1-rast_vi
er1<-er1^2
er1<-mean(er1)
er1<-sqrt(er1)

er2<-pred_rast2-rast_vi
er2<-er2^2
er2<-mean(er2)
er2<-sqrt(er2)

er3<-pred_rast3-rast_vi
er3<-er3^2
er3<-mean(er3)
er3<-sqrt(er3)

er4<-pred_rast4-rast_vi
er4<-er4^2
er4<-mean(er4)
er4<-sqrt(er4)

#percent error
pcter1<-mean(abs(pred_rast1-rast_vi)/rast_vi)
pcter2<-mean(abs(pred_rast2-rast_vi)/rast_vi)
pcter3<-mean(abs(pred_rast3-rast_vi)/rast_vi)
pcter4<-mean(abs(pred_rast4-rast_vi)/rast_vi)
park_boundary<-get_ext("BLCA",boundary = TRUE)

## for plotting
#read in value polygons (from park managers)
value_polys <- terra::vect(here("raw_data", "NCPN_PJ_Values_2024.gdb"))
value_polys<- subset(value_polys, value_polys$Park==park)
value_polys <- terra::project(value_polys,crs(ref_crs))

#plot percent error for JUST model 4 ()
pdf(here("results",park,"pcter_mod_2.pdf"))
  #plot(pcter4,  range = c(0,1.4), main = "mean proportion error")
  #lines(park_boundary)
  ggplot()+theme_void()+
    geom_spatraster(data = pcter2, aes(fill = mean)) + 
    scale_fill_continuous(type = "viridis", name = "MPE") +
    geom_spatvector(data = value_polys, aes(col = ValueRate), fill = NA,
                    lwd = 1,show.legend = FALSE) + 
    scale_color_manual(values=c("orange", "red", "#56B4E9"))+
    geom_spatvector(data = park_boundary, fill = NA, col = "lightgray",  lwd = 1)
  
dev.off()

#############################
# pj  emphasized plot

nps_class<-terra::vect(here("raw_data","vegmaps","blcageodata",
                            "blcageodata.gdb"),
                       layer = "fcl_Veg_Polys")
nps_class <- terra::project(nps_class,crs(ref_crs))
nps_class<-terra::intersect(nps_class, park_extent)
pj_classes<-subset(nps_class,  nps_class$Map_Unit_ID %in% 
                     c(69,71,72,74,78,87))

sum(expanse(terra::intersect(pj_classes, park_boundary), 
            unit = 'ha'))
#4742 hectares of PJ classification out of 
sum(expanse(park_boundary, unit = 'ha'))
#12637.5  total hectares
# accounting  for  approximately 
sum(expanse(terra::intersect(pj_classes, park_boundary), 
            unit = 'ha'))/sum(expanse(park_boundary, unit = 'ha'))
#37.5% of the land in BLCA (approx 3/8)

pcter2_pj_masked<-mask(pcter2, pj_classes)
pdf(here("results",park,"pcter_pj_mod_2.pdf"))

ggplot()+theme_void()+
  geom_spatraster(data = pcter2_pj_masked, aes(fill = mean),alpha = 1)+
  geom_spatraster(data = pcter2, aes(fill = mean), alpha = 0.4)+
  guides(alpha = "none") + 
  scale_fill_continuous(type = "viridis", name = "MPE") +
  geom_spatvector(data = value_polys, aes(col = ValueRate), fill = NA,lwd = 1,show.legend = FALSE) + 
  scale_color_manual(values=c("orange", "red", "#56B4E9"))+
  geom_spatvector(data = park_boundary, fill = NA, col = "lightgray",  lwd = 1)+
  geom_spatvector_text(aes(label = poly_num),data = value_polys,
                       fontface = "bold",
                       color = "white")
dev.off()

# plot rmse of ALL 4 models
pdf(here("results",park,"figures","rmse_4_mods.pdf"))
ggplot()+theme_bw()+
 geom_spatraster(data = er1, aes(fill = mean)) + 
  scale_fill_continuous(type = "viridis") +
 geom_spatvector(data = value_polys, aes(col = ValueRate), fill = NA,lwd = 1) + 
 scale_color_manual(values=c("orange", "red", "#56B4E9"))+
 geom_spatvector(data = park_boundary, fill = NA, col = "lightgray",  lwd = 1)
ggplot()+theme_bw()+
  geom_spatraster(data = er2, aes(fill = mean)) + 
  scale_fill_continuous(type = "viridis") +
  geom_spatvector(data = value_polys, aes(col = ValueRate), fill = NA,lwd = 1) + 
  scale_color_manual(values=c("orange", "red", "#56B4E9"))+
  geom_spatvector(data = park_boundary, fill = NA, col = "lightgray",  lwd = 1)
ggplot()+theme_bw()+
  geom_spatraster(data = er3, aes(fill = mean)) + 
  scale_fill_continuous(type = "viridis") +
  geom_spatvector(data = value_polys, aes(col = ValueRate), fill = NA,lwd = 1) + 
  scale_color_manual(values=c("orange", "red", "#56B4E9"))+
  geom_spatvector(data = park_boundary, fill = NA, col = "lightgray",  lwd = 1)
ggplot()+theme_bw()+
  geom_spatraster(data = er4, aes(fill = mean)) + 
  scale_fill_continuous(type = "viridis") +
  geom_spatvector(data = value_polys, aes(col = ValueRate), fill = NA,lwd = 1) + 
  scale_color_manual(values=c("orange", "red", "#56B4E9"))+
  geom_spatvector(data = park_boundary, fill = NA, col = "lightgray",  lwd = 1)

dev.off()

global(er1, rms)

##################################################3
# plot predicted versus real
yr = sample(1:nlyr(rast_vi),1)
par(mfrow = c(1,3))

pdf(here("results",park,"figures","rm_fit_rmse_no_xy.pdf"))
par(mfrow = c(1,3))
#true
ggplot()+theme_bw()+
  geom_spatraster(data = subset(rast_vi,yr),
                  show.legend = NA)+
  labs(title="Growing Season EVI")+
  scale_fill_continuous(limits = c(10,120)) +
  geom_spatvector(data = value_polys, aes(col = ValueRate), 
                  fill = NA,lwd = 1,
                  show.legend = FALSE) + 
  scale_color_manual(values=c("orange", "red", "#56B4E9"))+
  geom_spatvector(data = park_boundary,  col = "lightgray",
                  fill = NA,  lwd = 1)+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),      
        axis.ticks=element_blank())
#predicted
ggplot()+theme_bw()+
  geom_spatraster(data = subset(pred_rast,yr),
                  show.legend = FALSE,)+
  labs(title="RF Predicted EVI")+
  scale_fill_continuous(limits = c(10,120)) +
  geom_spatvector(data = value_polys, aes(col = ValueRate), 
                  fill = NA,lwd = 1,
                  show.legend = FALSE) + 
  scale_color_manual(values=c("orange", "red", "#56B4E9"))+
  geom_spatvector(data = park_boundary,  col = "lightgray",
                  fill = NA,  lwd = 1)+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),      
        axis.ticks=element_blank())

#rmse
ggplot()+theme_bw()+
  geom_spatraster(data = er, aes(fill = mean)) + 
  labs(title="RMSE")+
  scale_fill_continuous(type = "viridis") +
  geom_spatvector(data = value_polys, aes(col = ValueRate), 
                  fill = NA,lwd = 1,
                  show.legend = FALSE) + 
  scale_color_manual(values=c("orange", "red", "#56B4E9"))+
  geom_spatvector(data = park_boundary, fill = NA, col = "lightgray",  lwd = 1)+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),      
        axis.ticks=element_blank())

#mean percent  error
ggplot()+theme_bw()+
  geom_spatraster(data = er, aes(fill = mean)) + 
  labs(title="MPE")+
  scale_fill_continuous(type = "viridis") +
  geom_spatvector(data = value_polys, aes(col = ValueRate), 
                  fill = NA,lwd = 1,
                  show.legend = FALSE) + 
  scale_color_manual(values=c("orange", "red", "#56B4E9"))+
  geom_spatvector(data = park_boundary, fill = NA, col = "lightgray",  lwd = 1)+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),      
        axis.ticks=element_blank())

dev.off()



##################################################################################
# Predict with future data
# Save each GCM predictions
##################################################################################
library(here)
library(terra)
library(stringr)
library(ggplot2)
library(tidyr)
library(tidyterra)
library(dplyr)
library(randomForest)
source(here("src","get_ext.R"))

park <- "BLCA"
rast_vi<-terra::rast(list.files(here("results","CONUS","integrated_vi_gs","window_14days","cutoff_2.5","bilinear",vi,"tiled"),
                                pattern = paste0("*_",xoff,"_",yoff,".tif"), full.names = TRUE),  win = park_extent)


modlist<-readRDS(modlist, file=here("results","BLCA","RF_mods.RData"))
class(modlist[[2]])<-"randomForest"
#read in future models of wb  variables
rcp = "rcp85"
wbs<-"accumswe"
mods<-list.files(here("results","CONUS","wb_proj","futures","tiled"),
                 pattern = paste0("*_",rcp,"_",wbs[1],"_",xoff,"_",yoff,".tif"))
mods<-str_split_i(mods,"_",-5)
#mods<-mods[1]
for(m in mods){
  print(m)
  mod_files <- list.files(here("results","CONUS","wb_proj","futures","tiled"),pattern = paste0("V_1_5_annual_",m,"_",rcp,"_*"), full.names = TRUE)
  mod_files<- mod_files[str_ends(mod_files, paste0("_",xoff,"_",yoff,".tif"))]
  modrast<-lapply(mod_files, terra::rast, win = park_extent)
  modrast<-lapply(modrast,function(x) x/10)
  wbs<-str_split_i(mod_files,"_",-3)
  names(modrast)<-wbs
  moddf<-as.data.frame(modrast[[wbs[1]]], xy = TRUE)
  moddf<-moddf%>%pivot_longer(
      cols = contains(wbs[1]),
      names_to = "year",
      names_prefix = paste0("V_1_5_annual_",m,"_",rcp,"_",wbs[1],"_",xoff,"_",yoff,"_"),
      values_to = wbs[[1]],
      values_drop_na = FALSE)
  if(m=="BNU-ESM" & rcp=="rcp45"){
    moddf$year<-as.numeric(moddf$year)-24
  }
  for(wb in wbs[-1]){
    print(wb)
    tp<-as.data.frame(modrast[[wb]], xy = TRUE)
    #print(ncol(tp)-2)
    moddf<-merge(moddf,tp%>%pivot_longer(
        cols = contains(wb),
        names_to = "year",
        names_prefix = paste0("V_1_5_annual_",m,"_",rcp,"_",wb,"_",xoff,"_",yoff,"_"),
        values_to = wb,
        values_drop_na = FALSE
      ), all = TRUE)
  }
  moddf$year<-as.numeric(moddf$year)+2029
  #moddf$tavg <-NA
  #moddf$tmin<-NA
  #moddf$tmax<-NA
  #moddf$precip<-NA
  print("predicting")
  #use random forest trained above to predict
  moddf$fut_pred <- predict(modlist[[2]], moddf)
  yrs <-unique(moddf$year)
  yr_pred_df <- split(moddf[,c("x","y","fut_pred")], moddf$year)
  yr_pred_rast <- lapply(yr_pred_df, \(i) rast(i, type="xyz", crs = crs(rast_vi))) 
  pred_rast <- rast(yr_pred_rast)
  print("Saving")
  #save predicted futures
  writeRaster(pred_rast, here("results",park,"RF_future_fit",paste0("withwbxy_",m,"_",rcp,".tif")), 
              overwrite = TRUE)
}

##################
## compare to linear output
##################
linear_comp<-FALSE
if(linear_comp){
  piv_out<-terra::rast(here("results","CONUS", "pivot_pts","window_14days","cutoff_2.5","bilinear","tiled",
                            "EVI",paste0("pivs_",vi,"_",wb,"_",xoff,"_",yoff,".tif")), 
                       win = park_extent)
  #the wb values used to fit the model
  wb_rast<-terra::rast(here("results","CONUS", "wb_proj","wateryear","bilinear","tiled",
                            paste0("V_1_5_wateryears_gridmet_historical_",wb,"_resampled_",xoff,"_",yoff,".tif")), 
                       win = park_extent)
  wb_rast<-wb_rast/10000
  piv_pred_rast<-piv_out['slope']*wb_rast+piv_out['y-intercept']+ mean(rast_vi)
  piv_er<-piv_pred_rast-rast_vi
  piv_er<-piv_er^2
  piv_er<-mean(piv_er)
  piv_er<-sqrt(piv_er)
  ggplot()+theme_bw()+
    geom_spatraster(data = piv_er, aes(fill = mean)) +
    scale_fill_continuous(type = "viridis") +
    geom_spatvector(data = value_polys, aes(col = ValueRate), fill = NA,lwd = 1) + 
    scale_color_manual(values=c("orange", "red", "#56B4E9"))+
    geom_spatvector(data = park_boundary, fill = NA, col = "lightgray",  lwd = 1)
  
  #plot the difference between the RF rmse and the LM rmse
  ggplot()+theme_bw()+
    geom_spatraster(data = er-piv_er, aes(fill = mean)) +
    scale_fill_continuous(type = "viridis", limits = c(-100,0)) +
    geom_spatvector(data = value_polys, aes(col = ValueRate), fill = NA,lwd = 1) + 
    scale_color_manual(values=c("orange", "red", "#56B4E9"))+
    geom_spatvector(data = park_boundary, fill = NA, col = "lightgray",  lwd = 1)
  
  ##########
  #plot the three side by side of a random year
  yr = sample(1:nlyr(rast_vi),1)
  par(mfrow = c(1,3))
  plot(subset(rast_vi,yr), range = c(10,130), axes = FALSE, main = "Growing Season EVI")
  lines(park_boundary, col = "lightgray")
  plot(subset(pred_rast,yr), range = c(10,130), axes = FALSE, main = "RF Predicted EVI")
  lines(park_boundary, col = "lightgray")
  plot(subset(piv_pred_rast,yr), range = c(10,130), axes = FALSE, main = "LM Predicted EVI")
  lines(park_boundary, col = "lightgray")
  ###
  # or plot absolute error?
  #plot(abs(subset(piv_pred_rast,yr)-subset(rast_vi,yr)), range = c(0,80), axes = FALSE, main = "RF Predicted EVI")
  #lines(park_boundary, col = "lightgray")
  #plot(abs(subset(pred_rast,yr)-subset(rast_vi,yr)), range = c(0,80), axes = FALSE, main = "RF Predicted EVI")
  #lines(park_boundary, col = "lightgray")
}

