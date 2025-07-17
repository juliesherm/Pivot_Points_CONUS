
library(here)
library(terra)

VI = "EVI"
wbs = c("accumswe","aet","deficit","pet","runoff","soilwater")

############
#filter out developed land, open water, ...
land_class_proj<-terra::rast(here("raw_data","nlcd_2021","nlcd_2021_proj","land_cover_2021_proj.tif"))
# filter by land class = 0 = "unclassified"
# filter by land class = 11 = "open water"
# filter by land class = 12 = "perennial snow/ice" ...
#veg_mask<-land_class_proj$`NLCD Land Cover Class`>40 & (land_class_proj$`NLCD Land Cover Class`<80 | land_class_proj$`NLCD Land Cover Class`>89 ) 

#sig_pvals<-terra::rast(here("results","CONUS", "pivot_pts","sig_pvals.tif"))

#sig_pvals<-mask(sig_pvals, veg_mask, maskvalue = FALSE, updatevalue = NA)
#writeRaster(sig_pvals,here("results","CONUS", "pivot_pts","sig_pvals_masked.tif"))

sig_pvals<-terra::rast(here("results","CONUS", "pivot_pts",
                            "sig_pvals_masked.tif"))
sig_pvals_wb<-subset(sig_pvals,paste0("EVI_",wbs,"_pval"))

######Analysis over ALL wb variables
#calculate the proportion of sig pixels for ALL wb variables

# count how many significant relationships per pixel
# for any of the variables
num_sig_any<-sum(sig_pvals, na.rm = TRUE)
# for just the wb variables (not gridmet)
num_sig_wb<-sum(sig_pvals_wb, na.rm = TRUE)

# count how many have at least one significant relationship
#### for any variable - 0.7982557
any_sig_any<-num_sig_any>0
global(any_sig_any,  mean, na.rm = TRUE)
# by land_class
zonal(any_sig_any,  land_class,mean, na.rm = TRUE)
#### for just the wb variables - 0.7551616
any_sig_wb<-num_sig_wb>0
global(any_sig_wb,  mean, na.rm = TRUE)
# by land_class
zonal(any_sig_wb,land_class,  mean, na.rm = TRUE)


##################################################################
#table 1
#calculate the proportion of sig pixels for EACH wb variable by land type
#store in a dataframe that can be copied  into LaTeX table
descs = c("Unclassified","Open Water","Perennial Snow/Ice",
          "Developed, Open Space","Developed, Low Intensity","Developed, Medium Intensity",
          "Developed, High Intensity","Barren Land",
          "Deciduous Forest","Evergreen Forest","Mixed Forest",
          "Shrub/Scrub","Grassland / Herbaceous","Pasture / Hay",
          "Cultivated Crops","Woody Wetlands","Emergent Herbaceous Wetlands",
          "All")
perc_sig_df<-data.frame(#land_class_num = c(unique(land_class)$`NLCD Land Cover Class`,100), 
                        land_class_description = descs)
for(wb in names(sig_pvals)){
    print(wb)
    wb_sig<-subset(sig_pvals,wb)
    sig_all_land<-global(wb_sig,  mean, na.rm = TRUE)$mean
    perc_sig_df<-cbind(perc_sig_df,
                       c(zonal(wb_sig,  land_class_proj, mean,na.rm = TRUE)[,2],
                         sig_all_land))
#    print(head(perc_sig_df))
}
perc_sig_df<-cbind(perc_sig_df,
                   c(zonal(any_sig_wb,land_class_proj,mean, na.rm = TRUE)[,2],
                            global(any_sig_wb,mean,na.rm = TRUE)$mean))
perc_sig_df<-cbind(perc_sig_df,
                   c(zonal(any_sig_any,land_class_proj,mean, na.rm = TRUE)[,2],
                     global(any_sig_any,mean,na.rm = TRUE)$mean))
names(perc_sig_df)<-paste0("{",c("Land Class",str_split_i(names(sig_pvals),"_",2),
                      "Any WB", "Any"),"}")
write.table(perc_sig_df, file = here("results","CONUS","pivot_pts",
                                     "percent_sig_table_masked.txt"), 
            sep = "& ",eol = " \\\\ ", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)



##################################################################
#for individual parks
#calculate the proportion of significant pixels for each wb variable 
if(~park==="CONUS"){
  park<-"COLM"
  park<-"BLCA"
  park_extent<-get_ext(park, boundary = FALSE)
  park_boundary<-get_ext(park, boundary = TRUE)
  sig_pvals<-terra::rast(here("results","CONUS", "pivot_pts","sig_pvals.tif"),
                         win = park_extent)
  sig_pvals_park<-mask(sig_pvals,  park_boundary)
  land_class<-terra::rast(here("raw_data","nlcd_2021","nlcd_2021_proj","land_cover_2021_proj.tif"),
                          win = park_extent)
  land_class_park<-mask(land_class,  park_boundary)
  ## percent of pixels within each land class
  freq(land_class)
  freq(land_class_park)
  
}

perc_sig_df<-data.frame(wb = str_split_i(names(sig_pvals),"_",2), percent_sig = rep(NA, length(names(sig_pvals))), 
                        percent_sig_park = rep(NA, length(names(sig_pvals))))
for(wb in 1:nlyr(sig_pvals)){
  print(wb)
  wb_sig<-subset(sig_pvals,wb)
  wb_sig_park<-subset(sig_pvals_park,wb)
  perc_sig_df[wb,"percent_sig"]<-global(wb_sig,  mean, na.rm = TRUE)$mean
  perc_sig_df[wb,"percent_sig_park"]<-global(wb_sig_park,  mean, na.rm = TRUE)$mean
  #by the top two land classes
  #print(zonal(wb_sig_park,  land_class_park, na.rm = TRUE)[c(4,5),])
  print(zonal(wb_sig_park,  land_class_park, na.rm = TRUE)[c(2,3,4),])
}
perc_sig_df<-rbind(perc_sig_df,  list(wb = "any",  percent_sig = global(sum(sig_pvals)>0, mean, na.rm = TRUE)$mean,
                                      percent_sig_park = global(sum(sig_pvals)>0, mean, na.rm = TRUE)$mean))
perc_sig_df
## save the output table
# with format options making it easily latex-able
write.table(perc_sig_df, file = here("results",park,"pivot_pts","percent_sig_table.txt"),
            sep = "& ",eol = " \\\\ ", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)

#old
##################################################################
#see the proportion of significant pixels for each wb variable by land type
wb<-"accumswe"
wb_sig<-subset(sig_pvals,paste0("EVI_",wb,"_pval"))
perc_sig<-zonal(wb_sig,  land_class, mean,na.rm = TRUE)
for(wb in wbs[-1]){
    print(wb)
    wb_sig<-subset(sig_pvals,paste0("EVI_",wb,"_pval"))
    perc_sig<-cbind(perc_sig,zonal(wb_sig,  land_class, na.rm = TRUE)[,2])
    #print(zonal(wb_sig,  land_class, na.rm = TRUE))
}
perc_sig<-cbind(perc_sig,zonal(any_sig_wb,land_class,  mean, na.rm = TRUE)[,2])
names(perc_sig)<-c("land classification",wbs,"any")
write.table(perc_sig, file = here("results","CONUS","pivot_pts","percent_sig_table.txt"), sep = "& ",eol = " \\\\ ", row.names = FALSE, col.names = FALSE)

for(wb in wbs){
    print(wb)
    wb_sig<-subset(sig_pvals,paste0("EVI_",wb,"_pval"))
    perc_sig<-cbind(perc_sig,zonal(wb_sig,  land_class, na.rm = TRUE)[,2])
    print(global(wb_sig, mean, na.rm = TRUE))
}



#filter out developed land, open water, ...
veg_mask<-land_class$`NLCD Land Cover Class`>40 & (land_class$`NLCD Land Cover Class`<80 | land_class$`NLCD Land Cover Class`>89 ) 
any_sig_wb_veg<-mask(any_sig_wb, veg_mask, maskvalues=FALSE, updatevalue=NA)
global(any_sig_wb_veg,  mean, na.rm = TRUE) #  0.7787224
