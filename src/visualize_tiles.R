getwd()
library(here)
library(terra)

tiles<-lapply(list.files(here("results","CONUS","integrated_vi_gs","window_14days","cutoff_2.5","EVI","tiled"),
                          full.names = T,
                          pattern = "*.tif$"),terra::rast, lyrs = c(1,6,11,16))
together<-Reduce(merge,tiles)
together<-together/10000
mv<-global(together, min, na.rm = TRUE)
together<-terra::mask(together,  together, mask.values = mv)


pdf(here("results","CONUS","integrated_vi_gs","window_14days","cutoff_2.5","EVI",
         "evi_CONUS_plot.pdf"))
par(mfrow = c(2,2))
plot(together,1,  range = c(0,250), main = "EVI 2000")
plot(together,2,  range = c(0,250), main = "EVI 2005")
plot(together,3,  range = c(0,250), main = "EVI 2010")
plot(together,4,  range = c(0,250), main = "EVI 2020")

dev.off()


## plot pivot points
pivs<-lapply(list.files(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","tiled"),
                         full.names = T,
                         pattern = "*deficit"),terra::rast, lyrs = c(3,4,5))
names(pivs)<-list.files(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","tiled"),
                        full.names = F,
                        pattern = "*deficit")
wb_vars<-unique(stringr::str_split_i(names(pivs),"_",3))
vi_vars<-unique(stringr::str_split_i(names(pivs),"_",2))
vi_wb_prs <- expand.grid(vi_vars, wb_vars)

sum(stringr::str_detect(names(pivs),  paste0(x$Var1, "_",x$Var2)))

conus<-apply(vi_wb_prs,1, function(x){
  subpivs<-pivs[stringr::str_detect(names(pivs),  paste0(x[1], "_",x[2]))]
  return(Reduce(merge,subpivs))
}) 


pdf(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5",
         "pivs_deficit_plot.pdf"))
par(mfrow = c(3,2))
plot(conus[[1]],1,  range = c(-2,2), main = "EVI pivot point")
plot(conus[[2]],1,  range = c(-2,2), main = "SAVI pivot point")
plot(conus[[1]],2,  range = c(0,1), main = "EVI rsq")
plot(conus[[2]],2,  range = c(0,1), main = "SAVI rsq")
plot(conus[[1]],3,  range = c(0,1), main = "EVI pval")
plot(conus[[2]],3,  range = c(0,1), main = "SAVI pval")

dev.off()
