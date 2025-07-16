#########################################
## VI processing
#########################################
#read in
vi_stk<-read_in_data(park,"MODIS")
cloud_stk<-read_in_data(park,"Rely")
vi_vars<-names(vi_stk)

#mask
vi_stk<-lapply(vi_vars, function(x){
  st<-Sys.time()
  vi_mask<-terra::mask(vi_stk[[x]], cloud_stk,  maskvalues = 2, updatevalue = 500)
  vi_mask<-terra::mask(vi_mask, cloud_stk,  maskvalues = 3, updatevalue = NA)  
  print(Sys.time()-st)
  return(vi_mask)
}) 
names(vi_stk)<-vi_vars

#smooth
source(here("src","smooth_ts.R"))
vi_smth<-terra::app(vi_stk[["EVI"]], function(x){
  smooth_ts(x, dates = vi_dates)
}, filename = here("results",park,"vi_smoothed","EVI_smth.tif"),
overwrite = TRUE)

vi_smth<-terra::app(vi_stk[["SAVI"]], function(x){
  smooth_ts(x, dates = vi_dates)
}, filename = here("results",park,"vi_smoothed","SAVI_smth.tif"),
overwrite = TRUE)

vi_stk<-lapply(vi_vars,function(x) terra::rast(here("results",park,"vi_smoothed",paste0(x,"_smth.tif"))))
names(vi_stk)<-vi_vars

#growing season summary
source(here("src","integrate_gs.R"))
gs_dates<-terra::rast(here::here("results","CONUS","gs_dates","gs_all_doys_proj.tif"))
gs_dates<-terra::crop(gs_dates, vi_stk[[1]])
yrs<-2000:2021
EVI_stk<-vi_stk[["EVI"]]
lapply(yrs, function(yr){
  gs_rast<-c(gs_dates[paste0("first",yr)],gs_dates[paste0("last",yr)],
             terra::subset(EVI_stk,yr == year(vi_dates)))
  terra::app(gs_rast,integrate_gs,ydays = yday(vi_dates[yr == year(vi_dates)]),
             filename = here("results",park,"vi_gs_summary",paste0("EVI",yr,".tif")))
})


SAVI_stk<-vi_stk[["SAVI"]]
lapply(yrs, function(yr){
  gs_rast<-c(gs_dates[paste0("first",yr)],gs_dates[paste0("last",yr)],
             terra::subset(EVI_stk,yr == year(vi_dates)))
  terra::app(gs_rast,integrate_gs,ydays = yday(vi_dates[yr == year(vi_dates)]),
             filename = here("results",park,"vi_gs_summary",paste0("SAVI",yr,".tif")))
})

vi_stk<-lapply(vi_vars,function(x) {
  terra::rast(list.files(here("results",park,"vi_gs_summary"), pattern = x, 
                               full.names = TRUE))
})

#scale and center
vi_gs_anom <- lapply(vi_stk, terra::app, scale, center = TRUE, scale = FALSE)
vi_gs_anom <- lapply(vi_stk, function(x) x/10000)
names(vi_gs_anom) <- vi_vars

