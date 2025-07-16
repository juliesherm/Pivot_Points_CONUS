RF_future_plot<-function(RF, rcp, era, quant,
                         pj_shade = TRUE){
  fut_evis<-list.files(here("results",park,"RF_future_fit"),
                       pattern = regex(paste0(RF,".+",rcp,".tif")),full.names = TRUE)
  fut<-lapply(fut_evis,rast)
  ## mask out  cultivate land (82), water (11), and barren (31)
  land_class_proj<-terra::rast(here("raw_data","nlcd_2021","nlcd_2021_proj","land_cover_2021_proj.tif"),
                               win = park_extent)
  veg_mask<-land_class_proj$`NLCD Land Cover Class`>31 & (land_class_proj$`NLCD Land Cover Class`<82) 
  
  fut<-lapply(fut, function(x) mask(x, veg_mask,maskvalues=FALSE,updatevalue=NA))
  ## turn into half-centuries
  fut[[1]]<-subset(fut[[1]], as.character(2030:2099))
  halfcent_mag<-lapply(fut, function(x) tapp(x,c(rep(TRUE, 35),rep(FALSE, 35)),sum))
  era_mag<-switch(era,
                  midcent = rast(lapply(halfcent_mag, function(x) subset(x, 1)/35)),
                  endcent = rast(lapply(halfcent_mag, function(x) subset(x, 2)/35)))
  #calculate quantile over  GCMS
  era_diff<-quantile(era_mag, quant ,na.rm = TRUE)/vi_mean
  names(era_diff)<-"quant"
  eranm<-switch(era,
                midcent="2030-2064",
                endcent="2065-2099",
                none = "")
  qnm<-switch(as.character(quant),
              "0"="minimum",
              "1"="maximum",
              "0.5" ="median",
              quant)
  if(pj_shade){
    #############################
    # pj  emphasized plot
    nps_class<-terra::vect(here("raw_data","vegmaps","blcageodata",
                                "blcageodata.gdb"),
                           layer = "fcl_Veg_Polys")
    raw_vi<-lapply(list.files(here("raw_data","MODIS","CONUS","VI_16Days_250m_v61","EVI","tiled"),
                              pattern = paste0("EVI","_2000_*"), full.names = T),
                   terra::rast)
    
    nps_class <- terra::project(nps_class,crs(raw_vi[[1]]))
    rm(raw_vi)
    nps_class<-terra::intersect(nps_class, park_extent)
    pj_classes<-subset(nps_class,  nps_class$Map_Unit_ID %in% 
                         c(69,71,72,74,78,87))
    
    era_diff_pj_masked<-mask(era_diff, pj_classes)
    
    plt<-ggplot()+
      geom_spatraster(data = era_diff_pj_masked, aes(fill = quant),
                      alpha = 1)+
      geom_spatraster(data = era_diff, aes(fill = quant), alpha = 0.4)+
      guides(alpha = "none") 
      
    
  }else{
    plt<-ggplot()+geom_spatraster(aes(fill = quant),data = era_diff)
    
  }
  # add color scheme and value polys
  plt<-plt+scale_fill_gradient2(low="burlywood4", mid="white", high="forestgreen", midpoint = 1, 
                                name = "Relative veg",oob = squish,
                                limit = c(0.5,1.5),
                                breaks = seq(0.5, 1.5, length.out = 5), 
                                labels = c("<0.5", "0.75", "1", "1.25", ">1.5")) +
    geom_spatvector(data = value_polys, aes(col = ValueRate),show.legend = FALSE, fill = NA,lwd = 1) + 
    scale_color_manual(values=c("orange", "red", "#56B4E9"))+
    geom_spatvector(data = park_boundary, fill = NA, col = "darkgray",  lwd = 1)+
    theme_void()+labs(subtitle = paste(rcp,qnm))
  
    
  return(plt)
}
