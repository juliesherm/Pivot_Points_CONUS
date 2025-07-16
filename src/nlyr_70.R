library(here)
library(terra)
fils<-list.files(here("results","CONUS","wb_proj","futures","tiled"), full.names = TRUE)
nms<-list.files(here("results","CONUS","wb_proj","futures","tiled"), full.names = FALSE)
rasts<-lapply(fils, rast)
for(i in 1:length(fils)){
  if(nlyr(rasts[[i]])==94){
    print(nms[i])
    subset(rasts[[i]],  25:96, overwrite = TRUE, filename  = here("results","CONUS","wb_proj","futures",paste0("new",nms[i])))
    
  }
}