
###plot basfor outputs###

require(data.table) 
library(stringr)
fileLocs<-"/work/scratch-nopw/drcameron/outBASFOR"
filenames<-list.files(path = fileLocs, pattern = "\\.csv$", full.names = F, 
                      recursive = T)

##basfor combine function
readComb<-function(file){
  print(file)
  ff<-read.csv(file)
  fName<- substr(file,8,22)
  ff$fName<-sub(".csv","",file)
  return(ff)
}
ans = rbindlist(lapply(filenames, readComb),fill=T)
saveRDS(ans,"/home/users/aaronm7/BASFOR_spatOut_04_04.RDS")

library(data.table)
library(stringr)
#SoNWal combine function
fileLocs<-"/gws/nopw/j04/uknetzero/aarons_files/chess_manips/spatOutput_85_01_unc/"
filenames<-list.files(path = fileLocs, pattern = "\\.RDS$", full.names = F, 
                      recursive = T)
readComb<-function(file){
  print(file)
  ff<-readRDS(file)
  fName<- substr(file,8,22)
  ff$fName<-sub(".RDS","",fName)
  return(ff)
}
ans = rbindlist(lapply(filenames, readComb),fill=T)
saveRDS(ans,"/home/users/aaronm7/SoNWal_spatOut_85_01_unc.RDS")


###plot SoNWaL spatial outputs###

outSpatHist<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\spatial\\historical\\SoNWal_spatOut_27_04_22.RDS")
outSpat<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\spatial\\rcp26\\SoNWal_spatOut_26_01_unc.RDS")
#soilDataLocs<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_soil_data\\soilDataLocs.RDS")
#outSpat<-merge(lkList,outSpat,by.x="fName",by.y="fName")
library(grid)
outSpat<-outSpat[is.na(outSpat$Year)==F,]
#saveRDS(lkListSon,"lkListSon.RDS")
#lkListSon<-readRDS("lkListSon.RDS")
#outSpat$id<-lkListSon$fName
outSpat$GPPcsum <- ave(outSpat$GPP_value*365, outSpat$grid_id, FUN=cumsum)
outSpat$suit<-ifelse(outSpat$yc_value>summary(outSpatHist$yc_value)[5],"Highly Suitable", "Suitable")
outSpat$suit<-ifelse(outSpat$yc_value<summary(outSpatHist$yc_value)[4], "Unsuitable",outSpat$suit)
outSpat$suitU<-ifelse(outSpat$yc_q95>summary(outSpatHist$yc_value)[5],"Highly Suitable", "Suitable")
outSpat$suitU<-ifelse(outSpat$yc_q95<summary(outSpatHist$yc_value)[4], "Unsuitable",outSpat$suit)
outSpat$suitL<-ifelse(outSpat$yc_q05>summary(outSpatHist$yc_value)[5],"Highly Suitable", "Suitable")
outSpat$suitL<-ifelse(outSpat$yc_q05<summary(outSpatHist$yc_value)[4], "Unsuitable",outSpat$suit)

outSpat$YCrange<-outSpat$yc_q95-outSpat$yc_q05


years<-c(1961:2079)
ggList = as.list(years)
fg<-NULL
for(i in c(1:length(years))){
  outSpat<-as.data.frame(outSpat)
  outSpatY<-as.data.frame(outSpat%>%dplyr::filter(Year==years[i]|is.na(Year)==T))
  outSpatY<-outSpatY[!(is.na(outSpatY$Year) == T & is.na(outSpatY$Wsbr_q05) == F ), ]
 # fg2<-data.frame(Year=years[i],Total_GPP=sum(outSpatY$GPP_value*7.14,na.rm=T),Total_NPP=sum(outSpatY$NPP_value*7.14,na.rm=T),Total_NEE=sum(outSpatY$NEE_value*7.14,na.rm=T),AVG_LAI=mean(outSpatY$LAI_value,na.rm=T))
  
  #outSpatY$R<-ifelse(outSpatY$pH==0,0,outSpatY$R)
  
  #basOut3<-filter(basOut3,y>8e+05)
  #basOut3<-filter(basOut3,y<8.5e+05)
  
  g1<-ggplot() +
    geom_raster(data = outSpatY , aes(x = x, y = y, fill = ((R*7.14))))+
    theme_bw()+
    ggtitle("YC range")+
    theme(
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      panel.grid.major = element_blank(),
      plot.background = element_rect(fill="white", color = NA),
      panel.background = element_rect(fill="white", color = NA),
      panel.grid.minor = element_blank())+
    coord_equal() + 
    scale_fill_viridis_c(limits=c(-2,2),expression(paste("risk   ", sep="")), option = "turbo", na.value = "white")+ 
    theme(legend.background = element_rect(fill = "white"),legend.text=element_text(color="black",size=15),
          plot.title = element_text(color = "black", size=20),
          legend.title = element_text(color="black", size=20))
  
  g2<-ggplot() +
    geom_raster(data = outSpatY , aes(x = x, y = y, fill = ((LAI_value))))+
    theme_bw()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(),
          plot.background = element_rect(fill="white", color = NA),
          panel.background = element_rect(fill="white", color = NA),
          legend.key = element_rect(fill = "white"),
          panel.grid.minor = element_blank())+
    coord_equal() + 
    scale_fill_viridis_c(limits=c(0,12),expression(paste("LAI",sep="")),option = "turbo",na.value="white")+
    theme(legend.background = element_rect(fill = "white"),legend.text=element_text(color="black"),
          plot.title = element_text(color = "black"),
          legend.title = element_text(color="black"))
  
  
  g3<-ggplot() +
    geom_raster(data = outSpatY , aes(x = x, y = y, fill = ((NPP_q05*7.14))))+
    theme_bw()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(),
          plot.background = element_rect(fill="white", color = NA),
          panel.background = element_rect(fill="white", color = NA),
          panel.grid.minor = element_blank())+
    coord_equal() + 
    scale_fill_viridis_c(limits=c(0,15),expression(paste("NPP   ",sep="")),option = "turbo",na.value="white")+
    theme(legend.background = element_rect(fill = "white"),legend.text=element_text(color="black"),
          plot.title = element_text(color = "black"),
          legend.title = element_text(color="black"))
  
  g4<-ggplot() +
    geom_raster(data = outSpatY , aes(x = x, y = y, fill = ((NEE_value*7.14))))+
    theme_bw()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(),
          plot.background = element_rect(fill="white", color = NA),
          panel.background = element_rect(fill="white", color = NA),
          panel.grid.minor = element_blank())+
    coord_equal() + 
    scale_fill_viridis_c(limits=c(-8,8),expression(paste("NEE",sep="")),option = "turbo",na.value="white")+
    theme(legend.background = element_rect(fill = "white"),legend.text=element_text(color="black"),
          plot.title = element_text(color = "black"),
          legend.title = element_text(color="black"))
  
  f<-gridExtra::grid.arrange(g1,g2,g3,g4,top=textGrob(years[[i]],gp=gpar(fontsize=20,font=3)))
  ##4E5D6C
  g2<-cowplot::ggdraw(f) + 
    theme(plot.background = element_rect(fill="white"))
  # plot(g2)
  ggsave(paste0("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\markdownScrpts\\animations\\GPP_",years[i],".png"))
  # fg<-rbind(fg,fg2)
}
ggarrange(ggList[[1]],ggList[[2]],ggList[[3]],ggList[[4]],ggList[[5]],ggList[[6]],ggList[[7]],ggList[[8]],ggList[[9]],ggList[[10]],ggList[[11]],ggList[[12]],common.legend = T,legend="right")



imgs <- list.files("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\markdownScrpts\\animations\\rcp65_01", full.names = TRUE)
img_list <- lapply(imgs, image_read)
img_joined <- image_join(img_list)
img_animated <- image_animate(img_joined, fps = 2)
image_write(image = img_animated,
            path = "C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\gpp_growth.gif")



filenames<-unique(outSpat$fName)
readComb<-function(file){
  ff<-filter(outSpat,fName==file)
  fName<- file#substr(file,7,20)
  ff$fName<-paste0(fName,"_",rep(1:(nrow(ff)/46),each=46))
  return(ff)
}
ans = rbindlist(lapply(filenames, readComb))
lkList<-data.frame(fName=ans$fName,x=ans$x,y=ans$y)



