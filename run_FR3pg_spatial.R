library("tidyr")
library("tidyverse")
library("purrr")
library("BayesianTools")
library("sensitivity")
library("dplyr")
library("future")
library("ncdf4)" )
library("raster)" )
library("ggplot2")
library("httr")
library("furrr")
library("viridis")
library("fr3PGDN")
library("tibble")
library("miscTools")
library("parallel")
library("sf")
library("rgdal")
library(raster)


data_dir="/gws/nopw/j04/hydro_jules/data/uk/driving_data/chess/chess-met/daily"#Location of input data
output_dir="/home/users/aaronm7/"#output file for bricked rasters
save_file="/home/users/aaronm7/spatialChunks/"#save file for spatial chunks of data

output_dir="C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\"
data_dir="C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\Prafor\\models\\spatial_met_data\\CHESSscape\\daily"#Location of input data
save_file="C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\Prafor\\spatialChunks\\"#save file for spatial chunks of data




spatSplitXScape(data_dir=data_dir,output_dir=output_dir,save_file=save_file,dates=1960:2100,numChunks=c(1:31))
#create chunks of approx 10,000 grid cells in parallel, saves to file
#1:35 is approx scotland
plan(multisession, workers = 8)

future_map(c(31), ~spat_dat_UKnc(chunk=.x,output_dir,save_file),.progress = T)





args=(commandArgs(TRUE))
print(args)
chunk=args[1]

fileLocs<-"C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\spatialChunks\\"
paramsFile<-("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\weekly_1_T.RDS")
files <- list.files(path = fileLocs, pattern = "\\.RDS$", full.names = TRUE, 
                    recursive = T)

fileName<-sub("\\/.*", "",list.files(path = fileLocs, pattern = "\\.RDS$", full.names = F, 
                           recursive = T))[chunk]
fileName<-str_sub(fileName,14)

simDat<-readRDS(files[chunk])
clm_df_full<<-getClimDat("weekly")

#Read in MCMC output to get parameter draws - this code needs tidying a bit
param_drawX<-readRDS(paramsFile)
param_drawX<-as_tibble(getSample(param_drawX, start = 100, coda = TRUE, thin = 1,numSamples = 25 )[[1]])#,colMedians(getSample(param_drawX, start = 1000, coda = TRUE, thin = 1 )[[2]])))
param_draw<-as_tibble(1:nrow(param_drawX))
param_draw$pars<- tibble(list(split(param_drawX, 1:NROW(param_drawX))))%>%
  unnest_legacy()
names(param_draw)<-c("mcmc_id","pars")
param_draw<<-param_draw
param_draw$pars<-unname(param_draw$pars$`list(split(param_drawX, 1:NROW(param_drawX)))`)
simDat$site<-list(data.frame(from=(paste0( round(1971),"-01-01")),to=paste0( round(2016),"-30-12")))
soilLocVals<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_soil_data\\soilDataLocs.RDS")
#names(soilLocVals)<-c("vals","latitude","longitude","x","y","soil","soilDepth")



simDat<-as_tibble(merge(simDat,soilLocVals,by.x=c("x","y"),by.y=c("x","y")))
#run grid squares in parallel - as running over grid squares should be very scaleable
#outTemp <-pmap(simDat$site, simDat$clm,as.list(simDat$grid_id), .f=~FR3PG_spat_run(site=..1, clm=..2,grid_id=..3,param_draw=param_draw),.progress = T)
simDat$wp[is.na(simDat$wp)==T]<-0
simDat$fc[is.na(simDat$fc)==T]<-0
simDat$sp[is.na(simDat$sp)==T]<-0
simDat$soilCond[is.na(simDat$soilCond)==T]<-0
simDat$soilTex[is.na(simDat$soilTex)==T]<-0
simDat$soilDepth[is.na(simDat$soilDepth)==T]<-0

outTemp<-mapply(FR3PG_spat_run, site = simDat$site, clm = simDat$clm,soil=simDat$soilTex,wp=simDat$wp,fc=simDat$fc,sp=simDat$sp,cond=simDat$soilCond,soilDepth=simDat$soilDepth,grid_id=as.list(simDat$grid_id),MoreArgs = list(param_draw=param_draw),SIMPLIFY = F)

#outTemp<-do.call(rbind,outTemp)

#bind into a single tibble
out<-as_tibble(data.table::rbindlist(outTemp,fill=T))
#re-add grid_id values
#out$grid_id<-simDat$grid_id
grF<-simDat[,c(1,2,3,7,8,9)]
out<-merge(out,grF,by.x="grid_id",by.y = "grid_id")

saveRDS(out,paste0("/work/scratch-nopw/alm/spatOutput/","SoNWal_",fileName))










#function to find yield class - this has been changed - see updated function in shiny code etc!
YC.finder <- function(HT,AGE=59) 
{
  if(is.na(HT)==F){
  YC.site = ((HT/(1-exp(-0.033329*AGE))^1.821054)-14.856317)/1.425397
  if(YC.site>24) 24
  else if (YC.site<4) 4
  else unlist(sub("\\]","",unlist(strsplit(as.character(cut(YC.site,breaks=c(6,8,10,12,14,16,18,20,22,24),right=T)),split="\\,"))[2]))
} 
else return (NA)
  }

out$yc_value<-(out$yc_value-8)/(24-8)
out$yc_q95<-(out$yc_q95-8)/(24-8)
out$yc_q05<-(out$yc_q05-8)/(24-8)
#
#
out$yc_value<-ifelse(out$yc_value<0,0,out$yc_value)
out$yc_q95<-ifelse(out$yc_q95<0,0,out$yc_q95)
out$yc_q05<-ifelse(out$yc_q05<0,0,out$yc_q05)



#Plot results
g1<-out %>%
  mutate( range = yc_q95 - yc_q05) %>%
  dplyr::select( grid_id, mean = yc_value) %>%
  gather( variable, yc_value, -grid_id) %>%
  inner_join( ., simDat, by = 'grid_id') %>%
  ggplot( aes(x, y, fill = yc_value )) + #*100 to get hectares from 1km grid squares
  geom_raster()+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle("Suitability")+
  coord_equal() + 
 scale_fill_viridis_c( "",option = "plasma", na.value="gray70",limits=c(0,1),
                       breaks=c(0,0.5,1),labels=c("Unsuitable","Suitable","V. Suitable"))


g2<-out %>%
  mutate( range = GPP_q95 - GPP_q05) %>%
  dplyr::select( grid_id, mean = GPP_value) %>%
  gather( variable, GPP_value, -grid_id) %>%
  inner_join( ., simDat, by = 'grid_id') %>%
  ggplot( aes(x, y, fill = GPP_value )) + #*100 to get hectares from 1km grid squares
  geom_raster()+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle("GPP")+
  coord_equal() + 
  scale_fill_viridis_c( expression(paste("GPP [tDM"," ",ha^-1,"]",sep="")),option = "plasma", na.value="gray70")

g3<-out %>%
  mutate( range = NPP_q95 - NPP_q05) %>%
  dplyr::select( grid_id, mean = NPP_value) %>%
  gather( variable, NPP_value, -grid_id) %>%
  inner_join( ., simDat, by = 'grid_id') %>%
  ggplot( aes(x, y, fill = NPP_value )) + #*100 to get hectares from 1km grid squares
  geom_raster()+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle("NPP")+
  coord_equal() + 
  scale_fill_viridis_c( expression(paste("NPP [tDM"," ",ha^-1,"]",sep="")),option = "plasma", na.value="gray70",
                        limits=c(1,4))

g4<-out %>%
  mutate( range = NEE_q95 - NEE_q05) %>%
  dplyr::select( grid_id, mean = NEE_value) %>%
  gather( variable, NEE_value, -grid_id) %>%
  inner_join( ., simDat, by = 'grid_id') %>%
  ggplot( aes(x, y, fill = NEE_value )) + #*100 to get hectares from 1km grid squares
  geom_raster()+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle("NEE")+
  coord_equal() + 
  scale_fill_viridis_c( expression(paste("NEE [tDM"," ",ha^-1,"]",sep="")),option = "plasma", na.value="gray70")
  
  
  
  
  ggarrange(g1,g2,g3,g4)
  
  
  
  
  

    outSpatY<-outSpat%>%filter(Year==year|is.na(Year)==T)
    #outSpatY<-outSpatY[duplicated(outSpatY$grid_id),]
    
    outSpatY$Chunk<-substr(outSpatY$fName,1,2)
    outSpatY$Chunk<-sub("_", "",  outSpatY$Chunk)
    outSpatY$Chunk[is.na(outSpatY$Wsbr_q05)]<-NA
    
    
    g1<-  ggplot() +
        geom_raster(data = outSpatY , aes(x = x, y = y, fill = as.numeric(Chunk)))+
        theme_bw()+
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())+
        ggtitle("Big chunks")+
        coord_equal() + 
        scale_fill_viridis_c("Chunk",option = "plasma")
    
  
    

    outSpatX<-filter(outSpatY,Chunk>(40))
    outSpatX<-filter(outSpatX,Chunk<(49))
    
    outSpatX$SmallChunk<-substr(outSpatX$fName,4,6)
    #fg<-data.frame(sc=unique(outSpatX$SmallChunk))
    #fg$numVal<-1:nrow(fg)
   # outSpatX2<-merge(outSpatX,fg,by.x="SmallChunk",by.y="sc")
    
    
 g2<-  ggplot() +
        geom_raster(data = outSpatX , aes(x = x, y = y, fill = as.numeric(SmallChunk)))+
        theme_bw()+
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())+
        ggtitle("Small chunks")+
        coord_equal() + 
        scale_fill_viridis_c("Smaller split",option = "plasma", na.value="gray70")

ggarrange(g1,g2,common.legend = T,legend="bottom")



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
saveRDS(ans,"/home/users/aaronm7/BASFOR_spatOut_28_10.RDS")

library(data.table)
#SoNWal combine function
fileLocs<-"/work/scratch-nopw/alm/spatOutput"
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
saveRDS(ans,"/home/users/aaronm7/SoNWal_spatOut_16_11.RDS")








outSpat<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\SoNWal_spatOut_26_10.RDS")
#soilDataLocs<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_soil_data\\soilDataLocs.RDS")
#outSpat<-merge(lkList,outSpat,by.x="fName",by.y="fName")
library(grid)
outSpat<-outSpat[is.na(outSpat$Year)==F,]
#saveRDS(lkListSon,"lkListSon.RDS")
#lkListSon<-readRDS("lkListSon.RDS")
outSpat$id<-lkListSon$fName
outSpat$GPPcsum <- ave(outSpat$GPP_value*365, outSpat$id, FUN=cumsum)
years<-c(1971:2016)
ggList = as.list(years)
fg<-NULL
for(i in c(1:length(years))){
  
  outSpatY<-as.data.frame(outSpat%>%filter(Year==years[i]|is.na(Year)==T))
  outSpatY<-outSpatY[!(is.na(outSpatY$Year) == T & is.na(outSpatY$Wsbr_q05) == F ), ]
  fg2<-data.frame(Year=years[i],Total_GPP=sum(outSpatY$GPP_value*7.14,na.rm=T),Total_NPP=sum(outSpatY$NPP_value*7.14,na.rm=T),Total_NEE=sum(outSpatY$NEE_value*7.14,na.rm=T),AVG_LAI=mean(outSpatY$LAI_value,na.rm=T))
  #basOut3<-filter(basOut3,y>8e+05)
  #basOut3<-filter(basOut3,y<8.5e+05)
  
  g1<-ggplot() +
  geom_raster(data = outSpatY , aes(x = x, y = y, fill = ((GPP_value*7.14))))+
  theme_bw()+
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
  scale_fill_viridis_c(limits=c(0,12),expression(paste("GPP      ",sep="")),option = "turbo",na.value="white")+ 
    theme(legend.background = element_rect(fill = "white"),legend.text=element_text(color="black"),
          plot.title = element_text(color = "black"),
          legend.title = element_text(color="black"))
  
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
    scale_fill_viridis_c(limits=c(0,10),expression(paste("LAI",sep="")),option = "turbo",na.value="white")+
    theme(legend.background = element_rect(fill = "white"),legend.text=element_text(color="black"),
          plot.title = element_text(color = "black"),
          legend.title = element_text(color="black"))
  
  
  g3<-ggplot() +
    geom_raster(data = outSpatY , aes(x = x, y = y, fill = ((NPP_value*7.14))))+
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
    scale_fill_viridis_c(limits=c(0,6),expression(paste("NPP   ",sep="")),option = "turbo",na.value="white")+
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
    scale_fill_viridis_c(limits=c(-4,4),expression(paste("NEE",sep="")),option = "turbo",na.value="white")+
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



imgs <- list.files("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\markdownScrpts\\animations\\", full.names = TRUE)
img_list <- lapply(imgs, image_read)
img_joined <- image_join(img_list)
img_animated <- image_animate(img_joined, fps = 2)
image_write(image = img_animated,
            path = "C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\gpp_growth.gif")





basOut<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\BASFOR_spatOut_28_10.RDS")




filenames<-unique(outSpat$fName)
readComb<-function(file){
  ff<-filter(outSpat,fName==file)
  fName<- file#substr(file,7,20)
  ff$fName<-paste0(fName,"_",rep(1:(nrow(ff)/46),each=46))
  return(ff)
}
ans = rbindlist(lapply(filenames, readComb))
lkList<-data.frame(fName=ans$fName,x=ans$x,y=ans$y)







outSpat<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\SoNWal_spatOut_26_10.RDS")
outSpat$xy<-paste0(outSpat$x,outSpat$y)
filenames<-unique(outSpat$fName)
readComb<-function(file){
  ff<-filter(outSpat,fName==file)
  ff<-ff[!duplicated(ff$xy),]
  fName<- file#substr(file,7,20)
  ff$loc<-paste0(fName,"_",1:nrow(ff))
  return(ff)
}
#ans = rbindlist(lapply(filenames, readComb))
#lkList<-data.frame(loc=ans$loc,x=ans$x,y=ans$y)
lkList<-readRDS("locNames.RDS")

basOut<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\BASFOR_spatOut_28_10.RDS")
basOut$fName<-substr(basOut$fName,8,20)
basOut2<-basOut
basOut2$GPPcsumBas <- ave(basOut2$GPP_gCm2d*365, basOut2$fName, FUN=cumsum)

#lkListD<-lkList[duplicated(lkList$fName)==F,]
basOut2<-merge(lkList,basOut2,by.x="loc",by.y="fName",all=F)

years<-c(2005:2016)
ggList = as.list(years)
for(i in c(1:length(years))){
  
  basOut3<-filter(basOut2,Year==years[[i]])
  basOut3<-merge(outSpatY,basOut3,by.x=c("x","y"),by.y=c("x","y"),all=T)
  

  ggList[[i]]<-ggplot() +
  geom_raster(data = basOut3 , aes(x = x, y = y, fill = as.numeric(GPPcsumBas)))+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle("Cumulative difference BASFOR-SoNWal")+
  coord_equal() + 
  scale_fill_viridis_c(expression(paste("Cumulative GPP [gC"," ",cm^-2,"]",sep="")),option = "turbo")+
    theme(text = element_text(size = 20))

}

ggarrange(ggList[[1]],ggList[[2]],ggList[[3]],ggList[[4]],ggList[[5]],ggList[[6]],ggList[[7]],ggList[[8]],ggList[[9]],ggList[[10]],ggList[[11]],ggList[[12]],common.legend = T,legend="right")


##Sonwal
outSpatY<-as.data.frame(outSpat%>%filter(Year==1995|is.na(Year)==T))
outSpatY<-outSpatY[!(is.na(outSpatY$Year) == T & is.na(outSpatY$Wsbr_q05) == F ), ]
outSpatY<-outSpatY[is.na(outSpatY$Year)==F,]

ggSon<-ggplot() +
  geom_raster(data = outSpatY , aes(x = x, y = y, fill = (GPP_value*7.14)))+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle("SoNWal - 1995")+
  coord_equal() + 
  scale_fill_viridis_c(limits=c(0,12),expression(paste("GPP [gC"," ",cm^-2,"]",sep="")),option = "D")


##Soil plot
outSpatY$soil[outSpatY$soil==0]<-9
ggSoil<-ggplot() +
  geom_raster(data = outSpatY , aes(x = x, y = y, fill = as.factor(soil)))+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle("Soil type - corresponds to ESDB & soil hydraulic properties map")+
  coord_equal() + 
  scale_fill_viridis_d("Soil type",option = "D")

ggarrange(ggSon,ggBas,ggSonX,ggBasX,common.legend = T,legend="right")
ggarrange(cGPPSon+ggtitle("SoNWal"),ggBasCm+ggtitle("BASFOR"),common.legend = T,legend="right")




#plot model gpp averages
outSpatX<-outSpat%>%group_by(Year)%>%summarise(meanGPP=mean(GPP_value*7.14))
outSpatX$Model<-rep("SoNWal",nrow(outSpatX))
outSpatX2<-basOut%>%filter(Year<2017)%>%group_by(Year)%>%summarise(meanGPP=mean(GPP_gCm2d))
outSpatX2$Model<-rep("BASFOR",nrow(outSpatX2))
outSpatX<-rbind(outSpatX,outSpatX2)

ggplot() +
  geom_line(data = outSpatX , aes(x = Year, y = meanGPP, col = Model),size=1.2)+
  theme_bw()+
  ylab(expression(paste("GPP [gC"," ",cm^-2,day^-1,"]",sep="")))+
  scale_color_viridis_d(option="D")+
  theme(text = element_text(size = 20))




outSpatY<-as.data.frame(outSpat%>%filter(Year==1995|is.na(Year)==T))
outSpatY<-outSpatY[!(is.na(outSpatY$Year) == T & is.na(outSpatY$Wsbr_q05) == F ), ]
outSpatY<-outSpatY[is.na(outSpatY$Year)==F,]


ggplot() +
  geom_raster(data = outSpatY , aes(x = x, y = y, fill = as.factor(round((yc_value)))))+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle(years[i])+
  coord_equal() + 
  scale_fill_viridis_d(expression(paste("GPP [gC"," ",cm^-2,year^-1,"]",sep="")),option = "B")
