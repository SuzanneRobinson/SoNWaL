# calculate vulnerability

library(multidplyr)
outSpatHist<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\spatial\\historical\\SoNWal_spatOut_27_04_22.RDS")

hzYrsLoc<-"C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\hazardData\\hazYrs_V2_2051-2080.RDS"
outSpatLoc<-"C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\spatial\\rcp65\\SoNWal_spatOut_60_01.RDS"
outSpatHzLoc<-"C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\CHESSscape\\daily\\huss\\chess-scape_rcp60_bias-corrected_01_huss_uk_1km_daily_20800901-20800930.nc"
prHazLoc<-"C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\hazardData\\pHaz_2051-2080.rds"
  
plotVuln<-function(hzYrsLoc,outSpatLoc,outSpatHzLoc){
# read in hazard years
hzYrs<-do.call(rbind,readRDS(hzYrsLoc))%>%nest_by(x,y)
names(hzYrs)<-c("x","y","hzYrs")

hzPr<-readRDS(prHazLoc)
timePer<-as.numeric(sub('-.*', '', sub('.*hazYrs_V2_', '', hzYrsLoc)))-1

# read in map data
outSpat<-readRDS(outSpatLoc)


outSpatHz<-#as.data.frame(raster::coordinates(raster::brick(outSpatHzLoc))) %>%
  tibble(hzYrs) %>%
  left_join(outSpat) %>%
  mutate(hazPeriod=ifelse(Year-timePer > 0, Year - timePer, 0)) 


uncertainty_vals<-vuln_years%>%
  group_by(grid_id) %>%
  group_map(~uqFunc(.x$goodNPP, .x$badNPP, .y)) %>%
  bind_rows()


outSpatHzTmp<-merge(vuln_years,uncertainty_vals, by.x = "grid_id",by.y = "grid_id")

outSpatHz<-outSpatHzTmp%>%
  group_by(grid_id) %>%
  summarise(goodYears = mean(as.numeric(goodNPP), na.rm=T), 
            badYears=mean(as.numeric(badNPP), na.rm=T), 
            x=first(as.numeric(x)), 
            y=first(as.numeric(y)), 
            s_R=first(s_R),
            s_V=first(s_V),
           YC = mean(as.numeric(yc_value), na.rm = T),
            prHz=first(as.numeric(prHz)))%>%
  mutate(vuln=(((goodYears * 7.14) - (badYears * 7.14))))%>%
  mutate(risk = vuln * prHz)



#outSpatHz<-outSpatHz%>%mutate(vuln=ifelse(YC <= 16, -8, vuln))
#outSpatHist2<-outSpatHist[is.na(outSpatHist$Year)==F,]
#outSpatHist2<-outSpatHist2%>%group_by(x,y)%>%summarise(yc_value=first(yc_value))
#outSpatHz<-outSpatHz%>%left_join(outSpatHist2,by=c("x"="x","y"="y"))

outSpatHzXMask<-na.omit(filter(outSpatHz,YC<=16))

#outSpatHz<-outSpatHz%>%
#  mutate(vulnBins=cut(vuln,breaks=c(-10,-5,-4,0,0.5,1,2,4), labels =c("Unsuitable","none","none","low", "med", "high", "very high")))
outSpatHz$suit<-ifelse(outSpatHz$YC>summary(outSpatHz$yc_value)[5],"Highly Suitable", "Suitable")
outSpatHz$suit<-ifelse(outSpatHz$YC<summary(outSpatHz$yc_value)[4], "Unsuitable",outSpatHz$suit)


cols <- c("grey","#440154","#3b528b", "#21918c", "#5ec962","#fde725")

g1<-ggplot() +
  geom_tile(data = outSpatHz , aes(x = x, y = y, fill =  suit))+
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
    panel.background = element_rect(fill="white",color = NA),
    panel.grid.minor = element_blank())+
  coord_equal() + 
  scale_fill_viridis_d(expression(paste("Suitability",sep="")),option="A")+ 
  theme(legend.background = element_rect(fill = "white"),legend.text=element_text(color="black"),
        plot.title = element_text(color = "black"),
        legend.title = element_text(color="black"))+
  ggtitle(sub('\\..*', '', sub('.*hazYrs_V2_', '', hzYrsLoc)))


g2<-ggplot() +
  geom_tile(data = outSpatHz , aes(x = x, y = y, fill =  risk+(2*s_R)))+
 geom_tile(data=outSpatHzXMask,aes(x = x, y = y),fill = "grey")+
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
    panel.background = element_rect(fill="white",color = NA),
    panel.grid.minor = element_blank())+
  coord_equal() + 
  scale_fill_viridis_c(limits=c(-1.5,2),expression(paste("Risk",sep="")),na.value="white", option="turbo")+ 
  theme(legend.background = element_rect(fill = "white"),legend.text=element_text(color="black"),
        plot.title = element_text(color = "black"),
        legend.title = element_text(color="black"))+
  ggtitle(sub('\\..*', '', sub('.*hazYrs_V2_', '', hzYrsLoc)))

g3<-ggplot() +
  geom_tile(data = outSpatHz , aes(x = x, y = y, fill =  prHz))+
  geom_tile(data=outSpatHzXMask,aes(x = x, y = y),fill = "grey")+
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
    panel.background = element_rect(fill="white",color = NA),
    panel.grid.minor = element_blank())+
  coord_equal() + 
  scale_fill_viridis_c(limits=c(0,1),expression(paste("Probability of Hazard",sep="")),na.value="white", option="turbo")+ 
  theme(legend.background = element_rect(fill = "white"),legend.text=element_text(color="black"),
        plot.title = element_text(color = "black"),
        legend.title = element_text(color="black"))+
  ggtitle("2051-2080")+
  ggtitle(sub('\\..*', '', sub('.*hazYrs_V2_', '', hzYrsLoc)))


g4<-ggplot() +
  geom_tile(data = outSpatHz , aes(x = x, y = y, fill =  vuln))+
  geom_tile(data=outSpatHzXMask,aes(x = x, y = y),fill = "grey")+
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
    panel.background = element_rect(fill="white",color = NA),
    panel.grid.minor = element_blank())+
  coord_equal() + 
  scale_fill_viridis_c(limits=c(-2,3),expression(paste("Vulnerability",sep="")),na.value="white", option="turbo")+ 
  theme(legend.background = element_rect(fill = "white"),legend.text=element_text(color="black"),
        plot.title = element_text(color = "black"),
        legend.title = element_text(color="black"))+
  ggtitle(sub('\\..*', '', sub('.*hazYrs_V2_', '', hzYrsLoc)))


ggpubr::ggarrange(g1,g3,g4,g2)

}


plotVuln(hzYrsLoc,outSpatLoc,outSpatHzLoc)



outSpatHz<-as.data.frame(raster::coordinates(raster::brick(outSpatHzLoc))) %>%
  tibble(. ,hzYrs) %>%
  left_join(outSpat) %>%
  mutate(hazPeriod=ifelse(Year-timePer > 0, Year - timePer, 0)) %>%
  mutate(goodBad = unlist(pmap(outSpatHz[,c("hzYrs","hazPeriod","NPP_value")],checkFunc))) %>% 
  mutate(goodNPP=ifelse(goodBad=="good",NPP_value,NA),
         badNPP=ifelse(goodBad=="bad",NPP_value,NA),
         prHz=map_int(hzYrs,length)/30) %>%
  mutate(prHz=ifelse(lapply(hzYrs,function(x) is.na(x[1]))==T,0,prHz))%>%
  group_by(grid_id) %>%
  summarise(goodYears = mean(goodNPP, na.rm=T), 
            badYears=mean(badNPP, na.rm=T), 
            x=first(x), 
            y=first(y), 
            NPP=mean(NPP_value*7.14, na.rm = T), 
            YC = mean(yc_value, na.rm = T),
            prHz=first(prHz))%>%
  mutate(vuln=(((goodYears * 7.14) - (badYears * 7.14)))*hzPr) #%>%
# mutate(vuln=ifelse(YC <= 16, -8, vuln))



