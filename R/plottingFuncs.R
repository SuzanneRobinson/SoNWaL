###plotting functions for FR3pgn


## Function to plot the model against the Harwood data.
plotResults <- function(df,ShortTS=F){
  
  ##if plyr is loaded before dplyr can cause problems with groupby
  dt=12
  
  df <- df[c(2:nrow(df)),]
  df <- df %>% dplyr::group_by(Year)%>%mutate(cumGPP = cumsum(GPP),
                                              cumNPP = cumsum(NPP),
                                              timestamp = as.POSIXct(paste(sprintf("%02d",Year),sprintf("%02d",Month),sprintf("%02d",1),sep="-"),tz="GMT")) 
  
  
  if(ShortTS==T){
    df2<-NULL
    for(i in c(1:(nrow(df)-1))){
      
      if(df$Month[i]!=df$Month[i+1])
        df2<-rbind(df2,df[i,])
    }
    df2<-df2[-1,]
    df2$GPP<-aggregate(df$GPP~ df$Month+df$Year,FUN=sum)[-nrow(df2),3]
    df2$NPP<-aggregate(df$NPP~ df$Month+df$Year,FUN=sum)[-nrow(df2),3]
    df2$EvapTransp<-aggregate(df$EvapTransp~ df$Month+df$Year,FUN=sum)[-nrow(df2),3]
    df2$NEE<-aggregate(df$NEE~ df$Month+df$Year,FUN=sum)[-nrow(df2),3]
    df2$Reco<-aggregate(df$Reco~ df$Month+df$Year,FUN=sum)[-nrow(df2),3]
    df2$Rs<-aggregate(df$Rs~ df$Month+df$Year,FUN=sum)[-nrow(df2),3]
    df<-df2
  }  
  
  gpp1<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=timestamp,y=cumGPP,group=Year),colour="black",size=1)+
    geom_point(data=data,aes(x=timestamp,y=cumGPP),colour="red",size=2)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+
    labs(x="Year",y=expression(paste("GPP [tDM"," ",ha^-1,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  gpp2<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=timestamp,y=aggregate(df$GPP~ df$Month+df$Year,FUN=sum)[,3],group=Year),colour="black",size=1)+
    geom_point(data=data,aes(x=timestamp,y=gpp),colour="red",size=2)+
    ## geom_ribbon(data=data,aes(x=timestamp,ymin=gpp-gpp.sd,ymax=gpp+gpp.sd),alpha=0.3)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("GPP [tDM"," ",ha^-1,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  npp1<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=timestamp,y=cumNPP,group=Year),colour="black",size=1)+
    geom_point(data=data,aes(x=timestamp,y=cumNPP),colour="red",size=2)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+
    labs(x="Year",y=expression(paste("NPP [tDM"," ",ha^-1,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  npp2<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=timestamp,y=NPP,group=Year),colour="black",size=1)+
    geom_point(data=data,aes(x=timestamp,y=npp),colour="red",size=2)+
    ## geom_ribbon(data=data,aes(x=timestamp,ymin=npp-npp.sd,ymax=npp+npp.sd),alpha=0.3)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+
    labs(x="Year",y=expression(paste("NPP [tDM"," ",ha^-1,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  et<-ggplot()+theme_bw()+
    geom_point(data=data,aes(timestamp,et),size=2,colour="red")+
    geom_line(data=df,aes(timestamp,EvapTransp))+
    ## geom_ribbon(data=data,aes(x=timestamp,ymin=et-et.sd,ymax=et+et.sd),alpha=0.3)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+
    labs(x="Year",y=expression(paste(E[t]," [mm]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  nee<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=timestamp,y=NEE,group=Year),colour="black",size=1)+
    geom_point(data=data,aes(x=timestamp,y=nee),colour="red",size=2)+
    ## geom_ribbon(data=data,aes(x=timestamp,ymin=nee-nee.sd,ymax=nee+nee.sd),alpha=0.3)+
    geom_hline(yintercept=0,lty=3)+
    scale_y_continuous(limits=c(-5,2.5))+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("NEE [tDM"," ",ha^-1,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  reco<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=timestamp,y=Reco,group=Year),colour="black",size=1)+
    geom_point(data=data,aes(x=timestamp,y=reco),colour="red",size=2)+
    ## geom_ribbon(data=data,aes(x=timestamp,ymin=reco-reco.sd,ymax=reco+reco.sd),alpha=0.3)+
    geom_hline(yintercept=0,lty=3)+
    scale_y_continuous(limits=c(0,5))+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("Reco [tDM"," ",ha^-1,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  
  rs<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=timestamp,y=Rs,group=Year),colour="black",size=1)+
    geom_point(data=data,aes(x=timestamp,y=rs),colour="red",size=2)+
    ## geom_ribbon(data=data,aes(x=timestamp,ymin=rs-rs.sd,ymax=rs+rs.sd),alpha=0.3)+
    geom_hline(yintercept=0,lty=3)+
    scale_y_continuous(limits=c(0,2))+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("Rs [tDM"," ",ha^-1,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  
  leaf<-data.frame(rbind(
    c(dplyr::filter(df,Year==2015&Month==8)$Year+dplyr::filter(df,Year==2015&Month==8)$Month/dt,dplyr::filter(df,Year==2015&Month==8)$t.proj,5.7),
    c(dplyr::filter(df,Year==2018&Month==8)$Year+dplyr::filter(df,Year==2018&Month==8)$Month/dt,dplyr::filter(df,Year==2018&Month==8)$t.proj,5.56)
  )
  )
  names(leaf)<-c("year","age","lai")
  
  stemNo<-data.frame(rbind(
    c(dplyr::filter(df,Year==2018&Month==8)$Year+dplyr::filter(df,Year==2018&Month==8)$Month/dt,dplyr::filter(df,Year==2018&Month==8)$t.proj,1348)
  )
  )
  names(stemNo)<-c("year","age","stem")
  
  DBH<-data.frame(rbind(
    c(dplyr::filter(df,Year==2018&Month==8)$Year+dplyr::filter(df,Year==2018&Month==8)$Month/dt,dplyr::filter(df,Year==2018&Month==8)$t.proj,24.1)
  )
  )
  names(DBH)<-c("year","age","dbh")
  
  Wr<-data.frame(rbind(
    c(dplyr::filter(df,Year==2015&Month==8)$Year+dplyr::filter(df,Year==2015&Month==8)$Month/dt,dplyr::filter(df,Year==2015&Month==8)$t.proj,4.88)
  )
  )
  names(Wr)<-c("year","age","wr")
  
  difRoots<-data.frame(rbind(
    c(dplyr::filter(df,Year==2015&Month==8)$Year+dplyr::filter(df,Year==2015&Month==8)$Month/dt,dplyr::filter(df,Year==2015&Month==8)$t.proj,0.53)
  )
  )
  names(difRoots)<-c("year","age","difroots")
  
  totC<-data.frame(rbind(
    c(dplyr::filter(df,Year==2018&Month==7)$Year+dplyr::filter(df,Year==2018&Month==7)$Month/dt,dplyr::filter(df,Year==2018&Month==7)$t.proj,429.52)
  )
  )
  names(totC)<-c("year","age","totC")
  
  totN<-data.frame(rbind(
    c(dplyr::filter(df,Year==2018&Month==7)$Year+dplyr::filter(df,Year==2018&Month==7)$Month/dt,dplyr::filter(df,Year==2018&Month==7)$t.proj,14.30)
  )
  )
  names(totN)<-c("year","age","totN")
  
  pLAI<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=Year+Month/dt,y=LAI))+
    geom_point(data=leaf,aes(x=year,y=lai),shape=16,size=3,colour="red")+
    labs(x="Year",y=expression(paste("L"," ","[",m^-2," ",m^-2,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  pStemNo<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=Year+Month/dt,y=N))+
    geom_point(data=stemNo,aes(x=year,y=stem),shape=16,size=3,colour="red")+
    labs(x="Year",y=expression(paste("N"," ","[",ha^-1,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  pDBH<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=Year+Month/dt,y=dg))+
    geom_point(data=DBH,aes(x=year,y=dbh),shape=16,size=3,colour="red")+
    labs(x="Year",y="DBH [cm]")+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  pWr<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=Year+Month/dt,y=Wr))+
    geom_point(data=Wr,aes(x=year,y=wr),shape=16,size=3,colour="red")+
    labs(x="Year",y="Wr [tDM/ha]")+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  pdifRoots<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=Year+Month/dt,y=difRoots))+
    geom_point(data=difRoots,aes(x=year,y=difroots),shape=16,size=3,colour="red")+
    labs(x="Year",y="difRoots [tDM/ha]")+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))    
  
  ptotC<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=Year+Month/dt,y=totC))+
    geom_point(data=totC,aes(x=year,y=totC),shape=16,size=3,colour="red")+
    labs(x="Year",y="totC [t/ha]")+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  ptotN<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=Year+Month/dt,y=totN))+
    geom_point(data=totN,aes(x=year,y=totN),shape=16,size=3,colour="red")+
    labs(x="Year",y="totN [t/ha]")+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  return(list(gpp1,gpp2,npp1,npp2,nee,reco,rs,pLAI,pStemNo,pDBH,pWr,et,ptotC,ptotN))
}

## Function to plot the key output pools and nutritional modifier of the model
plotModel<-function(output.df){
  
  cls<-iprior::gg_colour_hue(2)
  
  flux<-ggplot(data=output.df)+theme_bw()+geom_line(aes(x=t.proj,y=GPP,colour=cls[1]))+
    geom_line(aes(x=t.proj,y=NPP,colour=cls[2]))+scale_colour_manual(values=cls,name="Fluxes",labels=c("NPP","GPP"))+theme(panel.grid=element_blank())+labs(x="Year",y="Flux")
  
  
  lai<-ggplot(data=output.df)+theme_bw()+geom_line(aes(x=t.proj,y=LAI,colour=cls[1]))+
    scale_colour_manual(values=cls,name="Fluxes",labels=c("LAI"))+theme(panel.grid=element_blank())+labs(x="Year",y="Flux")
  
  
  g<-ggplot(data=output.df)+theme_bw()+geom_line(aes(x=t.proj,y=G,colour=cls[1]))+
    scale_colour_manual(values=cls,name="Fluxes",labels=c("G"))+theme(panel.grid=element_blank())+labs(x="Year",y="G")
  
  
  dg<-ggplot(data=output.df)+theme_bw()+geom_line(aes(x=t.proj,y=dg,colour=cls[1]))+
    scale_colour_manual(values=cls,name="Fluxes",labels=c("DG"))+theme(panel.grid=element_blank())+labs(x="Year",y="Dg")
  
  
  hdom<-ggplot(data=output.df)+theme_bw()+geom_line(aes(x=t.proj,y=hdom,colour=cls[1]))+
    scale_colour_manual(values=cls,name="Fluxes",labels=c("HDOM"))+theme(panel.grid=element_blank())+labs(x="Year",y="H")
  
  
  rad<-ggplot(data=output.df)+theme_bw()+geom_line(aes(x=t.proj,y=RAD,colour=cls[1]))+
    scale_colour_manual(values=cls,name="Fluxes",labels=c("RAD"))+theme(panel.grid=element_blank())+labs(x="Year",y="Rg")
  
  fN<-ggplot(data=output.df)+theme_bw()+geom_line(aes(x=t.proj,y=fN,colour=cls[1]))+
    scale_colour_manual(values=cls,name="Fluxes",labels=c("fNn"))+theme(panel.grid=element_blank())+labs(x="Year",y="fNn")
  
  nav<-ggplot(data=output.df)+theme_bw()+geom_line(aes(x=t.proj,y=Nav,colour=cls[1]))+
    scale_colour_manual(values=cls,name="Fluxes",labels=c("Nav"))+theme(panel.grid=element_blank())+labs(x="Year",y="Nav")
  
  
  totc<-ggplot(data=output.df)+theme_bw()+geom_line(aes(x=t.proj,y=totC,colour=cls[1]))+
    scale_colour_manual(values=cls,name="Fluxes",labels=c("totC"))+theme(panel.grid=element_blank())+labs(x="Year",y="TotC")
  
  
  totn<-ggplot(data=output.df)+theme_bw()+geom_line(aes(x=t.proj,y=totN,colour=cls[1]))+
    scale_colour_manual(values=cls,name="Fluxes",labels=c("totN"))+theme(panel.grid=element_blank())+labs(x="Year",y="TotN")
  
  
  return(list(flux,lai,rad,g,dg,hdom,fN,nav,totc,totn))
}

## Function to find the yield class of the stand based on stand dominant height (assuming it is the same as top height)
YC.finder <- function(HT,AGE) 
{
  YC.site = ((HT/(1-exp(-0.033329*AGE))^1.821054)-14.856317)/1.425397
  if(YC.site>24) 24
  else if (YC.site<4) 4
  else unlist(sub("\\]","",unlist(strsplit(as.character(cut(YC.site,breaks=c(6,8,10,12,14,16,18,20,22,24),right=T)),split="\\,"))[2]))
}

############################################################################################
