##########################################
## CUSTOM FUNCTIONS FOR PLOTING RESULTS ##
##########################################

## Function to plot the model against the Harwood data.
plotResults <- function(df){
  
  df <- df[c(2:nrow(df)),]
  df <- df %>% group_by(Year) %>% mutate(cumGPP = cumsum(GPP),
                                         cumNPP = cumsum(NPP),
                                         timestamp = as.POSIXct(paste(sprintf("%02d",Year),sprintf("%02d",Month),sprintf("%02d",1),sep="-"),tz="GMT")) 
  
  gpp1<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=timestamp,y=cumGPP,group=Year),colour="black",size=1)+
    geom_point(data=data,aes(x=timestamp,y=cumGPP),colour="red",size=2)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+
    labs(x="Year",y=expression(paste("GPP [tDM"," ",ha^-1,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  gpp2<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=timestamp,y=GPP,group=Year),colour="black",size=1)+
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
    c(dplyr::filter(df,Year==2015&Month==8)$Year+dplyr::filter(df,Year==2015&Month==8)$Month/12,dplyr::filter(df,Year==2015&Month==8)$t.proj,5.7),
    c(dplyr::filter(df,Year==2018&Month==8)$Year+dplyr::filter(df,Year==2018&Month==8)$Month/12,dplyr::filter(df,Year==2018&Month==8)$t.proj,5.56)
  )
  )
  names(leaf)<-c("year","age","lai")
  
  stemNo<-data.frame(rbind(
    c(dplyr::filter(df,Year==2018&Month==8)$Year+dplyr::filter(df,Year==2018&Month==8)$Month/12,dplyr::filter(df,Year==2018&Month==8)$t.proj,1348)
  )
  )
  names(stemNo)<-c("year","age","stem")
  
  DBH<-data.frame(rbind(
    c(dplyr::filter(df,Year==2018&Month==8)$Year+dplyr::filter(df,Year==2018&Month==8)$Month/12,dplyr::filter(df,Year==2018&Month==8)$t.proj,24.1)
  )
  )
  names(DBH)<-c("year","age","dbh")
  
  Wr<-data.frame(rbind(
    c(dplyr::filter(df,Year==2015&Month==8)$Year+dplyr::filter(df,Year==2015&Month==8)$Month/12,dplyr::filter(df,Year==2015&Month==8)$t.proj,4.88)
  )
  )
  names(Wr)<-c("year","age","wr")
  
  difRoots<-data.frame(rbind(
    c(dplyr::filter(df,Year==2015&Month==8)$Year+dplyr::filter(df,Year==2015&Month==8)$Month/12,dplyr::filter(df,Year==2015&Month==8)$t.proj,0.53)
  )
  )
  names(difRoots)<-c("year","age","difroots")
  
  totC<-data.frame(rbind(
    c(dplyr::filter(df,Year==2018&Month==7)$Year+dplyr::filter(df,Year==2018&Month==7)$Month/12,dplyr::filter(df,Year==2018&Month==7)$t.proj,86.7)
  )
  )
  names(totC)<-c("year","age","totC")
  
  totN<-data.frame(rbind(
    c(dplyr::filter(df,Year==2018&Month==7)$Year+dplyr::filter(df,Year==2018&Month==7)$Month/12,dplyr::filter(df,Year==2018&Month==7)$t.proj,2.16)
  )
  )
  names(totN)<-c("year","age","totN")
  
  pLAI<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=Year+Month/12,y=LAI))+
    geom_point(data=leaf,aes(x=year,y=lai),shape=16,size=3,colour="red")+
    labs(x="Year",y=expression(paste("L"," ","[",m^-2," ",m^-2,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  pStemNo<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=Year+Month/12,y=N))+
    geom_point(data=stemNo,aes(x=year,y=stem),shape=16,size=3,colour="red")+
    labs(x="Year",y=expression(paste("N"," ","[",ha^-1,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  pDBH<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=Year+Month/12,y=dg))+
    geom_point(data=DBH,aes(x=year,y=dbh),shape=16,size=3,colour="red")+
    labs(x="Year",y="DBH [cm]")+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  pWr<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=Year+Month/12,y=Wr))+
    geom_point(data=Wr,aes(x=year,y=wr),shape=16,size=3,colour="red")+
    labs(x="Year",y="Wr [tDM/ha]")+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  pdifRoots<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=Year+Month/12,y=difRoots))+
    geom_point(data=difRoots,aes(x=year,y=difroots),shape=16,size=3,colour="red")+
    labs(x="Year",y="difRoots [tDM/ha]")+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))    
  
  ptotC<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=Year+Month/12,y=totC))+
    geom_point(data=totC,aes(x=year,y=totC),shape=16,size=3,colour="red")+
    labs(x="Year",y="totC [t/ha]")+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  ptotN<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=Year+Month/12,y=totN))+
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

######################
## START THE SCRIPT ##
######################

## Load necessary packages
library(fr3PGDN,quietly=TRUE)
library(tidyverse,quietly=TRUE)

## Switch off annoying warnings
options(warn=-1)

## Years of data to use for calibration
startYear = 2015
endYear = 2018

## Met data
clm.df.full<-read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\3pg\\data\\clm_df_full.csv")

## Read Harwood data for Sitka spruce and mutate timestamp to POSIXct
data <- read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\3pg\\data\\harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))

#######################################################
## Parameters for Sitka spruce - Minunno et al. 2010 ##
#######################################################
sitka<-list(weather=clm.df.full,
            ## ~~ Initial pools ~~ ##
            Wl = 0.01,
            WlDormant = 0,
            Wr = 0.01,
            Wsbr = 0.1,
            Wlitt = 0,
            YrC = 42.95,
            YlC = 85.90,
            OC = 300.66,
            YrN = 1.43,
            YlN = 2.86,
            ON = 10.01,
            Nav = 8.76,
            ## ~~ Site ~~ ##
            N = 2000,
            rotation = 1,
            cycle = 1,
            rm.sprouts = F,
            nyears = 35,
            initial.month = 1,
            latitude = 57.06,
            soilclass = -1,
            ASW = 165,
            MaxASW = 400,
            MinASW = 0,
            CO2 = 400,
            ## ~~ Parameters ~~ ##
            pFS2 = 1.4,
            pFS20 = 0.8,
            aS = 0.138,
            nS = 2.3,
            pRx = 0.45,
            pRn = 0.3,
            Tmin = -5,
            Topt = 15,
            Tmax = 35,
            kF = 1,
            SWconst0 = 0.55,
            SWpower0 = 6,
            m0 = 0,
            MaxAge = 400,
            nAge = 4,
            rAge = 0.95,
            gammaFx = 0.01888,
            gammaF0 = 0.001,
            tgammaF = 36,
            Rttover = 0.017,
            MaxCond = 0.02,
            LAIgcx = 3.33,
            BLcond = 0.2,
            wSx1000 = 500,
            thinPower = 1.5,
            mF = 0.5,
            mR = 0.3,
            mS = 0.2,
            SLA0 = 5,
            SLA1 = 3,
            tSLA = 3,
            k = 0.5,
            fullCanAge = 15,
            MaxIntcptn = 0.15,
            LAImaxIntcptn = 5,
            alpha = 0.06,
            Y = 0.47,
            poolFractn = 0,
            e20 = 2.2,
            rhoAir = 1.2,
            lambda = 2460000,
            VPDconv = 0.000622,
            fracBB0 = 0.15,
            fracBB1 = 0.15,
            tBB = 10,
            rhoMin = 0.55,
            rhoMax = 0.55,
            tRho = 5,
            Qa = -90,
            Qb = 0.8,
            gDM_mol = 24,
            molPAR_MJ = 2.3,
            CoeffCond = 0.05,
            fCalpha700 = 1.433,
            fCg700 = 0.451,
            fCalphax = 2.33333333333333,
            fCg0 = 1.75,
            MinCond = 0.015,
            klmax = 0.01,
            krmax = 0.00423943,
            komax = 0.00045886,
            hc = 0.2,
            qir = 334.290515,
            qil = 49.0841127,
            qh = 23.6348669,
            qbc = 2.21427684,
            el = 0.24636719,
            er = 0.56122150,
            Nf = 0.00684,
            Navm = 0.01,
            Navx = 10,
            leaf.grow = 0,
            leaf.fall = 0,
            Wl.s = 0.01,
            Wsbr.s = 0.1,
            Wr.s = 0.01,
            pWl.sprouts = 0.5,
            pWsbr.sprouts = 0.9,
            cod.pred = "3PG",
            cod.clim = "Month" 
)
#######################################################


#######################################################
## Parameters Xenakis 2020 (unpublished) ##
#######################################################
sitka<-list(weather=clm.df.full,
            ## ~~ Initial pools ~~ ##
            Wl = 0.01,
            WlDormant = 0,
            Wr = 0.01,
            Wsbr = 0.1,
            Wlitt = 0,
            YrC = 42.95,
            YlC = 85.90,
            OC = 300.66,
            YrN = 1.43,
            YlN = 2.86,
            ON = 10.01,
            Nav = 8.76,
            ## ~~ Site ~~ ##
            N = 2000,
            rotation = 1,
            cycle = 1,
            rm.sprouts = F,
            nyears = 35,
            initial.month = 1,
            latitude = 57.06,
            soilclass = -1,
            ASW = 165,
            MaxASW = 400,
            MinASW = 0,
            CO2 = 400,
            ## ~~ Parameters ~~ ##
            pFS2 = 1.4,
            pFS20 = 0.8,
            aS = 0.138,
            nS = 2.3,
            pRx = 0.45,
            pRn = 0.3,
            Tmin = -5,
            Topt = 15,
            Tmax = 35,
            kF = 1,
            SWconst0 = 0.55,
            SWpower0 = 6,
            m0 = 0,
            MaxAge = 400,
            nAge = 4,
            rAge = 0.95,
            gammaFx = 0.01888,
            gammaF0 = 0.001,
            tgammaF = 36,
            Rttover = 0.017,
            MaxCond = 0.02,
            LAIgcx = 3.33,
            BLcond = 0.2,
            wSx1000 = 500,
            thinPower = 1.5,
            mF = 0.5,
            mR = 0.3,
            mS = 0.2,
            SLA0 = 5,
            SLA1 = 3,
            tSLA = 3,
            k = 0.5,
            fullCanAge = 15,
            MaxIntcptn = 0.15,
            LAImaxIntcptn = 5,
            alpha = 0.06,
            Y = 0.47,
            poolFractn = 0,
            e20 = 2.2,
            rhoAir = 1.2,
            lambda = 2460000,
            VPDconv = 0.000622,
            fracBB0 = 0.15,
            fracBB1 = 0.15,
            tBB = 10,
            rhoMin = 0.55,
            rhoMax = 0.55,
            tRho = 5,
            Qa = -90,
            Qb = 0.8,
            gDM_mol = 24,
            molPAR_MJ = 2.3,
            CoeffCond = 0.05,
            fCalpha700 = 1.433,
            fCg700 = 0.451,
            fCalphax = 2.33333333333333,
            fCg0 = 1.75,
            MinCond = 0.015,
            klmax = 0.01,
            krmax = 0.00423943,
            komax = 0.00045886,
            hc = 0.2,
            qir = 334.290515,
            qil = 49.0841127,
            qh = 23.6348669,
            qbc = 2.21427684,
            el = 0.24636719,
            er = 0.56122150,
            Nf = 0.00684,
            Navm = 0.01,
            Navx = 10,
            leaf.grow = 0,
            leaf.fall = 0,
            Wl.s = 0.01,
            Wsbr.s = 0.1,
            Wr.s = 0.01,
            pWl.sprouts = 0.5,
            pWsbr.sprouts = 0.9,
            cod.pred = "3PG",
            cod.clim = "Month",
            sigma_2R = 1
)
#######################################################


#Key of terminology 
#NPP - net primary production
#GPP - gross primary production
#NEE - Net ecosystem exchange
#Rs - soil respiration
#Reco - ecosystem respiration
codM<-getSample(out, start = 2000, coda = TRUE)
codM<-as.data.frame(codM[[1]])
codM<-transpose(data.frame(colMedians(codM)))
names(codM)<-nm
sitka[nm]<-codM

## Run the 3PGN model using the sitka parameters
output<-do.call(fr3PGDN,sitka)

## Plot model outputs
pOut <- plotModel(output)

## Plot the timeseries of model output vs data
results<-plotResults(output)
results[[3]]
## Calculate yield class from height
output <- output%>%mutate(
  yct = ((hdom/(1-exp(-0.033329*t.proj))^1.821054)-14.856317)/1.425397,
  YC = ifelse(yct>24,24,ifelse(yct<4,4,yct))
)


## jpeg(file="results_3PGN.jpg",width=1800,height=1000,res=100)
egg::ggarrange(plots=results,labels=c("a","b","c","d","e","f","g","h","i","j","k","l","m","n"),nrow=7,ncol=2)
## dev.off()

## jpeg(file="model_output.jpg",width=1800,height=1000,res=100)
egg::ggarrange(plots=pOut,labels=c("a","b","c","d","e","f","g","h","i","j"),nrow=5,ncol=2)
## dev.off()

model<-data.frame(
  cbind(
    filter(output,Year>=startYear&Year<=endYear)$Year,
    filter(output,Year>=startYear&Year<=endYear)$Month,
    filter(output,Year>=startYear&Year<=endYear)$GPP,
    filter(output,Year>=startYear&Year<=endYear)$NPP,
    filter(output,Year>=startYear&Year<=endYear)$NEE,
    filter(output,Year>=startYear&Year<=endYear)$Reco,
    filter(output,Year>=startYear&Year<=endYear)$Rs,
    filter(output,Year>=startYear&Year<=endYear)$EvapTransp,
    filter(output,Year>=startYear&Year<=endYear)$CanCond[2:nrow(data)]
  ))
names(model) <- c("Year","Month","GPP","NPP","NEE","Reco","Rs","EvapTransp","CanCond")
model <- model%>%mutate(timestamp=as.POSIXct(paste(sprintf("%02d",Year),sprintf("%02d",Month),sprintf("%02d",1),sep="-"),tz="GMT"))


## 1:1 plots
all <- full_join(data,model,by="timestamp")

gpp <- ggpubr::ggscatter(all,
                         x = "gpp",
                         y = "GPP",
                         add = "reg.line")+
  theme_bw()+
  geom_abline(aes(intercept=0,slope=1),lty=2)+
  ggpubr::stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"))
  )+
  scale_y_continuous(limits=c(-2,10),breaks=seq(-2,10,2),
                     name="Modelled GPP")+
  scale_x_continuous(limits=c(-2,10),breaks=seq(-2,10,2),
                     name="Observed GPP")+
  theme_bw()+
  geom_abline(aes(intercept=0,slope=1),lty=2)+
  theme(panel.spacing.x = unit(-2, "lines"),
        panel.grid=element_blank(),
        plot.margin=unit(c(0.5,0,0,0),"cm"),
        axis.text.x=element_text(size=16,colour="black"),
        axis.text.y=element_text(size=16,colour="black"),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16,colour="black"),
        strip.text=element_blank(),#element_text(size=16),
        strip.background=element_blank(),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="right")



gpp2 <- ggpubr::ggscatter(filter(all,year<2018),
                          x = "gpp",
                          y = "GPP",
                          add = "reg.line")+
  theme_bw()+
  geom_abline(aes(intercept=0,slope=1),lty=2)+
  ggpubr::stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"))
  )+
  scale_y_continuous(limits=c(-2,10),breaks=seq(-2,10,2),
                     name="Modelled GPP")+
  scale_x_continuous(limits=c(-2,10),breaks=seq(-2,10,2),
                     name="Observed GPP")+
  theme_bw()+
  annotate("text",x=40,y=80,label="Without 2018 data")+
  geom_abline(aes(intercept=0,slope=1),lty=2)+
  theme(panel.spacing.x = unit(-2, "lines"),
        panel.grid=element_blank(),
        plot.margin=unit(c(0.5,0,0,0),"cm"),
        axis.text.x=element_text(size=16,colour="black"),
        axis.text.y=element_text(size=16,colour="black"),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16,colour="black"),
        strip.text=element_blank(),#element_text(size=16),
        strip.background=element_blank(),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="right")


npp <- ggpubr::ggscatter(all,
                         x = "npp",
                         y = "NPP",
                         add = "reg.line")+
  theme_bw()+
  geom_abline(aes(intercept=0,slope=1),lty=2)+
  ggpubr::stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"))
  )+
  scale_y_continuous(limits=c(-2,10),breaks=seq(-2,10,2),
                     name="Modelled NPP")+
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,2),
                     name="Observed NPP")+
  theme_bw()+
  geom_abline(aes(intercept=0,slope=1),lty=2)+
  theme(panel.spacing.x = unit(-2, "lines"),
        panel.grid=element_blank(),
        plot.margin=unit(c(0.5,0,0,0),"cm"),
        axis.text.x=element_text(size=16,colour="black"),
        axis.text.y=element_text(size=16,colour="black"),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16,colour="black"),
        strip.text=element_blank(),#element_text(size=16),
        strip.background=element_blank(),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="right")


nee <- ggpubr::ggscatter(all,
                         x = "nee",
                         y = "NEE",
                         add = "reg.line")+
  theme_bw()+
  geom_abline(aes(intercept=0,slope=1),lty=2)+
  ggpubr::stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"))
  )+
  scale_y_continuous(limits=c(-6,2),breaks=seq(-6,2,2),
                     name="Modelled NEE")+
  scale_x_continuous(limits=c(-6,2),breaks=seq(-6,2,2),
                     name="Observed NEE")+
  theme_bw()+
  geom_abline(aes(intercept=0,slope=1),lty=2)+
  theme(panel.spacing.x = unit(-2, "lines"),
        panel.grid=element_blank(),
        plot.margin=unit(c(0.5,0,0,0),"cm"),
        axis.text.x=element_text(size=16,colour="black"),
        axis.text.y=element_text(size=16,colour="black"),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16,colour="black"),
        strip.text=element_blank(),#element_text(size=16),
        strip.background=element_blank(),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="right")

reco <- ggpubr::ggscatter(all,
                          x = "reco",
                          y = "Reco",
                          add = "reg.line")+
  theme_bw()+
  geom_abline(aes(intercept=0,slope=1),lty=2)+
  ggpubr::stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"))
  )+
  scale_y_continuous(limits=c(0,6),breaks=seq(0,6,2),
                     name="Modelled Reco")+
  scale_x_continuous(limits=c(0,6),breaks=seq(0,6,2),
                     name="Observed Reco")+
  theme_bw()+
  geom_abline(aes(intercept=0,slope=1),lty=2)+
  theme(panel.spacing.x = unit(-2, "lines"),
        panel.grid=element_blank(),
        plot.margin=unit(c(0.5,0,0,0),"cm"),
        axis.text.x=element_text(size=16,colour="black"),
        axis.text.y=element_text(size=16,colour="black"),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16,colour="black"),
        strip.text=element_blank(),#element_text(size=16),
        strip.background=element_blank(),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="right")


rs <- ggpubr::ggscatter(all,
                        x = "rs",
                        y = "Rs",
                        add = "reg.line")+
  theme_bw()+
  geom_abline(aes(intercept=0,slope=1),lty=2)+
  ggpubr::stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"))
  )+
  scale_y_continuous(limits=c(0,3),breaks=seq(0,3,1),
                     name="Modelled Rs")+
  scale_x_continuous(limits=c(0,3),breaks=seq(0,3,1),
                     name="Observed Rs")+
  theme_bw()+
  geom_abline(aes(intercept=0,slope=1),lty=2)+
  theme(panel.spacing.x = unit(-2, "lines"),
        panel.grid=element_blank(),
        plot.margin=unit(c(0.5,0,0,0),"cm"),
        axis.text.x=element_text(size=16,colour="black"),
        axis.text.y=element_text(size=16,colour="black"),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16,colour="black"),
        strip.text=element_blank(),#element_text(size=16),
        strip.background=element_blank(),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="right")


et <- ggpubr::ggscatter(all,
                        x = "et",
                        y = "EvapTransp",
                        add = "reg.line")+
  theme_bw()+
  geom_abline(aes(intercept=0,slope=1),lty=2)+
  ggpubr::stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"))
  )+
  scale_y_continuous(limits=c(30,100),breaks=seq(30,100,10),
                     name="Modelled ET")+
  scale_x_continuous(limits=c(30,100),breaks=seq(30,100,10),
                     name="Observed ET")+
  theme_bw()+
  annotate("text",x=40,y=80,label="With 2018 data")+
  geom_abline(aes(intercept=0,slope=1),lty=2)+
  theme(panel.spacing.x = unit(-2, "lines"),
        panel.grid=element_blank(),
        plot.margin=unit(c(0.5,0,0,0),"cm"),
        axis.text.x=element_text(size=16,colour="black"),
        axis.text.y=element_text(size=16,colour="black"),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16,colour="black"),
        strip.text=element_blank(),#element_text(size=16),
        strip.background=element_blank(),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="right")

et2 <- ggpubr::ggscatter(filter(all,year<2018),
                         x = "et",
                         y = "EvapTransp",
                         add = "reg.line")+
  theme_bw()+
  geom_abline(aes(intercept=0,slope=1),lty=2)+
  ggpubr::stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"))
  )+
  scale_y_continuous(limits=c(30,100),breaks=seq(30,100,10),
                     name="Modelled ET")+
  scale_x_continuous(limits=c(30,100),breaks=seq(30,100,10),
                     name="Observed ET")+
  theme_bw()+
  annotate("text",x=40,y=80,label="Without 2018 data")+
  geom_abline(aes(intercept=0,slope=1),lty=2)+
  theme(panel.spacing.x = unit(-2, "lines"),
        panel.grid=element_blank(),
        plot.margin=unit(c(0.5,0,0,0),"cm"),
        axis.text.x=element_text(size=16,colour="black"),
        axis.text.y=element_text(size=16,colour="black"),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16,colour="black"),
        strip.text=element_blank(),#element_text(size=16),
        strip.background=element_blank(),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="right")





## jpeg(file="modelled_observed.jpg",width=1200,height=1200,res=100)
egg::ggarrange(gpp,npp,nee,reco,rs,et,et2,
               labels=c("(a)","(b)","(c)","(d)","(e)","(f)","(g)"),nrow=4,ncol=2)
## dev.off()









## m <- as.matrix(out$chain)[1000001:1000100,1:(ncol(as.matrix(out$chain))-3)]

## multiRunModel <- function(par,names,p,var){
##     par[names] <- p
##     output <- do.call(fr3PGDN,par)
##     return(output[[var]])
## }

## gpp <- apply(m,1,FUN=multiRunModel,par=sitka,names=nm,var="GPP")

## matplot(gpp,type='l')

