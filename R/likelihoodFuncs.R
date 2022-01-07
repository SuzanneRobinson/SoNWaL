




## Extract simulated data for use in likelihood function for shorter time-steps
#daily (and therefore weekly) data is in grams per meter square per day, so need to convert 3pg output from tons per hectare
sampleOutputTS<-function(df,sY,eY){
  df<- filter(df,Year>=sY&Year<=eY)
  m<-c(aggregate(df$GPP~ df$Month+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$NPP~ df$Month+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$NEE~ df$Month+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$Reco~ df$Month+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$Rs~ df$Month+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$Etransp~ df$Month+df$Year,FUN=sum)[,3],
       filter(df,Year==2015&Month==8)$LAI[1],
       filter(df,Year==2018&Month==8)$LAI[1],
       filter(df,Year==2018&Month==8)$N[1],
       filter(df,Year==2018&Month==8)$dg[1],
     #  filter(df,Year==2015&Month==7)$Wr[1],
     #  filter(df,Year==2015&Month==7)$difRoots[1],
       filter(df,Year==2015&Month==7)$totC[1],
       filter(df,Year==2015&Month==7)$totN[1],
       aggregate(df$volSWC_rz~ df$Month+df$Year,FUN=mean)[,3]
  )
  m
  return(m)
}

## Extract simulated data for use in likelihood function for shorter time-steps
#monthly data is in tons per hectare per month which matches 3pg output so no need to convert
sampleOutputMonth<-function(df,sY,eY){
  #convert to average grams per m2 per day depending on timestep of model
  modif<- if(nrow(df)<600) 1.6 else 7.142857
  nDays<- if(nrow(df)<600) 30 else 7
  
     df<-df%>%
      filter(Year>=sY&Year<=eY)%>%
      mutate(GPP=GPP*modif)%>%
      mutate(NEE=NEE*modif)%>%
      mutate(Rs=Rs*modif)%>%
      mutate(Reco=Reco*modif)
       
  
  m<-c(#aggregate(df$Reco/df$Rs~ df$Year,FUN=mean)[,2],
       aggregate(df$GPP~ df$Month+df$Year,FUN=mean)[,3],
       #aggregate(df$NPP~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$NEE~ df$Month+df$Year,FUN=mean)[,3],
      # aggregate(df$Reco~ df$Month+df$Year,FUN=mean)[,3],
       #aggregate(df$Rs~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$EvapTransp/nDays~ df$Month+df$Year,FUN=mean)[,3],
       filter(df,Year==2015&Month==8)$LAI[1],
       filter(df,Year==2018&Month==8)$LAI[1],
       filter(df,Year==2018&Month==8)$N[1],
       filter(df,Year==2018&Month==8)$dg[1],
       #  filter(df,Year==2015&Month==7)$Wr[1],
       #  filter(df,Year==2015&Month==7)$difRoots[1],
       filter(df,Year==2015&Month==7)$totC[1],
       filter(df,Year==2015&Month==7)$totN[1],
       aggregate(df$volSWC_rz~ df$Month+df$Year,FUN=mean)[,3],
      aggregate(df$Rs~ df$Month+df$Year,FUN=mean)[,3]
      
  )
  m
  return(m)
}



## Extract simulated data for use in likelihood function for shorter time-steps
#IF USING WEEKLY/daily, YOU USE DAILY DATA WHICH IS IN gCm-2d-1, so need to convert prafor output from tons per hectare per week (timestep)
sampleOutputWeekly<-function(df,sY,eY){
  modif<- 7.142857
  
  df<- filter(df,Year>=sY&Year<=eY)
  df$week<-c(1:53)
  
  df<-df%>%
    mutate(GPP=GPP*modif)%>%
    mutate(NPP=NPP*modif)%>%
    mutate(NEE=NEE*modif)%>%
    mutate(Rs=Rs*modif)%>%
    mutate(Reco=Reco*modif)
  
  
  m<-c(aggregate(df$Reco/df$Rs~ df$Year,FUN=mean)[,2],
       aggregate(df$GPP~ df$week+df$Year,FUN=mean)[,3],
       #aggregate(df$NPP~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$NEE~ df$week+df$Year,FUN=mean)[,3],
       # aggregate(df$Reco~ df$Month+df$Year,FUN=mean)[,3],
       #aggregate(df$Rs~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$EvapTransp~ df$week+df$Year,FUN=mean)[,3]/7,
       filter(df,Year==2015&Month==8)$LAI[1],
       filter(df,Year==2018&Month==8)$LAI[1],
       filter(df,Year==2018&Month==8)$N[1],
       filter(df,Year==2018&Month==8)$dg[1],
       #  filter(df,Year==2015&Month==7)$Wr[1],
       #  filter(df,Year==2015&Month==7)$difRoots[1],
       filter(df,Year==2015&Month==7)$totC[1],
       filter(df,Year==2015&Month==7)$totN[1],
       aggregate(df$volSWC_rz~ df$week+df$Year,FUN=mean)[,3]
  )
  m
  return(m)
}






## Extract simulated data for use in likelihood function for shorter time-steps
sampleOutputTS_noHyd<-function(df,sY,eY){
  modif<- if(nrow(df)<600) 1.6
  modif<- if(nrow(df)>600) 7.142857
  df<-df%>%
    mutate(GPP=GPP*modif)%>%
    mutate(NPP=NPP*modif)%>%
    mutate(NEE=NEE*modif)%>%
    mutate(Rs=Rs*modif)%>%
    mutate(Reco=Reco*modif)
  
  df<- filter(df,Year>=sY&Year<=eY)
  
  m<-c(aggregate(df$GPP~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$NPP~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$NEE~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$Reco~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$Rs~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$Etransp~ df$Month+df$Year,FUN=sum)[,3],
       filter(df,Year==2015&Month==8)$LAI[1],
       filter(df,Year==2018&Month==8)$LAI[1],
       filter(df,Year==2018&Month==8)$N[1],
       filter(df,Year==2018&Month==8)$dg[1],
       #  filter(df,Year==2015&Month==7)$Wr[1],
       #  filter(df,Year==2015&Month==7)$difRoots[1],
       filter(df,Year==2015&Month==7)$totC[1],
       filter(df,Year==2015&Month==7)$totN[1]  )
  m
  return(m)
}


## Extract simulated data for use in likelihood function for shorter time-steps
sampleOutputTS_noHyd_weekly<-function(df,sY,eY){
  df<- filter(df,Year>=sY&Year<=eY)
  df$week<-c(1:53)
  
  m<-c(aggregate(df$GPP~ df$week+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$NPP~ df$week+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$NEE~ df$week+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$Reco~ df$week+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$Rs~ df$week+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$Etransp~ df$week+df$Year,FUN=sum)[,3],
       filter(df,Year==2015&Month==8)$LAI[1],
       filter(df,Year==2018&Month==8)$LAI[1],
       filter(df,Year==2018&Month==8)$N[1],
       filter(df,Year==2018&Month==8)$dg[1],
       #  filter(df,Year==2015&Month==7)$Wr[1],
       #  filter(df,Year==2015&Month==7)$difRoots[1],
       filter(df,Year==2015&Month==7)$totC[1],
       filter(df,Year==2015&Month==7)$totN[1]  )
  m
  return(m)
}



## Extract simulated data for use in likelihood function
#2007 is removed as a year due to some odd observed data values, basically there's no GPP or NEE etc.
sampleOutputPineWeek<-function(df,sY,eY){
  df<-df[-1,]
df$week<-c(1:52)
  
  m<-c(
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,week)%>%
           dplyr::summarise(mean=mean(GPP,na.rm=TRUE))%>%
           select(mean))*100/7,
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,week)%>%
           dplyr::summarise(mean=mean(NEE,na.rm=TRUE))%>%
           select(mean))*100/7,
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,week)%>%
           dplyr::summarise(mean=mean(Reco,na.rm=TRUE))%>%
           select(mean))*100/7,
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(LAI,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(dg,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year==1996)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(totC,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year==1996)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(totN,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(N,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1997&Year<=2014&Year!=2007)%>%
           group_by(Year,week)%>%
           dplyr::summarise(mean=mean(volSWC_rz,na.rm=TRUE))%>%
           select(mean))
    
    
  )
  return(m)
}

sampleOutputPine<-function(df,sY,eY){
  #convert to match hytialla data format

  modif<- if(nrow(df)<600) 1.6
  modif<- if(nrow(df)>600) 7.142857
  df<-df%>%
    mutate(GPP=GPP*modif)%>%
    mutate(NPP=NPP*modif)%>%
    mutate(NEE=NEE*modif)%>%
    mutate(Rs=Rs*modif)%>%
    mutate(Reco=Reco*modif)
    
  m<-c(
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,Month)%>%
           dplyr::summarise(mean=mean(GPP,na.rm=TRUE))%>%
           select(mean)),
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,Month)%>%
           dplyr::summarise(mean=mean(NEE,na.rm=TRUE))%>%
           select(mean)),
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,Month)%>%
           dplyr::summarise(mean=mean(Reco,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(LAI,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(dg,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year==1996)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(totC,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year==1996)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(totN,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(N,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1997&Year<2014)%>%
           group_by(Year,Month)%>%
           dplyr::summarise(mean=mean(volSWC_rz,na.rm=TRUE))%>%
           select(mean))
    
    
  )
  return(m)
}




## Extract simulated data for use in likelihood function
#2007 is removed as a year due to some odd observed data values, basically there's no GPP or NEE etc.
sampleOutputPineWeek_noHyd<-function(df,sY,eY){
  df<-df[-1,]
  df$week<-c(1:52)
  
  m<-c(
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,week)%>%
           dplyr::summarise(mean=mean(GPP,na.rm=TRUE))%>%
           select(mean))*100/7,
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,week)%>%
           dplyr::summarise(mean=mean(NEE,na.rm=TRUE))%>%
           select(mean))*100/7,
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,week)%>%
           dplyr::summarise(mean=mean(Reco,na.rm=TRUE))%>%
           select(mean))*100/7,
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(LAI,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(dg,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year==1996)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(totC,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year==1996)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(totN,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(N,na.rm=TRUE))%>%
           select(mean))
    
    
  )
  return(m)
}

sampleOutputPine_noHyd<-function(df,sY,eY){
  #convert to match hytialla data format
  modif<- if(nrow(df)<600) 1.6
  modif<- if(nrow(df)>600) 7.142857
  df<-df%>%
    mutate(GPP=GPP*modif)%>%
    mutate(NPP=NPP*modif)%>%
    mutate(NEE=NEE*modif)%>%
    mutate(Rs=Rs*modif)%>%
    mutate(Reco=Reco*modif)
  
  m<-c(
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,Month)%>%
           dplyr::summarise(mean=mean(GPP,na.rm=TRUE))%>%
           select(mean)),
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,Month)%>%
           dplyr::summarise(mean=mean(NEE,na.rm=TRUE))%>%
           select(mean)),
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,Month)%>%
           dplyr::summarise(mean=mean(Reco,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(LAI,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(dg,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year==1996)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(totC,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year==1996)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(totN,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(N,na.rm=TRUE))%>%
           select(mean))
    
  )
  return(m)
}


## Likelihood function
NLL<- function(p){
  sitka[.GlobalEnv$nm]<-p
  
  #sitka$Q10X<-0
  
   NlogLik <- tryCatch(
    {
    #  runMod<-function(parms){
      
      
      output<-   do.call(fr3PGDN,sitka)
      #}
      
    #  siteLst<-list(sitka,sitka,sitka,sitka,sitka,sitka,sitka)
      
    #  plan(multisession, workers = 7)
      
      #system.time(future_map(siteLst,runMod))
      #use sampleOutputTS if using smaller time-steps
      modelled <-sampleOutputMonth(output,.GlobalEnv$startYear,.GlobalEnv$endYear)
      
      
     #  NlogLik  <-   ifelse(any(is.na(modelled)==T),-Inf,sum(dnorm(modelled,mean=.GlobalEnv$observed,sd=.GlobalEnv$dev,log=T),na.rm = T))
     NlogLik  <-   ifelse(any(is.na(modelled)==T),-Inf,flogL(data=.GlobalEnv$observed,sims=modelled,data_s=.GlobalEnv$dev))
      
      NlogLik<-ifelse(max(output$LAI)>8,-Inf,NlogLik)
      NlogLik<-ifelse(min(output$totN)<1,-Inf,NlogLik)
      NlogLik<-ifelse(sitka$fieldCap<sitka$wiltPoint,-Inf,NlogLik)
      NlogLik<-ifelse(sitka$satPoint<sitka$fieldCap,-Inf,NlogLik)
    },
    error=function(cond) {
      return(-Inf)
    })
  return(NlogLik)
}


flogL <- function(sims,data,data_s)
{ 
  Ri         <- (sims - data) / data_s
  i0         <- which( abs(Ri)<1.e-08 )
  
  logLi      <- log(1-exp(-0.5*Ri^2)) - log(Ri^2) - 0.5*log(2*pi) - log(data_s)
  logLi[i0]  <- -0.5*log(2*pi) - log(2*data_s[i0])
  
  sum(logLi)
}


## Likelihood function
NLL_weekly<- function(p){
  sitka[.GlobalEnv$nm]<-p
  
  
  NlogLik <- tryCatch(
    {
      output<-do.call(fr3PGDN,sitka)
      #use sampleOutputTS if using smaller time-steps
      modelled <-sampleOutputWeekly(output,.GlobalEnv$startYear,.GlobalEnv$endYear)
      
      NlogLik  <-   ifelse(any(is.na(modelled)==T),-Inf,sum(dnorm(modelled,.GlobalEnv$observed,GlobalEnv$dev,log=T)))
      
NlogLik<-ifelse(max(output$LAI)>8,-Inf,NlogLik)
#        NlogLik<-ifelse(mean(tail(output$LAI,500))<1,-Inf,NlogLik)
      NlogLik<-ifelse(min(output$totN)<1,-Inf,NlogLik)
      
      
    },
    error=function(cond) {
      return(-Inf)
    })
  return(NlogLik)
}


## Likelihood function
NLL_noHYD<- function(p){
  sitka[.GlobalEnv$nm]<-p
  
  NlogLik <- tryCatch(
    {
      output<-do.call(fr3PGDN,sitka)
      #use sampleOutputTS if using smaller time-steps
      modelled <-sampleOutputTS_noHyd(output,.GlobalEnv$startYear,.GlobalEnv$endYear)
      
      NlogLik  <-   ifelse(any(is.na(modelled)==T),-Inf,sum(dnorm(.GlobalEnv$observed,mean=modelled,sd=.GlobalEnv$dev,log=T),na.rm = T))
      
      
    },
    error=function(cond) {
      return(-Inf)
    })
  return(NlogLik)
}


## Likelihood function
NLL_noHYDWeekly<- function(p){
  sitka[.GlobalEnv$nm]<-p
  
  NlogLik <- tryCatch(
    {
      output<-do.call(fr3PGDN,sitka)
      #use sampleOutputTS if using smaller time-steps
      modelled <-sampleOutputTS_noHyd_weekly(output,.GlobalEnv$startYear,.GlobalEnv$endYear)
      
      NlogLik  <-   ifelse(any(is.na(modelled)==T),-Inf,sum(dnorm(.GlobalEnv$observed,mean=modelled,sd=.GlobalEnv$dev,log=T),na.rm = T))
      
      
    },
    error=function(cond) {
      return(-Inf)
    })
  return(NlogLik)
}


## Likelihood function
pineLL <- function(p){
  pine[.GlobalEnv$nm]<-p
  
  
  NlogLik <- tryCatch(
    {
      output<-do.call(fr3PGDN,pine)
      modelled <-suppressWarnings(suppressMessages(sampleOutputPine(output,.GlobalEnv$startYear,.GlobalEnv$endYear)))
      NlogLik  <-   ifelse(any(is.na(modelled)==T),-Inf,(flogL(data=.GlobalEnv$observedPine,sims=modelled,data_s=.GlobalEnv$devPine)))
      
    },
    error=function(cond) {
      return(-Inf)
    })
  
  return(NlogLik)
}


## Likelihood function
pineLL_weekly <- function(p){
  pine[.GlobalEnv$nm]<-p
  
  
  NlogLik <- tryCatch(
    {
      output<-do.call(fr3PGDN,pine)
      modelled <-suppressWarnings(suppressMessages(sampleOutputPineWeek(output,.GlobalEnv$startYear,.GlobalEnv$endYear)))
      ifelse(any(is.na(modelled)==TRUE),-Inf,sum(dnorm(x=.GlobalEnv$observedPine,sd =.GlobalEnv$devPine, mean=modelled,log=T),na.rm = T))
      
    },
    error=function(cond) {
      return(-Inf)
    })
  
  return(NlogLik)
}




## Likelihood function
pineLL_noHyd <- function(p){
  pine[.GlobalEnv$nm]<-p
  
  
  NlogLik <- tryCatch(
    {
      output<-do.call(fr3PGDN,pine)
      modelled <-suppressWarnings(suppressMessages(sampleOutputPine_noHyd(output,.GlobalEnv$startYear,.GlobalEnv$endYear)))
      ifelse(any(is.na(modelled)==TRUE),-Inf,sum(dnorm(x=.GlobalEnv$observedPine,sd =.GlobalEnv$devPine, mean=modelled,log=T),na.rm = T))
      
    },
    error=function(cond) {
      return(-Inf)
    })
  
  return(NlogLik)
}


## Likelihood function
pineLL_weekly_noHyd <- function(p){
  pine[.GlobalEnv$nm]<-p
  
  
  NlogLik <- tryCatch(
    {
      output<-do.call(fr3PGDN,pine)
      modelled <-suppressWarnings(suppressMessages(sampleOutputPineWeek_noHyd(output,.GlobalEnv$startYear,.GlobalEnv$endYear)))
      ifelse(any(is.na(modelled)==TRUE),-Inf,sum(dnorm(x=.GlobalEnv$observedPine,sd =.GlobalEnv$devPine, mean=modelled,log=T),na.rm = T))
      
    },
    error=function(cond) {
      return(-Inf)
    })
  
  return(NlogLik)
}




## Likelihood function
NLL_Reg<- function(p){
 px<-data.frame(t(p))
  fail<-rep(-Inf,nrow(px))
 siteLst<-(unique(sitka$weather$site))
 siteVec<-rep(siteLst,each=nrow(px))
 px$chain<-rep(1:nrow(px))
 
 px<-px[rep(seq_len(nrow(px)), length(siteLst)), ]
 
 px$site<-siteVec


      runMod<-function(g){
       res<-tryCatch(
         { 
        sitka[.GlobalEnv$nm]<-g[1:16]
        sitka$weather<-filter(sitka$weather,site==g$site)
        
        if(sitka$weather$soilDepth[1]==1){ 
          sitka$V_nr<- 0.25
          sitka$maxRootDepth<- 0.25
        }
         
        if(sitka$weather$soilDepth[1]==2){
          sitka$V_nr<- 0.5
          sitka$maxRootDepth<- 0.5
          }
        if(sitka$weather$soilDepth[1]==3){
          sitka$V_nr<- 0.8
          sitka$maxRootDepth<- 0.8
         }
        if(sitka$weather$soilDepth[1]==4){
          sitka$V_nr<- 1
          sitka$maxRootDepth<- 1
          }
        if(sitka$weather$soilDepth[1]==5){
          sitka$V_nr<- 1.5
          sitka$maxRootDepth<- 1.5
       }
        
        if(sitka$weather$soilTex[1]!=0){
          sitka$wiltPoint<-sitka$weather$wp[1]
          sitka$fieldCap<-sitka$weather$fc[1]
          sitka$satPoint<-sitka$weather$sp[1]
          sitka$K_s<-sitka$weather$soilCond[1]
        sitka$K_drain<-sitka$weather$soilCond[1]
        }
        if(sitka$weather$soilTex[1]==1) {
         sitka$E_S1<-0.05
         sitka$E_S2<-0.3
         sitka$SWpower0<-8
         sitka$SWconst0<-0.65
       }
        
        if(sitka$weather$soilTex[1]==0) {
         sitka$E_S1<-0.05
         sitka$E_S2<-0.3
         sitka$SWpower0<-8
         sitka$SWconst0<-0.65
         }
        
        if(is.na(sitka$weather$soilTex[1])==T)  {
         sitka$E_S1<-0.05
         sitka$E_S2<-0.3
         sitka$SWpower0<-8
         sitka$SWconst0<-0.65
          }
        
        if(sitka$weather$soilTex[1]==4)  {
         sitka$E_S1<-0.1
         sitka$E_S2<-0.3
         sitka$SWpower0<-5
         sitka$SWconst0<-0.5
         }
        
        
        if(sitka$weather$soilTex[1]==9) {
          sitka$E_S1<-0.3
          sitka$E_S2<-0.6
          sitka$SWpower0<-5
          sitka$SWconst0<-0.5
          }
        
        
        
        observed<-filter(sitka$weather,is.na(mean_dbh_cm)==F)%>%group_by(Year)%>%summarise(dbh=median(mean_dbh_cm),dbhSD=median(dbhSD_cm))
        observed$dbhSD<-ifelse(observed$dbhSD==0,0.0001,observed$dbhSD)
        sY=min(observed$Year)
        eY=max(observed$Year)
        output<-   do.call(fr3PGDN,sitka)
        
        modelled<-output%>%group_by(Year)%>%summarise(dg=mean(dg))%>%filter(Year>=sY&Year<=eY)
        ifelse(any(is.na(modelled)==T),-Inf,flogL(data=observed$dbh,sims=modelled$dg,data_s=observed$dbhSD))
         },
        error=function(cond) {
          return(-Inf)
        })
  return(res)
      }
      
    g <- split(px, seq(nrow(px)))
     px$ll<-do.call(rbind,mclapply(g, FUN = runMod,mc.cores=16))
      NlogLik<-as.vector(px%>%group_by(chain)%>%summarise(ll=sum(ll,na.rm=T))%>%pull(ll))
    

  return(NlogLik)
}




## Likelihood function
NLL_Reg_mp<- function(px){
  px<-as.data.frame(px)
  #names(px)<-nm
  #sitka$Q10X<-0
  #matrix rows = chains, cols = params
  fail<-rep(-Inf,nrow(px))
  
  siteLst<-(unique(sitka$weather$site))
  siteVec<-rep(siteLst,each=nrow(px))
  px$chain<-rep(1:nrow(px))
  
  px<-px[rep(seq_len(nrow(px)), length(siteLst)), ]
  
  px$site<-siteVec
  #  
  #  cl <- startMPIcluster(count=47)
  #  registerDoMPI(cl)
  #  
  #
  #  foreach(i = c(1:nrow(px))) %dopar%{
  #library(dplyr)
  #    library(fr3PGDN)
  #        pt<-px[i,]
  #        sitkaX<-sitka
  #        sitkaX[.GlobalEnv$nm]<-as.numeric(pt[.GlobalEnv$nm])
  #        sitkaX$weather<-dplyr::filter(sitkaX$weather,site==pt$site)
  #   observed<-dplyr::filter(sitkaX$weather,is.na(mean_dbh_cm)==F)%>%group_by(Year)%>%summarise(dbh=median(mean_dbh_cm),dbhSD=median(dbhSD_cm))
  #   observed$dbhSD<-ifelse(observed$dbhSD==0,0.0001,observed$dbhSD)
  #   sY=min(observed$Year)
  #   eY=max(observed$Year)
  #   output<-   do.call(fr3PGDN,sitkaX)
  #   
  #   modelled<-output%>%dplyr::group_by(Year)%>%dplyr::summarise(dg=mean(dg))%>%dplyr::filter(Year>=sY&Year<=eY)
  #   ifelse(any(is.na(modelled)==T),-Inf,flogL(data=observed$dbh,sims=modelled$dg,data_s=observed$dbhSD))
  #    
  #  }

  p <- split(px, seq(nrow(px)))
  px$ll<-do.call(rbind,mclapply(p, FUN = runMod,mc.cores = 16))
  NlogLik<-as.vector(px%>%group_by(chain)%>%summarise(ll=sum(ll,na.rm=T))%>%pull(ll))
  
  
  return(NlogLik)
}


