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
pineLL <- function(p){
  pine[.GlobalEnv$nm]<-p
  
  
  NlogLik <- tryCatch(
    {
      output<-do.call(SoNWaL,pine)
      modelled <-suppressWarnings(suppressMessages(sampleOutputPine(output,.GlobalEnv$startYear,.GlobalEnv$endYear)))
      NlogLik  <-   ifelse(any(is.na(modelled)==T),-Inf,(flogL(data=.GlobalEnv$observedPine,sims=modelled,data_s=.GlobalEnv$devPine)))
      NlogLik<-ifelse(pine$fieldCap<pine$wiltPoint,-Inf,NlogLik)
      NlogLik<-ifelse(pine$satPoint<pine$fieldCap,-Inf,NlogLik)
      NlogLik<-ifelse(pine$fieldCap<pine$wiltPoint,-Inf,NlogLik)
      
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
      output<-do.call(SoNWaL,pine)
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
      output<-do.call(SoNWaL,pine)
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
      output<-do.call(SoNWaL,pine)
      modelled <-suppressWarnings(suppressMessages(sampleOutputPineWeek_noHyd(output,.GlobalEnv$startYear,.GlobalEnv$endYear)))
      ifelse(any(is.na(modelled)==TRUE),-Inf,sum(dnorm(x=.GlobalEnv$observedPine,sd =.GlobalEnv$devPine, mean=modelled,log=T),na.rm = T))
      
    },
    error=function(cond) {
      return(-Inf)
    })
  
  return(NlogLik)
}

