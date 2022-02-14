




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
#function relies on the parameter list (paramList) defined at the submission of the mcmc, this list is what is submitted to the model run function to run SoNWal
#it contains elements for each parameter to be updated from the MCMC chain proposals
#additionally it contains the paramList$weather dataframe - this contains the longitudinal climate data and site specific data for each of the 13 regional sites
#'@param p matrix of parameter proposals for each chain, matrix: ncols = parameter number, nrows = number of internal chains
#'@return vector of likelihood values equal to number of rows in input p matrix
NLL_Reg_sitka <- function(p) {
    #Sometimes number of rows in p matrix is less than number of chains (uncertain if algorithm design or bug)
  #this means that sometimes it comes in as a single row, basically a vector and needs transposing when converting to dataframe (or it does for the way i've done things)
  px <- if (is.null(nrow(p)) == F)
    data.frame(p) else
    data.frame(t(p))
  #nm is a vector of the names of parameter values being fitted  e.g.
  #nm<-c("V_nr","sigma_zR","shared_area","maxRootDepth","K_drain",
  #      "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR"...)
  names(px) <- nm
  paramList<-sitka
  
  siteSoil<-group_by(paramList$weather,site)%>%summarise(nutrients=median(SoilNutrientRegime))
  
  #paramList$weather dataframe includes climate and site specific info for all regional sites
  #There are currently two soil nutrient regimes in the data "poor" and "very poor", in SoNWaL we fit starting N and C
  #To reduce fitting parameters for each site (and thus having lots of extra params to fit), I assume that sites with matching soil nutrient regime share N and C start values
  #I then fit starting C and N for the "poor" sites and fit a single modifer "poorSoilMod" (val between 0.01-0.95) which reduces start C and N for "very poor" sites - could be better ways of doing this
  
  #create vector of -Inf vals to return if the model run function fails
  fail <- rep(-Inf, nrow(px))
  #Get list of site names from the paramList$weather dataframe
  siteLst <- (unique(paramList$weather$site))
  #Update dataframe with each chains proposed param values repeated for each site
  siteVec <- rep(siteLst, each = nrow(px))
  px$chain <-
    rep(1:nrow(px))  #add chain number for reference and aggregating later
  px <- px[rep(seq_len(nrow(px)), length(siteLst)),]
  px$site <- siteVec
  
  px<-merge(px,siteSoil,by.x="site",by.y="site")

  #split px dataframe into list for running in parallel lapply function
  splitParams <- split(px, seq(nrow(px)))
#  splitParamsX<<-splitParams
 # pxX<<-px
#  paramListX<<-paramList
  #update px dataframe with likelihood values for each site
  px$ll <-
    do.call(rbind, mcmapply(runMod,splitParams,MoreArgs = list(paramList),SIMPLIFY = F, mc.cores = 15))
  
  #sum likelihood values by chain and return vector of likelihood values, one for each chain being run
  NlogLik <-
    as.vector(px %>% group_by(chain) %>% summarise(ll = sum(ll)) %>%
                pull(ll))
  
  NlogLik[is.na(NlogLik)==T]<--Inf
  return(NlogLik)
  
}


#runMod function
#'@param newParams list element containing dataframe of proposed parameter values for single chain for single site
runMod <- function(newParams,paramListX) {
  res <- tryCatch({
    #update parameters (-45, as poorSoilMod is not a model parameter..this will depend on your nm param vector of course)
    nmX<-c("sigma_zR","shared_area",
           "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
           "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
           "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","Qa","Qb","MaxIntcptn","k","startN","startC")
    
    paramListX[nmX] <- newParams[nmX]
    #filter paramList$weather dataframe by site
    paramListX$weather <- filter(paramListX$weather, site == newParams$site)
    #modify start C and start N by soil modifier if soil regime is "very poor"
    paramListX$startC <- paramListX$startC * ifelse(newParams$nutrients=="Poor",1, newParams$poorSoilMod)
    paramListX$startN <- paramListX$startN * ifelse(newParams$nutrients=="Poor",1, newParams$poorSoilMod)
    
    #update fixed parameter values based on site specific values contained in paramList$weather dataframe
    if (paramListX$weather$soilDepth[1] == 1) {
      paramListX$V_nr <- 0.25
      paramListX$maxRootDepth <- 0.25
    }
    
    if (paramListX$weather$soilDepth[1] == 2) {
      paramListX$V_nr <- 0.5
      paramListX$maxRootDepth <- 0.5
    }
    if (paramListX$weather$soilDepth[1] == 3) {
      paramListX$V_nr <- 0.8
      paramListX$maxRootDepth <- 0.8
    }
    if (paramListX$weather$soilDepth[1] == 4) {
      paramListX$V_nr <- 1
      paramListX$maxRootDepth <- 1
    }
    if (paramListX$weather$soilDepth[1] == 5) {
      paramListX$V_nr <- 1.5
      paramListX$maxRootDepth <- 1.5
    }
    
    if (paramListX$weather$soilTex[1] != 0) {
      paramListX$wiltPoint <- paramListX$weather$wp[1]
      paramListX$fieldCap <- paramListX$weather$fc[1]
      paramListX$satPoint <- paramListX$weather$sp[1]
      paramListX$K_s <- paramListX$weather$soilCond[1]
      paramListX$K_drain <- paramListX$weather$soilCond[1]
    }
    if (paramListX$weather$soilTex[1] == 1) {
      paramListX$E_S1 <- 0.05
      paramListX$E_S2 <- 0.3
      paramListX$SWpower0 <- 8
      paramListX$SWconst0 <- 0.65
    }
    if (paramListX$weather$soilTex[1] == 0) {
      paramListX$E_S1 <- 0.05
      paramListX$E_S2 <- 0.3
      paramListX$SWpower0 <- 8
      paramListX$SWconst0 <- 0.65
    }
    if (is.na(paramListX$weather$soilTex[1]) == T)  {
      paramListX$E_S1 <- 0.05
      paramListX$E_S2 <- 0.3
      paramListX$SWpower0 <- 8
      paramListX$SWconst0 <- 0.65
    }
    
    if (paramListX$weather$soilTex[1] == 4)  {
      paramListX$E_S1 <- 0.1
      paramListX$E_S2 <- 0.3
      paramListX$SWpower0 <- 5
      paramListX$SWconst0 <- 0.5
    }
    if (paramListX$weather$soilTex[1] == 9) {
      paramListX$E_S1 <- 0.3
      paramListX$E_S2 <- 0.6
      paramListX$SWpower0 <- 5
      paramListX$SWconst0 <- 0.5
    }
    
    #get observed values (also contained in paramListX$weather dataframe) for site being fitted
    observed <-
      filter(paramListX$weather, is.na(mean_dbh_cm) == F) %>% group_by(Year) %>% summarise(dbh =
                                                                                             median(mean_dbh_cm),
                                                                                           dbhSD = median(dbhSD_cm))
    observed$dbhSD <-
      ifelse(observed$dbhSD == 0, 0.0001, observed$dbhSD)
    sY = min(observed$Year)
    eY = max(observed$Year)
    
    #filter climate data by planting year so model runs from planting year
    paramListX$weather <-
      filter(paramListX$weather, Year >= paramListX$weather$plantingYear[1])
    
    #run model
    output <- do.call(fr3PGDN, paramListX)
    
    #filter simulated data to match observed data format
    modelled <-
      output %>% group_by(Year) %>% summarise(dg = mean(dg)) %>% filter(Year >=
                                                                          sY & Year <= eY)
    
    #run likelihood function of observed vs simulated and get liklelihood value
    ifelse(
      any(is.na(modelled) == T),
      -Inf,
      flogL(
        data = observed$dbh,
        sims = modelled$dg,
        data_s = observed$dbhSD
      )
    )
  },
  error = function(cond) {
    #return -Inf if something goes wrong with param proposals
    return(-Inf)
  })
  return(res)
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






##############likelihood for mixed prior regional calibrations####################



## Likelihood function
#function relies on the parameter list (paramList) defined at the submission of the mcmc, this list is what is submitted to the model run function to run SoNWal
#it contains elements for each parameter to be updated from the MCMC chain proposals
#additionally it contains the paramList$weather dataframe - this contains the longitudinal climate data and site specific data for each of the 13 regional sites
#'@param p matrix of parameter proposals for each chain, matrix: ncols = parameter number, nrows = number of internal chains
#'@return vector of likelihood values equal to number of rows in input p matrix
NLL_Reg_sitka_mixedP <- function(p) {
  #Sometimes number of rows in p matrix is less than number of chains (uncertain if algorithm design or bug)
  #this means that sometimes it comes in as a single row, basically a vector and needs transposing when converting to dataframe (or it does for the way i've done things)
  px <- if (is.null(nrow(p)) == F)
    data.frame(p) else
      data.frame(t(p))

  nm_all<-c(paste0("wiltPoint_Si",unique(paramList$weather$site)),
             paste0("fieldCap_Si",unique(paramList$weather$site)),
             paste0("satPoint_Si",unique(paramList$weather$site)),
             paste0("K_s_Si",unique(paramList$weather$site)),
             paste0("V_nr_Si",unique(paramList$weather$site)),
             paste0("E_S1_Si",unique(paramList$weather$site)),
             paste0("E_S2_Si",unique(paramList$weather$site)),
             paste0("shared_area_Si",unique(paramList$weather$site)),
             paste0("maxRootDepth_Si",unique(paramList$weather$site)),
             paste0("K_drain_Si",unique(paramList$weather$site)),
             paste0("startN_Si",unique(paramList$weather$site)),
             paste0("startC_Si",unique(paramList$weather$site)),
             "pFS2","pFS20","gammaF0","tgammaF","Rttover","mF","mR",
             "mS","Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0",
             "SWpower0","sigma_zR"
             ,"aS","nS","pRx","pRn","gammaFx",
             "SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
             "k","Qa","Qb","MaxIntcptn")
  

  names(px) <- nm_all

  #create vector of -Inf vals to return if the model run function fails
 # fail <- rep(-Inf, nrow(px))
  #Get list of site names from the paramList$weather dataframe
  siteLst <- (unique(paramList$weather$site))
  #Update dataframe with each chains proposed param values repeated for each site
  siteVec <- rep(siteLst, each = nrow(px))
  px$chain <-
    rep(1:nrow(px))  #add chain number for reference and aggregating later
  px <- px[rep(seq_len(nrow(px)), length(siteLst)),]
  px$site <- siteVec
  
  
  
  #split px dataframe into list for running in parallel lapply function
  splitParams <- split(px, seq(nrow(px)))

  #update px dataframe with likelihood values for each site
  if(Sys.info()[1]!="Windows")
    { px$ll <-
    do.call(rbind, mcmapply(runMod_mixedP,splitParams,MoreArgs = list(paramList,nm_all),SIMPLIFY = F, mc.cores = 15))
  }else{
    px$ll <-
      do.call(rbind, mapply(runMod_mixedP,splitParams,MoreArgs = list(paramList,nm_all),SIMPLIFY = F))
  }
  #sum likelihood values by chain and return vector of likelihood values, one for each chain being run

   NlogLik <-
    as.vector(px %>% group_by(chain) %>% summarise(ll = sum(ll)) %>%
                pull(ll))

  
  
  NlogLik[is.na(NlogLik)==T]<--Inf
  return(NlogLik)
  
}





#runMod function
#'@param newParams list element containing dataframe of proposed parameter values for single chain for single site
runMod_mixedP <- function(newParams,paramListX,nm_all) {
  res <- tryCatch({
    
    #update parameter list with site specific params proposals
    siteSpecNm<-c(nm_all[grepl(paste0("\\_Si",newParams$site), nm_all)],
                  nm_all[!grepl("\\_Si", nm_all)])
    fitNm<-sub("_Si.*", "",siteSpecNm)
    
    
    paramListX[fitNm] <- newParams[siteSpecNm]
    #filter paramList$weather dataframe by site
    paramListX$weather <- filter(paramListX$weather, site == newParams$site)
    #get observed values (also contained in paramListX$weather dataframe) for site being fitted
    observed <-
      filter(paramListX$weather, is.na(mean_dbh_cm) == F) %>% group_by(Year) %>% summarise(dbh =
                                                                                             median(mean_dbh_cm),
                                                                                           dbhSD = median(dbhSD_cm))
    observed_N <-paramListX$weather$StemsPerHa[1]
    
    observed$dbhSD <-
      ifelse(observed$dbhSD == 0, 0.0001, observed$dbhSD)
    sY = min(observed$Year)
    eY = max(observed$Year)
    
    #filter climate data by planting year so model runs from planting year
    paramListX$weather <-
      filter(paramListX$weather, Year >= paramListX$weather$plantingYear[1])
    
    #run model
    output <- do.call(fr3PGDN, paramListX)
    
    #filter simulated data to match observed data format
    modelled <-
      output %>% group_by(Year) %>% summarise(dg = mean(dg)) %>% filter(Year >=
                                                                          sY & Year <= eY)
    modelled_N <-tail(output$N,1)
      
    
    #run likelihood function of observed vs simulated and get liklelihood value
    ifelse(
      any(is.na(modelled) == T),
      -Inf,
      flogL(
        data = c(observed$dbh,observed_N),
        sims = c(modelled$dg,modelled_N),
        data_s = c(observed$dbhSD,modelled_N*0.1)
      )
    )
  },
  error = function(cond) {
    #return -Inf if something goes wrong with param proposals
    return(-Inf)
  })
  
  res<-ifelse(max(output$LAI)>10,-Inf,res)
  return(res)
}


