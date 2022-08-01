
# Extract simulated data for use in likelihood function for shorter time-steps
# monthly data is in tons per hectare per month which matches 3pg output so no need to convert
sampleOutput<-function(sim,sY,eY,swc=T){
  #convert to average grams per m2 per day depending on timestep of model (using output length to quickly get time-step)
  modif<- if(nrow(sim)<1000) 1.6 else  7.142857
  nDays<- if(nrow(sim)<1000) 30 else  7
  
  sim<- filter(sim,Year>=sY&Year<=eY)
  
  m<-sim%>%
    mutate(GPP=GPP*modif, NPP=NPP*modif, NEE=NEE*modif, Rs=Rs*modif, Reco=Reco*modif, swc = volSWC_rz, EvapTransp = EvapTransp/nDays)%>%
    group_by(Year, Month) %>%
    summarise_all(mean)
  
  m<-c( m$Rs, m$GPP, m$NEE, m$EvapTransp,
        filter(sim,Year==2015&Month==8)$LAI[1],
        filter(sim,Year==2018&Month==8)$LAI[1],
        filter(sim,Year==2018&Month==8)$N[1],
        filter(sim,Year==2018&Month==8)$dg[1],
        filter(sim,Year==2015&Month==7)$totC[1],
        filter(sim,Year==2015&Month==7)$totN[1],
        m$swc)
  
  return(m)
}




# sivia likelihood calculation
flogL <- function(sims,data,data_s)
{ 
  Ri         <- (sims - data) / data_s
  i0         <- which( abs(Ri)<1.e-08 )
  
  logLi      <- log(1-exp(-0.5*Ri^2)) - log(Ri^2) - 0.5*log(2*pi) - log(data_s)
  logLi[i0]  <- -0.5*log(2*pi) - log(2*data_s[i0])
  
  sum(logLi)
}



## Likelihood function
LL_sitka<- function(p){
  p<-p*.GlobalEnv$param_scaler
  sitka[.GlobalEnv$nm]<-p
  
   NlogLik <- tryCatch(
    {
      output<-   do.call(SoNWaL,sitka)
      modelled <-sampleOutput(output,.GlobalEnv$startYear,.GlobalEnv$endYear,swc=sitka$waterBalanceSubMods)
      NlogLik  <-   ifelse(any(is.na(modelled)==T),-Inf,flogL(data=.GlobalEnv$observed,sims=modelled,data_s=.GlobalEnv$dev))
     # NlogLik  <-   ifelse(any(is.na(modelled)==T),-Inf,sum(dnorm(modelled, .GlobalEnv$observed, .GlobalEnv$dev,log=T)))
      
      NlogLik<-ifelse(max(output$LAI)>10,-Inf,NlogLik)
      NlogLik<-ifelse(min(output$totN)<1,-Inf,NlogLik)
      NlogLik<-ifelse(sitka$fieldCap<sitka$wiltPoint,-Inf,NlogLik)
      NlogLik<-ifelse(sitka$satPoint<sitka$fieldCap,-Inf,NlogLik)
    },
    error=function(cond) {
      return(-Inf)
    })
  return(NlogLik)
}


