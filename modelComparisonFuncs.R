##model comparisons for sitka


bicFunc<-function(n,k,ll){
  k*log(n)-2*(ll)
}


#bicFunc(n=294,k=36,ll=-7425.341)

sitkaBICcomp<-function(obsData=data,startY=2018,endY=2018,pNum=47,params=sitka,timeStep="monthly"){
  
  obsData<- filter(obsData,year>=startY & year<=endY)
  
  
  if(timeStep =="weekly"){
    obsData<-obsData%>%
      mutate(grp=week(as.Date(yday, origin = paste0(year,"-01-01"))))
  }
  
  if(timeStep =="monthly"){
    obsData<-obsData%>%
      mutate(grp=month(as.Date(yday, origin = paste0(year,"-01-01"))))
  }
  
  
  obsSitka <- c(pull(obsData%>% 
                       group_by(year,grp) %>%
                       dplyr::summarise(gpp=mean(gpp))%>%
                       select(gpp)),                ## GPP
                pull(obsData%>%
                       group_by(year,grp) %>%
                       dplyr::summarise(npp=mean(npp))%>%
                       select(npp)),                ## NPP
                pull(obsData%>%
                       group_by(year,grp) %>%
                       dplyr::summarise(nee=mean(nee))%>%
                       select(nee)),                ## NEE
                pull(obsData%>%
                       group_by(year,grp) %>%
                       dplyr::summarise(reco=mean(reco))%>%
                       select(reco)),           ## Rs
                pull(obsData%>%
                       group_by(year,grp) %>%
                       dplyr::summarise(et=sum(et))%>%
                       select(et))             ## Etransp
                
  )
  
  
  devSitka <- c(rep(.3,length(pull(obsData%>%
                                group_by(year,grp) %>%
                                dplyr::summarise(gpp=mean(gpp))%>%
                                select(gpp)))),
           rep(.3,length( pull(obsData%>%
                                 group_by(year,grp) %>%
                                 dplyr::summarise(npp=mean(npp))%>%
                                 select(npp)))),
           rep(.3,length(pull(obsData%>%
                                group_by(year,grp) %>%
                                dplyr::summarise(nee=mean(nee))%>%
                                select(nee)))),
           rep(.3,length(pull(obsData%>%
                                group_by(year,grp) %>%
                                dplyr::summarise(reco=mean(reco))%>%
                                select(reco)))),
           rep(4,length(pull(obsData%>%
                               group_by(year,grp) %>%
                               dplyr::summarise(et=sum(et))%>%
                               select(et))))
           
  )
  
  

## Extract simulated data for use in likelihood function
sampleOutputSitka<-function(df,sY,eY){
  
  if(nrow(df)<600){
    df<-df%>%
      mutate(GPP=GPP*100/days_in_month(Month))%>%
      mutate(NPP=NPP*100/days_in_month(Month))%>%
      mutate(NEE=NEE*100/days_in_month(Month))%>%
      mutate(Reco=Reco*100/days_in_month(Month))
  } else {
    df<-df%>%
      mutate(GPP=GPP*100/7)%>%
      mutate(NPP=NPP*100/7)%>%
      mutate(NEE=NEE*100/7)%>%
      mutate(Reco=Reco*100/7)
  }
  
  df<- filter(df,Year>=sY&Year<=eY)
  
  m<-c(aggregate(df$GPP~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$NPP~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$NEE~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$Reco~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$Etransp~ df$Month+df$Year,FUN=sum)[,3]
  )
  m
  return(m)
}

output<-do.call(fr3PGDN,params)
modelled <-sampleOutputSitka(output,startY,endY)
ll<-ifelse(any(is.na(modelled)==T),-Inf,sum(dnorm(x=obsSitka,sd =devSitka, mean=modelled,log=T),na.rm = T))
print(paste0("Log likelihood = ",ll))
bicFunc(n=342,pNum,ll)

}

sitkaBICcomp(obsData = flxdata_daily,startY=2015,endY=2018,pNum=36,params=sitka,timeStep="monthly")


#bicFunc(n=294,k=36,ll=-7425.341)

pineBICcomp<-function(obsData=data,startY=2018,endY=2018,pNum=47,params=sitka){
  
  observedPine <- c(GPP$GPP,                ## GPP - monthly avg
                    NEE$NEE,                ## NPP - monthly avg
                    reco$reco,              ## NEE - monthly avg
                    LAI[,2],                    ##LAI -yearly average 
                    dbh[,2],                    ## DBH
                    totC,                   ## totC, 1995-1996
                    totN                    ## totN, 1995-1996
  )
  
  
  startYear<-1996
  endYear<-2014
  
  devPine <- c(rep(.3,nrow(dplyr::filter(GPP,year>=startYear&year<=endYear&year!=2007))),
               rep(.3,nrow(dplyr::filter(NEE,year>=startYear&year<=endYear&year!=2007))),
               rep(.3,nrow(dplyr::filter(reco,year>=startYear&year<=endYear&year!=2007))),
               rep(.1,nrow(LAI)),
               rep(.3,nrow(dbh)),
               20,
               5
               
  )
  
  sampleOutputPine<-function(df,sY,eY){
    m<-c(
      pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
             group_by(Year,Month)%>%
             summarise(sum=sum(GPP,na.rm=TRUE))%>%
             select(sum)),
      pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
             group_by(Year,Month)%>%
             summarise(sum=sum(NEE,na.rm=TRUE))%>%
             select(sum)),
      pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
             group_by(Year,Month)%>%
             summarise(sum=sum(Reco,na.rm=TRUE))%>%
             select(sum)),
      
      pull(filter(df,Year>=1995&Year<=2011)%>%
             group_by(Year)%>%
             summarise(mean=mean(LAI,na.rm=TRUE))%>%
             select(mean)),
      
      pull(filter(df,Year>=1995&Year<=2011)%>%
             group_by(Year)%>%
             summarise(mean=mean(dg,na.rm=TRUE))%>%
             select(mean)),
      
      pull(filter(df,Year==1996)%>%
             group_by(Year)%>%
             summarise(mean=mean(totC,na.rm=TRUE))%>%
             select(mean)),
      
      pull(filter(df,Year==1996)%>%
             group_by(Year)%>%
             summarise(mean=mean(totN,na.rm=TRUE))%>%
             select(mean))
      
      
    )
    return(m)
  }
  
  output<-do.call(fr3PGDN,params)
  modelled <-sampleOutputPine(output,startY,endY)
  ll<-ifelse(any(is.na(modelled)==T),-Inf,sum(dnorm(x=obsSitka,sd =devSitka, mean=modelled,log=T),na.rm = T))
  print(paste0("Log likelihood = ",ll))
  bicFunc(n=length(obsSitka),pNum,ll)
  
}


