##model comparisons for sitka


bicFunc<-function(n,k,ll){
  k*log(n)-2*(ll)
}


#bicFunc(n=294,k=36,ll=-7425.341)

sitkaBICcomp<-function(obsData=data,startY=2018,endY=2018,pNum=47,params=sitka){
obsSitka <- c(
  pull(obsData%>%filter(year>=startY&year<=endY)%>%select(gpp)),                ## GPP
  pull(obsData%>%filter(year>=startY&year<=endY)%>%select(nee)),                ## NEE
  pull(obsData%>%filter(year>=startY&year<=endY)%>%select(reco)),               ## Reco
  pull(obsData%>%filter(year>=startY&year<=endY)%>%select(et))                 ## Etransp
          
)

devSitka <- c(rep(.3,nrow(dplyr::filter(obsData,year>=startY&year<=endY))),
         rep(.3,nrow(dplyr::filter(obsData,year>=startY&year<=endY))),
         rep(.3,nrow(dplyr::filter(obsData,year>=startY&year<=endY))),
         rep(6,nrow(dplyr::filter(obsData,year>=startY&year<=endY)))
         
)


## Extract simulated data for use in likelihood function
sampleOutputSitka<-function(df,sY,eY){
  df<- filter(df,Year>=sY&Year<=eY)
  m<-c(aggregate(df$GPP~ df$Month+df$Year,FUN=sum)[,3],
       aggregate(df$NEE~ df$Month+df$Year,FUN=sum)[,3],
       aggregate(df$Reco~ df$Month+df$Year,FUN=sum)[,3],
       aggregate(df$Etransp~ df$Month+df$Year,FUN=mean)[,3]
       
  )
  return(m)
}

output<-do.call(fr3PGDN,params)
modelled <-sampleOutputSitka(output,startY,endY)
ll<-ifelse(any(is.na(modelled)==T),-Inf,sum(dnorm(x=obsSitka,sd =devSitka, mean=modelled,log=T),na.rm = T))
print(paste0("Log likelihood = ",ll))
bicFunc(n=length(obsSitka),pNum,ll)

}

