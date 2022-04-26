#' createPriors create priors for MCMC calibration for sitka
#' @param sitka current sitka parameters to base priors on
#' @return dataframe of priors
#' @export
createPriors_sitka<-function(sitka,sd=F){
  
  nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","E_S1","E_S2","shared_area","maxRootDepth","K_drain",
        "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
        "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
        "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0","SWpower0","Qa","Qb","MaxIntcptn","k","startN","startC","Q10")
  
  
  
  f.decrease <- c(
    0.06, #wiltPoint
    0.15, #fieldCap
    0.3,#satPoint
    0.05, #K_s
    0.3, #V_nr
    0.1, #sigma_zR
    0.0001, #E_S1
    0.001, #E_S2
    1, #shared_area
    0.2, #maxRootDepth
    0.05, #K_drain
    0.4, #pFS2
    0.2, #pFS20
    0.001, #aS
    2, #nS
    0.4, #pRx
    0.2, #pRn
    0.008, #gammaFx
    0.00005, #gammaF0
    20, #tgammaF
    0.001, #Rttover
    0.001, #mF
    0.1, #mR
    0.1, #mS
    4, #SLA0
    2, #SLA1
    2, #tSLA
    0.03, #alpha
    0.43, #Y
    0.001, #m0
    0.02, #MaxCond
    0.6, #LAIgcx
    0.01, #CoeffCond
    0.05, #BLcond
    0.0001, #Nf
    0.0001, #Navm
    0.5, #Navx
    0.0002, #klmax
    0.0002, #krmax
    0.00002, #komax
    0.01, #hc
    90, #qir
    1, #qil
    3, #qh
    0.6, #qbc
    0.01, #el
    0.01, #er
    0.1,
    2,
    -95,
    0.7,
    0.05,
    0.4,
    1,
    10,
    1.1) 
  
  f.increase <-
    c(
      0.26,#wiltPoint
      0.33,#fieldCap
      0.6,#satPoint
      1,#K_s
      1.5,#V_nr
      2,#sigma_zR
      1,#E_S1
      1,#E_S2
      6, #shared_area
      2, #maxRootDepth
      1, #K_drain
      1,#pFS2
      1,#pFS20
      0.25,#aS
      3,#nS
      1,#pRx
      0.4,#pRn
      0.07,#gammaFx
      0.0025,#gammaF0
      110,#tgammaF
      0.16,#Rttover
      0.8,#mF
      0.8,#mR
      0.5,#mS
      9,#SLA0
      6,#SLA1
      8,#tSLA
      0.06,#alpha
      0.6,#Y
      1,#m0
      0.05,#MaxCond
      4,#LAIgcx
      0.2,#CoeffCond
      0.4,#BLcond
      0.1,#Nf
      0.5,#Navm
      20,#Navx
      0.2,#klmax
      0.05,#krmax
      0.01,#komax
      0.8,#hc
      600,#qir
      60,#qil
      50,#qh
      20,#qbc
      0.9,#el
      0.9, #er
      0.6,
      6,
      -5,
      0.9,
      0.2,
      0.6,
      15,
      1200,
      10
      
    )
  
  ##Need to check what priors we are using!

  pValues <- as.vector(unlist(sitka[nm]))
  
  pMaxima <- f.increase*1.5
  pMinima<- f.decrease*0.5
  sdVals<-(pMaxima-pMinima)*0.8
  
  pMinima[50]<--95 #needs manual adjusting as negative value
  
  #pMaxima[56]<-1
  #pMinima[56]<-0.001
  sdVals[55]<-20
  
  
 # pValues<-pValues[-c(7,8)]
#  sdVals<-sdVals[-c(7,8)]
#  pMaxima<-pMaxima[-c(7,8)]
#  pMinima<-pMinima[-c(7,8)]
  
  
  sc<-rowMeans( abs( cbind(pMinima,pMaxima) ) )
  
  priorVals <- createTruncatedNormalPrior(mean = pValues/sc, sd=sdVals/sc,
                                          lower = pMinima/sc, upper = pMaxima/sc)
  

  
  ifelse(sd==F,return(list(priorVals,sc)),return(sdVals))
  
}



createPriors_pine<-function(pine,sd=F){
  
  
  nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","E_S1","E_S2","shared_area","maxRootDepth","K_drain",
        "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
        "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
        "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0","SWpower0","Qa","Qb","MaxIntcptn","k","startN","startC","kF")
  
  
  
  f.decrease <- c(
    0.06, #wiltPoint
    0.15, #fieldCap
    0.3,#satPoint
    0.05, #K_s
    0.3, #V_nr
    0.1, #sigma_zR
    0.0001, #E_S1
    0.001, #E_S2
    1, #shared_area
    0.2, #maxRootDepth
    0.05, #K_drain
    0.4, #pFS2
    0.2, #pFS20
    0.001, #aS
    2, #nS
    0.4, #pRx
    0.2, #pRn
    0.008, #gammaFx
    0.00005, #gammaF0
    20, #tgammaF
    0.001, #Rttover
    0.001, #mF
    0.1, #mR
    0.1, #mS
    4, #SLA0
    2, #SLA1
    2, #tSLA
    0.03, #alpha
    0.43, #Y
    0.001, #m0
    0.02, #MaxCond
    0.6, #LAIgcx
    0.01, #CoeffCond
    0.05, #BLcond
    0.0001, #Nf
    0.0001, #Navm
    0.5, #Navx
    0.0002, #klmax
    0.0002, #krmax
    0.00002, #komax
    0.01, #hc
    90, #qir
    1, #qil
    3, #qh
    0.6, #qbc
    0.01, #el
    0.01, #er
    0.1,
    2,
    -95,
    0.7,
    0.05,
    0.4,
    1,
    10,
    0)#kF
  
  f.increase <-
    c(
      0.45,#wiltPoint
      0.55,#fieldCap
      0.6,#satPoint
      1,#K_s
      1.5,#V_nr
      2,#sigma_zR
      2,#E_S1
      2,#E_S2
      6, #shared_area
      2, #maxRootDepth
      1, #K_drain
      1,#pFS2
      1,#pFS20
      0.25,#aS
      3,#nS
      1,#pRx
      0.4,#pRn
      0.07,#gammaFx
      0.0025,#gammaF0
      110,#tgammaF
      0.16,#Rttover
      0.8,#mF
      0.8,#mR
      0.5,#mS
      9,#SLA0
      6,#SLA1
      8,#tSLA
      0.06,#alpha
      0.6,#Y
      0.4,#m0
      0.05,#MaxCond
      4,#LAIgcx
      0.2,#CoeffCond
      0.4,#BLcond
      0.1,#Nf
      0.1,#Navm
      20,#Navx
      0.2,#klmax
      0.05,#krmax
      0.01,#komax
      0.8,#hc
      300,#qir
      40,#qil
      50,#qh
      20,#qbc
      0.9,#el
      0.9, #er
      0.6,
      6,
      -5,
      0.9,
      0.2,
      0.6,
      15,
      1200,
      1#kF
    )
  
  ##Need to check what priors we are using!
  
  pValues <- as.vector(unlist(sitka[nm]))
  
  pMaxima <- f.increase*1.5
  pMinima<- f.decrease*0.5
  sdVals<-(pMaxima-pMinima)*0.8
  
  pMinima[50]<--95 #needs manual adjusting as negative value
  pMaxima[56]<-1
  pMinima[56]<-0.001
  
  pValues<-pValues[-c(7,8)]
  sdVals<-sdVals[-c(7,8)]
  pMaxima<-pMaxima[-c(7,8)]
  pMinima<-pMinima[-c(7,8)]
  
  priorVals <- createTruncatedNormalPrior(mean = pValues, sd=sdVals,
                                          lower = pMinima, upper = pMaxima)
  
  
  ifelse(sd==F,return(priorVals),return(sdVals))
  
}















createPriors_sitka_rg<-function(sitka,sd=F){
  
  
  nm<-c("sigma_zR","shared_area",
        "Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","startN","startC")
  
  
  
  f.decrease <- c(
    
    0.1, #sigma_zR
    1, #shared_area
    0.001, #Navm
    4, #Navx
    0.009, #klmax
    0.001, #krmax
    0.0001, #komax
    0.1, #hc
    90, #qir
    10, #qil
    10, #qh
    1, #qbc
    0.1, #el
    0.1, #er
    1,
    10)
  
  f.increase <-
    c(
      3,#sigma_zR
      6, #shared_area
      0.5,#Navm
      15,#Navx
      0.1,#klmax
      0.05,#krmax
      0.01,#komax
      0.4,#hc
      600,#qir
      60,#qil
      50,#qh
      20,#qbc
      0.6,#el
      0.6, #er
      15,
      500
    )
  
  ##Need to check what priors we are using!
  
  pValues <- as.vector(unlist(sitka[nm]))
  
  pMaxima <- f.increase*1.5
  pMinima<- f.decrease*0.5
  
  
  sdVals<-(pMaxima-pMinima)*0.3
  
  
  priorVals <- createTruncatedNormalPrior(mean = pValues, sd=sdVals,
                                          lower = pMinima, upper = pMaxima)
  
  ifelse(sd==F,return(priorVals),return(sdVals))
  
}



createPriors_sitka_rg_all<-function(paramList,sd=F,originalFit,makePlots=F){
  
  #firstly get the posterior values from the original harwood calibrations ("originalFit" is the harwood calibration)
  codM<-as.data.frame(mergeChains(originalFit$chain))[,-c(54:56)] #remove the likelihood values
  #add param names to dataframe
  colnames(codM)<- c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","shared_area","maxRootDepth","K_drain",
                     "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
                     "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
                     "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0","SWpower0","Qa","Qb","MaxIntcptn","k","startN","startC")

  
  #name the "tight" species specific priors - these will be based on the posteriors of the original harwood calibration
  nm_species<-c("sigma_zR"
                ,"aS","nS","pRx","pRn","gammaFx",
                "SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
                "k","Qa","Qb","MaxIntcptn")
  
  
  #name priors that will be site specific or have wider distributions
  #first set (ones with paste0) as site specific, essentially these are repeated and site names applied (i.e. there will be 13 wiltPoints_Si[siteName])
  nm_site<-c(paste0("wiltPoint_Si",unique(paramList$weather$site)),
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
             #the next set are not site specific priors but are considered not to be species specific, so a single value is fit for all sites but with a 
             #wide prior distribution
             "pFS2","pFS20","gammaF0","tgammaF","Rttover","mF","mR",
             "mS","Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0",
             "SWpower0")
  
  #get a vector with the names for all priors
  nmAll<-c(nm_site,nm_species)
  
  
  #define lower limit of priors for truncated gaussian distributions
  f.decrease_site <- c(
    #site specific priors
    rep(0.06,length(unique(paramList$weather$site))), #wiltPoint
    rep(0.1,length(unique(paramList$weather$site))), #fieldCap
    rep(0.2,length(unique(paramList$weather$site))),#satPoint
    rep(0.05,length(unique(paramList$weather$site))), #K_s
    rep(0.3,length(unique(paramList$weather$site))), #V_nr
    rep(0.03,length(unique(paramList$weather$site))), #E_S1
    rep(0.2,length(unique(paramList$weather$site))), #E_S2
    rep(1,length(unique(paramList$weather$site))), #shared_area
    rep(0.2,length(unique(paramList$weather$site))), #maxRootDepth
    rep(0.05,length(unique(paramList$weather$site))), #K_drain
    rep(3,length(unique(paramList$weather$site))), #startN
    rep(50,length(unique(paramList$weather$site))), #startC
    
    #wide distribution priors
    0.4, #pFS2
    0.2, #pFS20
    0.0005, #gammaF0
    40, #tgammaF
    0.003, #Rttover
    0.1, #mF
    0.1, #mR
    0.1, #mS
    
    #species specific priors
    0.001, #Nf
    0.001, #Navm
    4, #Navx
    0.009, #klmax
    0.001, #krmax
    0.0001, #komax
    0.1, #hc
    90, #qir
    10, #qil
    10, #qh
    1, #qbc
    0.1, #el
    0.1, #er
    0.1,#SWconst0
    2#SWpower0
    
  )
  
  #define upper limit of priors for truncated gaussian distributions
  f.increase_site <-
    c(
      #site specific priors
      rep(0.35,length(unique(paramList$weather$site))),#wiltPoint
      rep(0.4,length(unique(paramList$weather$site))),#fieldCap
      rep(0.8,length(unique(paramList$weather$site))),#satPoint
      rep(1,length(unique(paramList$weather$site))),#K_s
      rep(4,length(unique(paramList$weather$site))),#V_nr
      rep(1,length(unique(paramList$weather$site))),#E_S1
      rep(1,length(unique(paramList$weather$site))),#E_S2
      rep(6,length(unique(paramList$weather$site))), #shared_area
      rep(2,length(unique(paramList$weather$site))), #maxRootDepth
      rep(1,length(unique(paramList$weather$site))), #K_drain
      rep(20,length(unique(paramList$weather$site))), #startN
      rep(1200,length(unique(paramList$weather$site))), #startC
      
      #wide distribution priors
      1,#pFS2
      1,#pFS20
      0.0025,#gammaF0
      110,#tgammaF
      0.16,#Rttover
      0.5,#mF
      0.5,#mR
      0.5,#mS
      
      #species specific priors
      2,#Nf
      0.5,#Navm
      15,#Navx
      0.1,#klmax
      0.05,#krmax
      0.01,#komax
      0.4,#hc
      600,#qir
      60,#qil
      50,#qh
      20,#qbc
      0.6,#el
      0.6, #er
      0.6,#SWconst0
      8#SWpower0
      
      
    )
  
  
  codM$E_S1<-0.3
  codM$E_S2<-0.6
  
  #get the mode from the original harwood fit to define the mode of the new site specific and wide dist. priors
  pValues_site<-as.vector(sapply(codM[sub("_Si.*", "", nm_site)],DescTools::Mode))
  #set upper and lower limits
  pMaxima_site <- as.vector(f.increase_site)
  pMinima_site<- as.vector(f.decrease_site)
  #set SD - currently using difference between upper and lower limits of the prior *2.5, 
  #this gives a prior specific SD value that leads to wide gaussian dist. - essentially they are almost uniform priors
  sdVals_site<-as.vector((pMaxima_site-pMinima_site)*2.5)
  #add sd for es1 and es2 as these were not originally fixed params
  sdVals_site[66:91]<-0.2
  
  #species priors
  #get mode from original cals
  pValues_species<-as.vector(sapply(codM[nm_species],DescTools::Mode))
  #use SD from original marginal posteriors
  sdVals_species<-as.vector(matrixStats::colSds(as.matrix(codM[c("sigma_zR"
                                                           ,"aS","nS","pRx","pRn","gammaFx",
                                                           "SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
                                                           "k","Qa","Qb","MaxIntcptn")])))
  
  #set upper and lower limits of the truncated pior as +/-50% (this is generally much higher than the tails of the posterior distribution, 
  #so distributions should be as in the calibration)
  #but makes life easier if all values can be put into the same trunctated normal dist creation function which requires upper and lower limits
  pMaxima_species <- as.vector(pValues_species*1.5)
  pMinima_species<- as.vector(pValues_species*0.5)
  
  #the values for this prior need swapping as negative values
  QaMax<-pMaxima_species[18]
  QaMin<-pMinima_species[18]
  pMaxima_species[18]<-QaMin
  pMinima_species[18]<-QaMax
  
  #create priors
  priorVals <- createTruncatedNormalPrior(mean = c(pValues_site,pValues_species), sd=c(sdVals_site,sdVals_species),
                                          lower = c(pMinima_site,pMinima_species), upper = c(pMaxima_site,pMaxima_species))
  
  
  
##Ignore this part - just for quickly plotting all the priors
  if(makePlots==T){
 nmAllred<-unique(sub("_Si.*", "", nmAll))
 prDatB<-data.frame(names=nmAll,mean = c(pValues_site,pValues_species), sd=c(sdVals_site,sdVals_species),
                   lower = c(pMinima_site,pMinima_species), upper = c(pMaxima_site,pMaxima_species))
 prDatB$names<-sub("_Si.*", "", prDatB$names)
 prDatB<-prDatB[!duplicated(prDatB$names),]
 
 priorValsX <- createTruncatedNormalPrior(mean = prDatB$mean, sd=prDatB$sd,
                                        lower = prDatB$lower, upper =prDatB$upper)

priPlotFunc<-function(prDatBS){
  pp<-createTruncatedNormalPrior(mean = prDatBS$mean, sd=prDatBS$sd,
                             lower = prDatBS$lower, upper =prDatBS$upper)
  
  ppDat<-data.frame(pp=(pp$sampler(10000)))
return(
  ggplot(data=ppDat,aes(pp))+
    geom_histogram(bins=50,col="black")+
    ggtitle(paste0(prDatBS$names,"10k prior samps"))
)
  
}
prDatBX<-split(prDatB,seq(nrow(prDatB)))

plots<-lapply(prDatBX,priPlotFunc)

ggarrange(plotlist=plots)
}
  
  #return prior vals (or for report markdown script tables the sd vals)
  ifelse(sd==F,return(priorVals),return(c(sdVals_site,sdVals_species)))
  
}
