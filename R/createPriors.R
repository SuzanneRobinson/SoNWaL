
createPriors_sitka<-function(sitka,sd=F){
  
  nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","E_S1","E_S2","shared_area","maxRootDepth","K_drain",
        "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
        "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
        "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0","SWpower0","Qa","Qb","MaxIntcptn","k","startN","startC","Q10","Q10X","Topt","Nleach_r")
  
  
  
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
    0.01, #aS
    2, #nS
    0.4, #pRx
    0.2, #pRn
    0.02, #gammaFx
    0.0005, #gammaF0
    40, #tgammaF
    0.003, #Rttover
    0.1, #mF
    0.1, #mR
    0.1, #mS
    4, #SLA0
    3, #SLA1
    4, #tSLA
    0.03, #alpha
    0.43, #Y
    0.001, #m0
    0.02, #MaxCond
    2.5, #LAIgcx
    0.01, #CoeffCond
    0.1, #BLcond
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
    0.1,
    2,
    -95,
    0.7,
    0.05,
    0.4,
    1,
    10,
  0.01,
  1.2,
  5,#Topt
  0.01)
  
  f.increase <-
    c(
      0.26,#wiltPoint
      0.33,#fieldCap
      0.6,#satPoint
      1,#K_s
      4,#V_nr
      1,#sigma_zR
      2,#E_S1
      2,#E_S2
      6, #shared_area
      1.5, #maxRootDepth
      1, #K_drain
      1,#pFS2
      1,#pFS20
      0.25,#aS
      3,#nS
      1,#pRx
      0.4,#pRn
      0.03,#gammaFx
      0.0025,#gammaF0
      110,#tgammaF
      0.16,#Rttover
      0.5,#mF
      0.5,#mR
      0.5,#mS
      9,#SLA0
      6,#SLA1
      8,#tSLA
      0.06,#alpha
      0.49,#Y
      1,#m0
      0.03,#MaxCond
      4,#LAIgcx
      0.2,#CoeffCond
      0.4,#BLcond
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
      0.6,
      6,
      -5,
      0.9,
      0.2,
      0.6,
      15,
      500,
      1,
      6,
      35,
      1
    )
  
  ##Need to check what priors we are using!

  pValues <- as.vector(unlist(sitka[nm]))
  
  pMaxima <- f.increase*1.5
  pMinima<- f.decrease*0.5
  
  pMaxima[[30]]<-0.5
  pMinima[[31]]<-0.01
  
  pMaxima[[48]]<-0.65
  pMinima[[48]]<-0.2
  
  pMaxima[[49]]<-7
  pMinima[[49]]<-2
  
  pMaxima[[34]]<-0.3
  pMinima[[34]]<-0.01
  
  sdVals<-(pMaxima-pMinima)*0.8
  #some manual adjustments
  sdVals[5]<-2
  sdVals[7]<-15
  sdVals[8]<-15
  pValues[7]<-10
  pValues[8]<-10
  
  sdVals[[48]]<-0.3
  sdVals[[49]]<-2
  
  sdVals[[50]]<-20
  sdVals[[51]]<-0.05
  pMaxima[[50]]<--5
  pMinima[[50]]<--95
  
  pMaxima[[51]]<-0.9
  pMinima[[51]]<-0.7
  
  pMaxima[[52]]<-0.2
  pMinima[[52]]<-0.05

  sdVals[[52]]<-0.01

  pMaxima[[33]]<-0.3
  pMinima[[33]]<-0.01
  
  sdVals[[21]]<-0.02
  
  sdVals[[53]]<-0.03
  
  pMaxima[[54]]<-25
  pMinima[[54]]<-1
  sdVals[[54]]<-10
  
  pMaxima[[55]]<-1500
  pMinima[[55]]<-50
  sdVals[[55]]<-800
  
  pMaxima[[56]]<-1
  pMinima[[56]]<-0.01
  sdVals[[56]]<-0.5
  
  pMaxima[[57]]<-5
  pMinima[[57]]<-1.2
  sdVals[[57]]<-3
  
  sdVals[[58]]<-10
  pMaxima[[58]]<-45
  pMinima[[58]]<-1
  
  sdVals[[59]]<-1
  pMaxima[[59]]<-1
  pMinima[[59]]<-0.01
  

  pMaxima[[21]]<-0.2
  pMinima[[21]]<-0.001
  sdVals[[21]]<-.2
 
  pValues<-pValues[-c(7,8)]
  sdVals<-sdVals[-c(7,8)]
  pMaxima<-pMaxima[-c(7,8)]
  pMinima<-pMinima[-c(7,8)]
  
  priorVals <- createTruncatedNormalPrior(mean = pValues, sd=sdVals,
                                          lower = pMinima, upper = pMaxima)
  
  ifelse(sd==F,return(priorVals),return(sdVals))
  
}



createPriors_pine<-function(pine,sd=F){
  
  
  nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","E_S1","E_S2","shared_area","maxRootDepth","K_drain",
        "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
        "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
        "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0","SWpower0","Qa","Qb","MaxIntcptn","k","startN","startC")
  
  
  
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
    0.01, #aS
    2, #nS
    0.4, #pRx
    0.2, #pRn
    0.02, #gammaFx
    0.0005, #gammaF0
    40, #tgammaF
    0.003, #Rttover
    0.1, #mF
    0.1, #mR
    0.1, #mS
    4, #SLA0
    3, #SLA1
    4, #tSLA
    0.03, #alpha
    0.43, #Y
    0.001, #m0
    0.02, #MaxCond
    2.5, #LAIgcx
    0.01, #CoeffCond
    0.1, #BLcond
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
    0.1,
    2,
    -95,
    0.7,
    0.05,
    0.6,
    1,
    10
  )
  
  f.increase <-
    c(
      0.26,#wiltPoint
      0.33,#fieldCap
      0.6,#satPoint
      1,#K_s
      4,#V_nr
      1,#sigma_zR
      2,#E_S1
      2,#E_S2
      6, #shared_area
      1.5, #maxRootDepth
      1, #K_drain
      1,#pFS2
      1,#pFS20
      0.25,#aS
      3,#nS
      1,#pRx
      0.4,#pRn
      0.03,#gammaFx
      0.0025,#gammaF0
      110,#tgammaF
      0.16,#Rttover
      0.5,#mF
      0.5,#mR
      0.5,#mS
      9,#SLA0
      6,#SLA1
      8,#tSLA
      0.06,#alpha
      0.49,#Y
      1,#m0
      0.03,#MaxCond
      4,#LAIgcx
      0.2,#CoeffCond
      0.4,#BLcond
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
      0.6,
      6,
      -5,
      0.9,
      0.2,
      0.4,
      15,
      500
    )
  
  pValues <- as.vector(unlist(sitka[nm]))
  
  pMaxima <- f.increase*1.5
  pMinima<- f.decrease*0.5
  
  pMaxima[[30]]<-0.5
  pMinima[[31]]<-0.01
  
  pMaxima[[48]]<-0.65
  pMinima[[48]]<-0.2
  
  pMaxima[[49]]<-7
  pMinima[[49]]<-2
  
  pMaxima[[34]]<-0.3
  pMinima[[34]]<-0.01
  
  sdVals<-(pMaxima-pMinima)*0.8
  #some manual adjustments
  sdVals[5]<-2
  sdVals[7]<-15
  sdVals[8]<-15
  pValues[7]<-10
  pValues[8]<-10
  
  sdVals[[48]]<-0.3
  sdVals[[49]]<-2
  
  sdVals[[50]]<-20
  sdVals[[51]]<-0.05
  pMaxima[[50]]<--5
  pMinima[[50]]<--95
  
  pMaxima[[51]]<-0.9
  pMinima[[51]]<-0.7
  
  pMaxima[[52]]<-0.2
  pMinima[[52]]<-0.05
  
  sdVals[[52]]<-0.01
  
  pMaxima[[33]]<-0.3
  pMinima[[33]]<-0.01
  
  sdVals[[21]]<-0.02
  
  sdVals[[53]]<-0.03
  
  
  pMaxima[[21]]<-0.2
  pMinima[[21]]<-0.001
  sdVals[[21]]<-.2
  
  pValues<-pValues[-c(7,8)]
  sdVals<-sdVals[-c(7,8)]
  pMaxima<-pMaxima[-c(7,8)]
  pMinima<-pMinima[-c(7,8)]
  
  priorVals <- createTruncatedNormalPrior(mean = pValues, sd=sdVals,
                                          lower = pMinima, upper = pMaxima)
  
  
  ifelse(sd==F,return(priorVals),return(sdVals))
  
}





createPriors_Cont<-function(out){
  
  nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","E_S1","E_S2","shared_area","maxRootDepth","K_drain",
        "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover")
  
  library(matrixStats)
  nmc = nrow(out$chain[[1]])
  outSample   <- getSample(out,start=nmc/1.2,thin=5)
  
  sdVals<-t(as.data.frame(colSds(outSample)))
  names(sdVals)<-names(as.data.frame(outSample))
  sdVals<-sdVals[nm]
  
  codM<-miscTools::colMedians(as.data.frame(outSample))
  codM<-codM[nm]
  

  ##Need to check what priors we are using!
  pMaxima <- as.vector(codM+5*sdVals)
  pMinima <- as.vector(codM-5*sdVals)

  priorVals <- createTruncatedNormalPrior(mean = codM, sd=sdVals,
                                          lower = pMinima, upper = pMaxima)
  
  return(priorVals)
  
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





createPriors_sitka_rg_all_old<-function(sitka,sd=F,out){
  codM<-as.data.frame(mergeChains(out$chain))[,-c(54:56)]
  colnames(codM)<- c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","shared_area","maxRootDepth","K_drain",
                       "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
                       "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
                       "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0","SWpower0","Qa","Qb","MaxIntcptn","k","startN","startC")
  
  
  nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","E_S1","E_S2","shared_area","maxRootDepth","K_drain",
        "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
        "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
        "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0","SWpower0","Qa","Qb","MaxIntcptn","k","startN","startC","Q10","Q10X","Topt","Nleach_r","poorSoilMod")

  
  
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
    0.01, #aS
    2, #nS
    0.4, #pRx
    0.2, #pRn
    0.02, #gammaFx
    0.0005, #gammaF0
    40, #tgammaF
    0.003, #Rttover
    0.1, #mF
    0.1, #mR
    0.1, #mS
    4, #SLA0
    3, #SLA1
    4, #tSLA
    0.03, #alpha
    0.43, #Y
    0.001, #m0
    0.02, #MaxCond
    2.5, #LAIgcx
    0.01, #CoeffCond
    0.1, #BLcond
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
    0.1,
    2,
    -95,
    0.7,
    0.05,
    0.4,
    1,
    10,
    0.01,
    1.2,
    5,#Topt
    0.01,
    0.1#poorSoilMod
    )
  
  f.increase <-
    c(
      0.26,#wiltPoint
      0.33,#fieldCap
      0.6,#satPoint
      1,#K_s
      4,#V_nr
      3,#sigma_zR
      2,#E_S1
      2,#E_S2
      6, #shared_area
      1.5, #maxRootDepth
      1, #K_drain
      1,#pFS2
      1,#pFS20
      0.25,#aS
      3,#nS
      1,#pRx
      0.4,#pRn
      0.03,#gammaFx
      0.0025,#gammaF0
      110,#tgammaF
      0.16,#Rttover
      0.5,#mF
      0.5,#mR
      0.5,#mS
      9,#SLA0
      6,#SLA1
      8,#tSLA
      0.06,#alpha
      0.49,#Y
      1,#m0
      0.03,#MaxCond
      4,#LAIgcx
      0.2,#CoeffCond
      0.4,#BLcond
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
      0.6,
      6,
      -5,
      0.9,
      0.2,
      0.6,
      15,
      500,
      1,
      6,
      35,
      1,
      0.95
    )
  
  ##Need to check what priors we are using!
  
  pValues <- as.vector(unlist(sitka[nm]))
  
  pMaxima <- f.increase*1.5
  pMinima<- f.decrease*0.5
  
  pMaxima[[30]]<-0.5
  pMinima[[31]]<-0.01
  
  pMaxima[[48]]<-0.65
  pMinima[[48]]<-0.2
  
  pMaxima[[49]]<-7
  pMinima[[49]]<-2
  
  pMaxima[[34]]<-0.3
  pMinima[[34]]<-0.01
  
  sdVals<-(pMaxima-pMinima)*0.8
  
  #some manual adjustments
  sdVals[5]<-2
  sdVals[7]<-15
  sdVals[8]<-15
  pValues[7]<-10
  pValues[8]<-10
  
  sdVals[[48]]<-0.3
  sdVals[[49]]<-2
  
  sdVals[[50]]<-20
  sdVals[[51]]<-0.05
  pMaxima[[50]]<--5
  pMinima[[50]]<--95
  
  pMaxima[[51]]<-0.9
  pMinima[[51]]<-0.7
  
  pMaxima[[52]]<-0.2
  pMinima[[52]]<-0.05
  
  sdVals[[52]]<-0.01
  
  pMaxima[[33]]<-0.3
  pMinima[[33]]<-0.01
  
  sdVals[[21]]<-0.02
  
  sdVals[[53]]<-0.03
  
  pMaxima[[54]]<-25
  pMinima[[54]]<-1
  sdVals[[54]]<-10
  
  pMaxima[[55]]<-1500
  pMinima[[55]]<-50
  sdVals[[55]]<-800
  
  pMaxima[[56]]<-1
  pMinima[[56]]<-0.01
  sdVals[[56]]<-0.5
  
  pMaxima[[57]]<-5
  pMinima[[57]]<-1.2
  sdVals[[57]]<-3
  
  sdVals[[58]]<-10
  pMaxima[[58]]<-45
  pMinima[[58]]<-1
  
  sdVals[[59]]<-1
  pMaxima[[59]]<-1
  pMinima[[59]]<-0.01
  
  
  pMaxima[[21]]<-0.2
  pMinima[[21]]<-0.001
  sdVals[[21]]<-.2
  
  pValues<-pValues[-c(1:5,7,8,10,11,48,49,56:59)]
  
  
 # sdVals<-1.5*colSds(as.matrix(codM[c("sigma_zR","shared_area",
 #                                 "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
#                                  "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
#                                  "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","Qa","Qb","MaxIntcptn","k","startN","startC")]))
  
  sdVals[45]<-0.5

  pMaxima<-pMaxima[-c(1:5,7,8,10,11,48,49,56:59)]
  pMinima<-pMinima[-c(1:5,7,8,10,11,48,49,56:59)]
  pMaxima[45]<-0.95
  pMinima[45]<-0.1
  pValues[45]<-0.8
  
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
