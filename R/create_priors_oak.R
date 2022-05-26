#' @export
createPriors_oak<-function(oak,sd=F){
  
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
  
  pValues <- as.vector(unlist(oak[nm]))
  
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
