
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
    0.2, #sigma_zR
    0.0001, #E_S1
    0.001, #E_S2
    2, #shared_area
    0.3, #maxRootDepth
    0.05, #K_drain
    0.588503613257886, #pFS2
    0.752929538228874, #pFS20
    0.956131627577964, #aS
    0.050456035523466, #nS
    0.384021499609213, #pRx
    0.250229439327847, #pRn
    0.57408236899746, #gammaFx
    0.909666760291794, #gammaF0
    0.853276910139941, #tgammaF
    0.974961101217424, #Rttover
    1, #mF
    0.636422367959785, #mR
    0.732916669791679, #mS
    0.443930919848964, #SLA0
    0.741758519667562, #SLA1
    0.816463641720414, #tSLA
    0.221779786451702, #alpha
    0.303779963365252, #Y
    1, #m0
    0.00141038795075, #MaxCond
    0.730688961031379, #LAIgcx
    0.899808741360758, #CoeffCond
    0.024817372196732, #BLcond
    0.99632339563598, #Nf
    0.996373181003088, #Navm
    0.999649942946159, #Navx
    0.996388219783102, #klmax
    0.998203040988276, #krmax
    0.998245174258832, #komax
    0.97983098579238, #hc
    0.913069476259938, #qir
    0.961283723717706, #qil
    0.950056672692535, #qh
    0.893875965852296, #qbc
    0.991080780202615, #el
    0.990457295759556, #er
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
      3,#sigma_zR
      2,#E_S1
      2,#E_S2
      6, #shared_area
      4, #maxRootDepth
      1, #K_drain
      0.573973679288588,#pFS2
      0.235352308855631,#pFS20
      1.86098081013281,#aS
      0.374136113325978,#nS
      0.231957000781575,#pRx
      0.56202200140032,#pRn
      3.45793787115991,#gammaFx
      1.30349761255926,#gammaF0
      0.600615525746093,#tgammaF
      0.251944939128821,#Rttover
      0.768680943537667,#mF
      0.817888160201076,#mR
      0.335416651041606,#mS
      0.668207240453109,#SLA0
      0.549448881994627,#SLA1
      0.835363582795864,#tSLA
      0.03762695139773,#alpha
      0.218385064110809,#Y
      0.917925202458998,#m0
      2.5949226033773,#MaxCond
      1.15448831174897,#LAIgcx
      0.001912586392424,#CoeffCond
      5.82627839462287,#BLcond
      6.35320872803933,#Nf
      2.62681899691161,#Navm
      2.50057053840858,#Navx
      2.61178021689853,#klmax
      0.796959011723578,#krmax
      0.754825741168422,#komax
      1.01690142076198,#hc
      0.738610474801243,#qir
      0.935813814114711,#qil
      0.498299819223935,#qh
      1.12248068295408,#qbc
      0.783843959477034,#el
      0.90854084808886, #er
      0.6,
      6,
      -5,
      0.9,
      0.2,
      0.4,
      15,
      500
    )
  
  ##Need to check what priors we are using!
  pMaxima <- as.vector(unlist(pine[nm])*(1+(f.increase)))
  pMinima <- as.vector(unlist(pine[nm])*(1-(f.decrease)))
  pValues <- as.vector(unlist(pine[nm]))
  
  pMaxima[1:11] <- f.increase[1:11]
  pMinima[1:11] <- f.decrease[1:11]
  
  pMaxima[[30]]<-0.5
  pMinima[[31]]<-0.01
  
  pMaxima[[48]]<-0.65
  pMinima[[48]]<-0.2
  
  pMaxima[[49]]<-7
  pMinima[[49]]<-2
  
  pMaxima[[34]]<-0.3
  pMinima[[34]]<-0.01
  
  sdVals<-(pMaxima-pMinima)*0.6
  #ES1 and 2 on very different scales so set manually
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
  
  pMaxima[[54]]<-350
  pMinima[[54]]<-1
  sdVals[[54]]<-100
  pValues[[52]]<-50
  
  pMaxima[[55]]<-5000
  pMinima[[55]]<-10
  sdVals[[55]]<-2000
  pValues[[55]]<-2000
  
  
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





createPriors_sitka_rg_all<-function(sitka,sd=F,out){
  codM<-as.data.frame(mergeChains(out$chain))[,-c(54:56)]
  colnames(codM)<- c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","shared_area","maxRootDepth","K_drain",
                       "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
                       "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
                       "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0","SWpower0","Qa","Qb","MaxIntcptn","k","startN","startC")
  
  
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
  
  sdVals<-(pMaxima-pMinima)*0.1
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
  
  
  sdVals<-1.2*colSds(as.matrix(codM[c("sigma_zR","shared_area",
                                  "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
                                  "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
                                  "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","Qa","Qb","MaxIntcptn","k","startN","startC")]))
  
  
  pMaxima<-pMaxima[-c(1:5,7,8,10,11,48,49,56:59)]
  pMinima<-pMinima[-c(1:5,7,8,10,11,48,49,56:59)]
  
  priorVals <- createTruncatedNormalPrior(mean = pValues, sd=sdVals,
                                          lower = pMinima, upper = pMaxima)
  
  ifelse(sd==F,return(priorVals),return(sdVals))
  
}


