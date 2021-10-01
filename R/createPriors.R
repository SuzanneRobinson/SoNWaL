
createPriors_sitka<-function(sitka){
  
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
    0.01, #E_S1
    0.01, #E_S2
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
      5,#E_S1
      5,#E_S2
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
  pMaxima <- as.vector(unlist(sitka[nm])*(1+(f.increase)))
  pMinima <- as.vector(unlist(sitka[nm])*(1-(f.decrease)))
  pValues <- as.vector(unlist(sitka[nm]))
  
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
  
  sdVals<-(pMaxima-pMinima)*0.3
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
  sdVals[[52]]<-0.05
  
  pMaxima[[33]]<-0.3
  pMinima[[33]]<-0.01
  
  sdVals[[21]]<-0.02
  
  sdVals[[53]]<-0.03
  
  pMaxima[[54]]<-15
  pMinima[[54]]<-1
  sdVals[[54]]<-4
  
  pMaxima[[55]]<-500
  pMinima[[55]]<-10
  sdVals[[55]]<-200
  
  
  priorVals <- createTruncatedNormalPrior(mean = pValues, sd=sdVals,
                                          lower = pMinima, upper = pMaxima)
  
  return(priorVals)
  
}



createPriors_pine<-function(pine){
  
  
  nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","E_S1","E_S2","shared_area","maxRootDepth","K_drain",
        "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
        "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
        "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0","SWpower0","Qa","Qb","MaxIntcptn","k","startN","startC")
  
  
  
  
  f.decrease <- c(
    0.06, #wiltPoint
    0.26, #fieldCap
    0.3,#satPoint
    0.05, #K_s
    3, #V_nr
    0.3, #sigma_zR
    0.01, #E_S1
    0.01, #E_S2
    2, #shared_area
    .4, #maxRootDepth
    0.05, #K_drain
    0.588503613257886, #pFS2
    0.752929538228874, #pFS20
    0.956131627577964, #aS
    0.050456035523466, #nS
    0.384021499609213, #pRx
    0.250229439327847, #pRn
    0.87408236899746, #gammaFx
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
    1,
    10
    
  )
  
  f.increase <-
    c(
      0.25,#wiltPoint
      0.35,#fieldCap
      0.6,#satPoint
      1,#K_s
      5,#V_nr
      3,#sigma_zR
      50,#E_S1
      50,#E_S2
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
      0.90854084808886,#er
      0.6,
      6,
      -5,
      0.9,
      0.2,
      15,
      500
      
    )
  
  ##Need to check what priors we are using!
  pMaxima <- as.vector(unlist(sitka[nm])*(1+(f.increase)))
  pMinima <- as.vector(unlist(sitka[nm])*(1-(f.decrease)))
  pValues <- as.vector(unlist(sitka[nm]))
  
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
  
  sdVals<-(pMaxima-pMinima)
  #ES1 and 2 on very different scales so set manually
  sdVals[7]<-15
  sdVals[8]<-15
  
  
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
  
  pMaxima[[33]]<-0.3
  pMinima[[33]]<-0.01
  
  priorVals <- createTruncatedNormalPrior(mean = pValues, sd=sdVals,
                                          lower = pMinima, upper = pMaxima)
  
  
  return(priorVals)
  
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

