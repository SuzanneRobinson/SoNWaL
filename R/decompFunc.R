#dg<-0
#difRoots<-0
#YlC<-2
#YrC<-2
#OC<-2
#YlN<-2
#YrN<-2
#ON<-2
#qbc<-2.21427684
#el<-0.24636719
#hc<-0.2
#difWl<-0
#Nf<-0.00684
#Nav<-12
#krmax<-0.01
#klmax<-0.00423943
#komax<-0.00045886
#qh<-23.6348669
#er<-0.5612215
#difLitter<-0.2
#qir<-334.290515
#qil<-49.0841127
#Navm<-0.01
#Navx<-10
#t<-1e+05
#soilDecomp<-function(difLitter,dg,difRoots,YlC,YrC,OC,YlN,YrN,ON,qbc,el,er,hc,difWl,Nf,Nav,krmax,klmax,komax,qh,qir,qil){
#  fSW<-1
#  fT<-1
#  Wds<-0
#  Wdl<-0
#  Wdr<-0
#  ONres<-rep(0,length(t))
#  OCres<-rep(0,length(t))
#  for(i in c(1:t)){
#  ## Calculate the fineCoarseRootRatio based on DBH    
#  fineCoarseRatio <- 0.1276+1462.2671*exp(-1.7958*dg) ## Empirical relationship for Scots Pine - GX !!! Requires a revised method !!!
#  coarseDifRoots <- difRoots / (1 + fineCoarseRatio)
#  fineDifRoots <- difRoots - coarseDifRoots
#  ##Calculate the decomposition rate    
#  kr <- mean(output$kr)#krmax * fSW * fT
#  kl <- mean(output$kl)#klmax * fSW * fT
#  ko <- mean(output$ko)#komax * fSW * fT  
#  ##Calculate fluxes in, out and between carbon and nitrogen pools
#  ##Carbon fluxes
#  YlCflx <- kl * (1 - hc) * YlC
#  YrCflx <- kr * (1 - hc) * YrC
#  OCflx <- ko * OC   
#  ##Humification coefficients
#  hl <- kl * hc * YlC
#  hr <- kr * hc * YrC
#  ##Humification coefficients for nitrogen pools
#  hNl <- kl * hc * (YlN / qh)
#  hNr <- kr * hc * (YrN / qh)
#  ##Nitrogen fluxes
#  YlNflx <- kl * ((1 - hc) / (1 - el)) * (YlN - el * (YlC / qbc))
#  YlNflx <- max(0.,YlNflx)
#  YrNflx <- kr * ((1 - hc) / (1 - er)) * (YrN - er * (YrC / qbc))
#  YrNflx <- max(0.,YrNflx)
#  ONflx <- ko * ON
#  ONflx <- max(0.,ONflx)
#  ##Calculate heterotrophic respiration
#  Rs <- YlCflx + YrCflx + OCflx
#  ##Now calculate carbon and nitrogen pools
#  YrC <- YrC + ( (Wds + coarseDifRoots) / 2) - YrCflx - hr
#  YlC <- YlC + ((difLitter + fineDifRoots + Wdl + Wdr) / 2) - YlCflx - hl
#  OC <- OC + hl + hr - OCflx
#  YrN <- YrN + ((Wds + coarseDifRoots) / (2 * qir)) - YrNflx - hNr
#  YlN <- YlN + ((difLitter + fineDifRoots + Wdl + Wdr) / (2 * qil)) - YlNflx - hNl
#  ON <- ON + hNr + hNl - ONflx
#
#
#    ## Calculate the total pools
#  totC <- YrC + YlC + OC
#  totN <- YrN + YlN + ON
#  
#  ONres[i]<-totN
#  OCres[i]<-totC
#  ## Calculate the Fertility Rating
#  ## First calculate available nitrogen
#  Navflx <- YrNflx + YlNflx + ONflx
#  ## Calculate nitrogen uptake
#  Un <- difWl * Nf
#  ## Update available nitrogen pool
#  Nav <- Nav + Navflx - Un
#  Nav <- ifelse(is.na(max( Navm, Nav )),0,max( Navm, Nav))  ## A small change to accommodate Bayesian calibration and avoid NAs
#  if( Nav > Navx){
#    Nleach <- Nav - Navx
#    Nav <- Navx
#  }else{
#    Nleach <- 0
#  }
#  ## Now estimate fN
#  fN <- (Nav - Navm) / (Navx - Navm)
#  
#  }
#}
#
#OC+YlC+YrC
#ON+YlN+YrN
#
#YlN/ON
#YrN/ON
#YrN/totN
#
#OC/totC
#YlC/totC
#YrC/totC
#
#plot(ONres[seq(1, length(ONres), 1000)])
#plot(OCres[seq(1, length(ONres), 1000)])
#