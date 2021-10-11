

library(sensitivity)
library(ggplot2)
library(ggpubr)

###Morris testing
#Simple non-linear
testNonLin<-function(x){
  x[,1]^3
}

#Simple linear
testLin<-function(x){
  x[,1]*3
}

set.seed(10)
#run morris sensitivity analysis for linear function
moOutLin <- morris(
  model = testLin,
  factors = 1, 
  r = 10, 
  design = list(type = "oat", levels = 50, grid.jump=2), 
  binf = 1, 
  bsup = 20, 
  scale = TRUE)

#run morris sensitivity analysis for non-linear function
moOutNonLin <- morris(
  model = testNonLin,
  factors = 1, 
  r = 10, 
  design = list(type = "oat", levels = 50, grid.jump=2), 
  binf = 1, 
  bsup = 20, 
  scale = TRUE)


#tabulate and plot results
senseResLin<-data.frame(func_output=moOutLin$y,inputVal=moOutLin$X[,1])
senseResNonLin<-data.frame(func_output=moOutNonLin$y,inputVal=moOutNonLin$X[,1])

linPLot<-ggplot(senseResLin,aes(x=inputVal,y=func_output))+
  geom_smooth()+
  geom_point(col="red")+
  theme_bw()

nonLinPlot<-ggplot(senseResNonLin,aes(x=inputVal,y=func_output))+
  geom_smooth()+
  geom_point(col="red")+
  theme_bw()

ggarrange(linPLot,nonLinPlot)



###Run with actual model 

#get parameters (currently using monthly timestep for sense analysis)
sitka<-getParms(timeStp = 12)

#read in and update sitka params with current MCMC results
#out<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\monthly_outx_2021-02-265949.RDS")
codM<-as.data.frame(mergeChains(out$chain))

nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","E_S1","E_S2","shared_area","maxRootDepth","K_drain",
      "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
      "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
      "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0","SWpower0","Qa","Qb","MaxIntcptn")



names(codM)<-nm
codM<-miscTools::colMedians(as.data.frame(codM))
sitka[nm]<-codM[nm]

## Likelihood function
likelihoodFunc <- function(p){
  print(p)
sitka$weather[sitka$weather$Year>=2015,6]<-sitka$weather[sitka$weather$Year>=2015,6]*p
    NlogLik <- tryCatch(
    {
      output<-do.call(fr3PGDN,sitka)%>%
      filter(Year>=2015)%>%
        mutate(meanNPP = mean(LAI))
               
      mean(output$meanNPP,na.rm=T)
    },
    error=function(cond) {
      return(NA)
    })
  
  return(NlogLik)
}


##Set priors
priors<-createPriors_sitka(sitka=sitka)

#param names to test for sensitivity
nm<-c("tempMod")

#create priors/ranges for hydrological params, plus add some for rain and temp mods
pMaxima<-priors$upper[1:length(nm)]
pMinima<-priors$lower[1:length(nm)]
pMinima[c(1)]<-.2#up and down by 50% for rainfall
pMaxima[c(1)]<-1.8
pMinima[c(1)]<-0.8#up and down by 10% for temperature
pMaxima[c(1)]<-1.2
Uprior <- createPrior(lower = pMinima, upper = pMaxima)

#observed years to fit to
startYear = 2015
endYear = 2018

##setup likelihood and model running functions
morris_setup <- createBayesianSetup(
  likelihood = likelihoodFunc, 
  prior = Uprior, 
  names = nm)



#set.seed(50)
#run morris sensitivity analysis
morrisOut <- morris(
  model = morris_setup$likelihood$density,
  factors = nm, 
  r = 25, 
  design = list(type = "oat", levels = 50, grid.jump = 3), 
  binf = pMinima, 
  bsup = pMaxima, 
  scale = TRUE)


#tabulate and plot results
senseOutTempX<-data.frame(func_output=morrisOut$y,inputVal=morrisOut$X[,1])
senseOutRain<-data.frame(func_output=morrisOut$y,inputVal=morrisOut$X[,1])

senseOutRain2$Model<-"3PG"
senseOutTemp$Model<-"SonWal"
senseOutRainX<-rbind(senseOutRain,senseOutRain2)

g1<-ggplot(senseOutRain,aes(x=(inputVal*100)-100,y=func_output,col=Model))+
  geom_point(alpha=0.8)+
  geom_smooth()+
  ylab(expression(paste("NPP [gC"," ",cm^-2,"]",sep="")))+
  xlab("Rainfall change (%)")+
  theme_bw()+
  scale_color_viridis_d()

g2<-ggplot(senseOutTemp,aes(x=(inputVal*100)-100,y=func_output,col=Model))+
  geom_point(alpha=0.8)+ 
  geom_smooth()+
  ylab(expression(paste("NPP [gC"," ",cm^-2,"]",sep="")))+
  xlab("Temperature change (%)")+
  theme_bw()+
  scale_color_viridis_d()


ggarrange(g1,g2,common.legend = T, legend ="bottom")


resi<-NULL
parSamp<-getSample(out,start=10000,numSamples=200)
for(i in c(1:nrow(parSamp))){
  nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","E_S1","E_S2","shared_area","maxRootDepth","K_drain",
        "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
        "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
        "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0","SWpower0","Qa","Qb","MaxIntcptn","k","startN","startC")
  

  codM<-parSamp[i,]
  names(codM)<-nm
  
  sitka<-getParms(waterBalanceSubMods=T, timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365)
  sitka[nm]<-codM[nm]
  
  output<-do.call(fr3PGDN,sitka)

  ff<-data.frame(mean=mean(output$LAI),max=max(output$LAI))
  resi<-rbind(resi,ff)
}


for(i in c(1:ncol(parSamp))){
  
  plot(resi$max~parSamp[,i],main=nm[i])
  
}
