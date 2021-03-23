

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
out<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\monthly_outx_2021-02-265949.RDS")
codM<-as.data.frame(mergeChains(out$chain))
nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","E_S1","E_S2","shared_area","maxRootDepth","K_drain",
      "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
      "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
      "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er")
names(codM)<-nm
codM<-colMedians(as.data.frame(codM))
sitka[nm]<-codM[nm]

## Likelihood function
likelihoodFunc <- function(p){
  sitka[.GlobalEnv$nm]<-p
  NlogLik <- tryCatch(
    {
      output<-do.call(fr3PGDN,sitka)%>%
      filter(Year==2018)%>%
        mutate(meanNPP = mean(NPP))
               
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
nm<-c("rainMod")

#create priors/ranges for hydrological params, plus add some for rain and temp mods
pMaxima<-priors[[1]][1:length(nm)]
pMinima<-priors[[2]][1:length(nm)]
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
  model = morris_setup$posterior$density,
  factors = nm, 
  r = 25, 
  design = list(type = "oat", levels = 50, grid.jump = 3), 
  binf = pMinima, 
  bsup = pMaxima, 
  scale = TRUE)


#tabulate and plot results
senseOutTemp<-data.frame(func_output=morrisOut$y,inputVal=morrisOut$X[,1])
senseOutRain<-data.frame(func_output=morrisOut$y,inputVal=morrisOut$X[,1])

senseOutRain2$inputVal2<-"noHS"
senseOutRain$inputVal2<-"HS"
senseOutRainX<-rbind(senseOutRain,senseOutRain2)

g1<-ggplot(senseOutRainX,aes(x=(inputVal*100)-100,y=func_output,col=inputVal2))+
  geom_point(alpha=1)+
  geom_smooth()+
  ylab(expression(paste("NPP [tDM"," ",ha^-1,"]",sep="")))+
  xlab("Rainfall change (%)")+
  theme_bw()

g2<-ggplot(senseOutTemp,aes(x=(inputVal*100)-100,y=func_output))+
  geom_point(col="red",alpha=0.8)+ 
  geom_smooth()+
  ylab(expression(paste("NPP [tDM"," ",ha^-1,"]",sep="")))+
  xlab("Temperature change (%)")+
  theme_bw()

ggarrange(g1,g2)
