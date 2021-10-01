
##arrays ....AltCor(12),DayL(12),PotHorRad(12),Tb(12),IRad(12),MTemp(12),MHum(12),Mprec(12)
#
#AltCor=array(0,dim=12)
#DayL=array(0,dim=12)
#PotHorRad=array(0,dim=12)
#IRad=array(0,dim=12)
#MTemp=array(0,dim=12)
#MHum=array(0,dim=12)
#Mprec=array(0,dim=12)
#
#Sc = 2.0
#
#frac<-function(x){
#x - round(x,0)
#}
#
#
#
#
##L0 latitude
#HorRad<-function(J,L0){
#  sinD = 0.39785 * sin(4.868961 + 0.017203 * J + 0.033446 * sin(6.224111 + 0.017202 * J))
#  cosD = sqrt(1-(sinD^2))
#  COShs = -sin(L0)*sinD/(cosD*cos(L0))
#  HorRad = 458.37*Sc*(1 + 0.033 * cos(2*pi*0.002739726* J)) * (cos(L0)* cosD * sqrt(1-(COShs^2)) &
#                                                                 + acos(COShs)*180/(pi*57.269)*sin(L0)*sinD)
#  return(data.frame(HorRad=HorRad,sinD=sinD,cosD=cosD,COShs=COShs))
#}
#
######
#NumOfDays<-function(mth){
#  
#  if(mth==1) NumOfDays = 31
#  if(mth==2) NumOfDays = 28
#  if(mth==3) NumOfDays = 31
#  if(mth==4) NumOfDays = 30
#  if(mth==5) NumOfDays = 31
#  if(mth==6) NumOfDays = 30
#  if(mth==7) NumOfDays = 31
#  if(mth==8) NumOfDays = 31
#  if(mth==9) NumOfDays = 30
#  if(mth==10) NumOfDays = 31
#  if(mth==11) NumOfDays = 30
#  if(mth==12) NumOfDays = 31
#  
#  return(NumOfDays)
#  
#}
#
#
#
#DayLengthOf<-function(dy,L0){
#  Ampl = exp(7.42+0.045*L0*180/pi)/3600
#  
#  dl = Ampl*sin((dy-79)*0.1721)+12
#  if(dl<0) dl = 0
#  if(dl>24)dl = 24
#  
#  return(dl)
#}
#
###############
#
#convertTopDat<-function(L0,Is,Aw){
#if (L0>=0) Sign = 1 else Sign = -1.0
#
#L0 = (round(L0,0)+Sign*frac(L0)*0.1666667)*pi/180
#Is = (round(Is,0)+frac(Is)*0.1666667)*pi/180
#Aw = round(Aw,0)+frac(Aw)*0.1666667
#Aw = (180-Aw)*pi/180
#if(Is>(0.5*pi)||Is<0) print('Slope inclination angle out of range')
#
#return(data.frame(L0=L0,Is=Is,Aw=Aw))
#}
#
#
##############
#
#
##k-24
#
#SolarAI<-function(L0,SINa,Is,SINUSa,COSi){
#for (k in c(1:24)){
#  
#  h = (k-12)*15*pi/180
#SNa = sin(L0)*sinD + cos(L0)*cosD*cos(h)
#
#if(SNa>0) {
#  SINa = SINa + SNa
#  CSa = sqrt(1 - (SNa^2))
#  if (Aw == 0)
#    CSi = sin(L0 - Is) * sinD + cos(L0 - Is) * cosD * cos(h)
#  else{
#    As = asin(cosD * sin(-h) / CSa)
#    if (As > pi) {
#      As = pi - As
#    }
#    CSi = cos(As - Aw) * CSa * sin(Is) + SNa * cos(Is)
#  }
#  if (CSi > 0) {
#    COSi = COSi + CSi
#    SINUSa = SINUSa + SNa
#  }
#} else if (k >= 12) break
#
#}
#return(data.frame(SNa=SNa,CSa=CSa,As=As,CSi=CSi,COSi=COSi,SINUSa=SINUSa))
#
#}
#
##InitSolRad
#lat=53
#Incl=1
#ASP=2
#topDat<-convertTopDat(L0=lat,Is=Incl,Aw=ASP)
#L0 = topDat$L0
#Is = topDat$Is
#Aw = topDat$Aw
#
#InitSolRad<-function(L0,Is,Aw){
#  
#  Ampl = exp(7.42+0.045*L0*180/pi)/3600
#  Rcor = 1 - 1.3614*cos(abs(L0))
#  Arad = 32.9835 - 64.884*Rcor
#  Brad = 0.715 - 0.31831*Rcor
#  Td = (cos(Is*0.5))^2
#  
#  
#  #DayL = 0
#  #PotHorRad = 0
#  #Day = 0 
#  
#  for (month = c(1:12)){
#  SINa = 0
#  SINUSa = 0
#  COSi = 0
#  if(Is>0) Tb(month) = 0 else Tb(month) = 1
#  LastDay = NumOfDays(month)
#  for(numDay=c(1:LastDay)){
#    
#    day = day + 1
#    DayL(month) = DayL(month) + DayLengthOf(day)
#    PotHorRad(month) = PotHorRad(month) + HorRad(day)    
#    solAiOut<-SolarAI(L0,SINa,Is,SINUSa,COSi)
#    SNa=solAiOut$SNa
#    CSa=solAiOut$CSa
#    As=solAiOut$As
#    CSi=solAiOut$CSi
#    COSi=solAiOut$COSi
#    SINUSa=solAiOut$SINUSa
#  }
#  
#  SINa = SINa/DayL(month)
#  if(SINa>0) AltCor(month) = (1-exp(-0.014*(Alt-274)/(SINa*274))) else AltCor(month) = 0
#  DayL(month) = DayL(month)/LastDay
#  PotHorRad(month) = PotHorRad(month)/LastDay
#  
#  if(Is>0) { 
#  if(COSi>0) Tb(month) = COSi/SINUSa
#  }
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  }
#  
#}#