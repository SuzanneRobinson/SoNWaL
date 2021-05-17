fr3PGDN <-
function (weather, presc, t = 0, N = 2500, Wl = 0.01, WlDormant = 0, Wr = 0.01, Wsbr = 0.1, 
              Wlitt = 0, YrC= 18.2, YlC= 6.5, OC= 122.77, YrN= 0, YlN= 1.52, ON= 6.9,
              Nav = 3, rotation = 1, cycle = 1, rm.sprouts = F, nyears = 35,
              initial.month = 1, latitude = 57.06, soilclass = -1, ASW = 200,
              MaxASW = 300,  MinASW = 0, CO2 = 400, pFS2 = 0.76, pFS20 = 0.46, aS = 0.0572,
              nS = 2.4568, pRx = 0.48, pRn = 0.21, Tmin = -5, Topt = 15, Tmax = 35, kF = 1,
              SWconst0 = 0.7,  SWpower0 = 9, m0 = 0.018, MaxAge = 265.6,
              nAge = 3.545,  rAge = 0.796, gammaFx = 0.025, gammaF0 = 0.0022, tgammaF = 60,
              Rttover = 0.07, MaxCond = 0.02, LAIgcx = 3.33, BLcond = 0.2,  wSx1000 = 255,
              thinPower = 1.5, mF = 0.2, mR = 0.2, mS = 0.2,  SLA0 = 6, SLA1 = 4, tSLA = 6,
              k = 0.52, fullCanAge = 18,  MaxIntcptn = 0.15, LAImaxIntcptn = 5, alpha = 0.05,
              Y = 0.49,  poolFractn = 0, e20 = 2.2, rhoAir = 1.2, lambda = 2460000,  VPDconv = 0.000622,
              fracBB0 = 0.3, fracBB1 = 0.1, tBB = 10,  rhoMin = 0.39, rhoMax = 0.39, tRho = 5,
              Qa = -90, Qb = 0.8,  gDM_mol = 24, molPAR_MJ = 2.3, CoeffCond = 0.05, fCalpha700 = 1.433,
              fCg700 = 0.451, fCalphax = 2.33333333333333, fCg0 = 1.75, MinCond = 0.015,
              klmax = 0.031532358, krmax = 0.0042394358, komax = 0.000458861, hc = 0.23, qir = 334.29,
              qil = 49.08, qh = 23.63, qbc = 2.21, el = 0.25, er = 0.56, Nf = 0.0684, Navm = 0.01, Navx = 20,
              leaf.grow = 0, leaf.fall = 0, Wl.s = 0.01, Wsbr.s = 0.1, Wr.s = 0.01, pWl.sprouts = 0.5,
              pWsbr.sprouts = 0.9, cod.pred = "3PG", cod.clim = "Average" ,K_s=0.1,V_nr=5,sigma_zR =1,SWC_nr=800,E_S=100,
              E_S1 =50,E_S2 =0.01,waterBalanceSubMods=T,wiltPoint=0.1,fieldCap=0.3,MaxASW_state=50,timeStp=12,maxRootDepth=2,shared_area=5,
          K_drain=0.1,SWC_rz,satPoint=0.35,tempMod=1,rainMod=1)
        ## added WlDormant, leafgrow and leaffall, leafgrow and leaffall - Tom Locateli
        ## pfsPower = -0.522878745280338, pfsConst = 0.43104582317421,
{
    
  
  maxRootDepth<-min(maxRootDepth,V_nr) 
  
    parms <- c(pFS2, pFS20, aS, nS, pRx, pRn, Tmin, Topt, Tmax, kF, SWconst0, SWpower0, m0,
               MaxAge, nAge, rAge, gammaFx, gammaF0, tgammaF, Rttover, MaxCond, LAIgcx, BLcond,
               wSx1000, thinPower, mF, mR, mS, SLA0, SLA1, tSLA, k, fullCanAge, MaxIntcptn, LAImaxIntcptn,
               alpha, Y, poolFractn, e20, rhoAir, lambda, VPDconv, fracBB0, fracBB1, tBB, rhoMin, rhoMax,
               tRho, Qa, Qb, gDM_mol, molPAR_MJ, CoeffCond, fCalpha700, fCg700, fCalphax, fCg0,
               klmax, krmax, komax, hc, qir, qil, qh, qbc, el, er, Nf, Navm, Navx,
               MinCond, Wl.s, Wsbr.s, Wr.s, pWl.sprouts, pWsbr.sprouts, leaf.grow, leaf.fall,K_s,V_nr,sigma_zR, 
               E_S1,E_S2,waterBalanceSubMods,wiltPoint,fieldCap,timeStp,maxRootDepth,shared_area,K_drain,satPoint)

    names(parms)<-c("pFS2","pFS20","aS","nS","pRx","pRn",
                    "Tmin","Topt","Tmax","kF",
                    "SWconst0","SWpower0","m0", "MaxAge",
                    "nAge","rAge","gammaFx","gammaF0","tgammaF","Rttover",
                    "MaxCond","LAIgcx","BLcond","wSx1000","thinPower",
                    "mF","mR","mS","SLA0","SLA1","tSLA","k","fullCanAge",
                    "MaxIntcptn","LAImaxIntcptn","alpha","Y","poolFractn",
                    "e20","rhoAir","lambda","VPDconv","fracBB0","fracBB1",
                    "tBB","rhoMin","rhoMax","tRho","Qa","Qb","gDM_mol",
                    "molPAR_MJ","CoeffCond","fCalpha700","fCg700","fCalphax","fCg0",
                    "klmax","krmax","komax","hc","qir","qil","qh","qbc",
                    "el","er","Nf","Navm","Navx",
                    "MinCond","Wl.s","Wsbr.s","Wr.s","pWl.sprouts",
                    "pWsbr.sprouts","leaf.grow","leaf.fall","K_s","V_nr","sigma_zR",
                    "E_S1","E_S2","waterBalanceSubMods","wiltPoint","fieldCap","timeStp","maxRootDepth","shared_area","K_drain","satPoint")

    vars.ini <- c(t, N, Wl, WlDormant, Wr, Wsbr, Wlitt, YrC, YlC, OC, YrN, YlN, ON, Nf, Nav,
                  rotation, cycle, rm.sprouts, nyears, initial.month)
    
    names(vars.ini) <- c("t", "N", "Wl", "WlDormant", "Wr", "Wsbr", "Wlitt", "YrC", "YlC",
                         "OC", "YrN", "YlN", "ON", "Nf", "Nav", "rotation", "cycle",
                         "rm.sprouts", "nyears", "initial.month")

    site.info <- c(latitude, soilclass, ASW, MaxASW, MinASW, 
                   CO2,SWC_nr,E_S,MaxASW_state,E_S)

    names(site.info) <- c("latitude", "soilclass", "ASW", 
                          "MaxASW", "MinASW", "CO2","SWC_nr","E_S","MaxASW_state","SWC_rz")

    parms.general <- list(daysinmonth = c(Jan = 31, Feb = 28, 
                                          Mar = 31, Apr = 30,
                                          May = 31, Jun = 30,
                                          Jul = 31, Aug = 31, 
                                          Sep = 30, Oct = 31,
                                          Nov = 30, Dec = 31),
                          parms.soil = data.frame(soilclass.name = c("Sandy", "Sandy loam",
                                                                     "Clay loam", "Clay",
                                                                     "Non standard", "No effect of ASW"),
                                                  soilclass = c(1, 2, 3, 4, -1, 0),
                                                  SWconst = c(0.7, 0.6, 0.5, 0.4, parms[["SWconst0"]], 1),
                                                  SWpower = c(9, 7, 5, 3, parms[["SWpower0"]], 1)))

    ##For spatial increase sensitivity analysis, temp and rainfall mods
    weather$Rain<-weather$Rain*rainMod
    weather[,c("Tmean","Tmax","Tmin")]<-weather[,c("Tmean","Tmax","Tmin")]*tempMod
    
    
    proj <- RunModel(stand.init = vars.ini, weather = weather, 
                     site = site.info, general.info = parms.general, presc = presc, 
                     parms = parms, cod.pred = cod.pred, cod.clim = cod.clim)
    
    return(as.data.frame(proj))
}
