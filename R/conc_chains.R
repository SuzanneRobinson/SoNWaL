# concatonate site specific params from regional mcmc chains
#' @param mcmcReg regional mcmc output
#' @param plot whether to plot distributions
#' @return concChains concatonate chains
conc_chains<-function(mcmcReg,plot=F,clm_df_reg){

  #concatonate chains
  nm_all<-c(paste0("wiltPoint_Si",unique(clm_df_reg$siteName)),
            paste0("fieldCap_Si",unique(clm_df_reg$siteName)),
            paste0("satPoint_Si",unique(clm_df_reg$siteName)),
            paste0("K_s_Si",unique(clm_df_reg$siteName)),
            paste0("V_nr_Si",unique(clm_df_reg$siteName)),
            paste0("E_S1_Si",unique(clm_df_reg$siteName)),
            paste0("E_S2_Si",unique(clm_df_reg$siteName)),
            paste0("shared_area_Si",unique(clm_df_reg$siteName)),
            paste0("maxRootDepth_Si",unique(clm_df_reg$siteName)),
            paste0("K_drain_Si",unique(clm_df_reg$siteName)),
            paste0("startN_Si",unique(clm_df_reg$siteName)),
            paste0("startC_Si",unique(clm_df_reg$siteName)),
            "pFS2","pFS20","gammaF0","tgammaF","Rttover","mF","mR",
            "mS","Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0",
            "SWpower0","sigma_zR"
            ,"aS","nS","pRx","pRn","gammaFx",
            "SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
            "k","Qa","Qb","MaxIntcptn","llp","ll","pp")
  
  mChains<- if (is.null(nrow(mcmcReg))==T) as.data.frame((mcmcReg$chain[[1]])) else as.data.frame(mcmcReg)
  names(mChains)<-if (is.null(nrow(mcmcReg))==T) nm_all else nm_all[-c(200:202)]
  concChains<-mChains %>% 
    pivot_longer(cols = starts_with(c("wiltPoint","fieldCap","satPoint","K_s","V_nr","E_S1",
                                      "E_S2","shared_area","maxRootDepth",
                                      "K_drain","startN","startC")), 
                 names_to = c(".value", "wpKey"), names_sep = "_Si") %>% 
    dplyr::select(-wpKey)
  
  
  concChains<- if (is.null(nrow(mcmcReg))==F) distinct(concChains, pFS2, pFS20,gammaF0,tgammaF,Rttover, .keep_all=T ) else concChains
  
  if(plot==T){
    gpL<-list()  
    for(i in c(1:ncol(concChains))){
      
      concTmp<-data.frame(cc=concChains[i])
      names(concTmp)<-"cc"
      gpL[[i]]<-ggplot(data=concTmp,aes(cc))+
        geom_histogram(bins=100,col="black")+
        xlab(names(concChains[i]))
      
    }    
    ggarrange(plotlist=gpL)
  }
  return(concChains)
}
