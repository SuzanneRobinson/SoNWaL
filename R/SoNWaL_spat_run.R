

# Spatial run function which takes a single grid cell and associated climate variables, 
# runs simulations for a range of MCMC posterior values
# Calculates the mean and a 95% credible interval
# Currently there is some data manipulation going on to spread a single years data
# over multiple years until all data is downloaded
#' @param site site data
#' @param clm climate data
#' @param param_draw parameter draws
#' @export
SoNWaL_spat_run <-
  function(site,
           clm,
           param_draw,
           grid_id,
           soil_depth,
           wp,
           fc,
           sp,
           carbon,
           cond,
           N0,
           C0,
           plant_year,scape=T, 
           CN_ratio=10) {
    
   clmY<- clm%>%
      tibble::rownames_to_column()
    
    clmY<-clmY[,1]
    clm$pyear<-stringr::str_sub(clmY,start=2, end =5)
    clm<-filter(clm,pyear>=plant_year)
    

    library(lubridate)

    # update parameters with fixed values from soil data etc.
      param_draw$pars <- lapply(param_draw$pars, function(x) {
        x$V_nr <- soil_depth/100
        x$maxRootDepth <- soil_depth/100
        x$startC <- carbon
        x$startN <- carbon/CN_ratio
        x$wiltPoint <- wp
        x$fieldCap <- fc
        x$satPoint <- sp
        # currently K_s and drainage use same value from maps, could be improved?
        x$K_s <- cond
        x$K_drain <- cond
        x$K_drain_nrz <- cond
        #######################
        x$SWpower0<-N0
        x$SWconst0<-C0
        return(x)
      })
      
      
    # return empty cell output if there is no climate data (no land mass)
    if (is.na(clm[1, 2]) == T) {
      print("no climate data")
      return (tibble::as_tibble(
        data.frame(
          Year = NA,
          grid_id = grid_id,
          Wsbr_q05 = NA,
          Wsbr_q95 = NA,
          Wsbr_value = NA,
          Wsbr_var = NA,
          
          Rs_q05 =
            NA,
          Rs_q95 = NA,
          Rs_value = NA,
          Rs_var = NA,
          
          EvapTransp_q05 =
            NA,
          EvapTransp_q95 = NA,
          EvapTransp_value = NA,
          EvapTransp_var = NA,
          
          volSWC_rz_q05 =
            NA,
          volSWC_rz_q95 = NA,
          volSWC_rz_value = NA,
          volSWC_rz_var = NA,
          
          yc_q05 =
            NA,
          yc_q95 = NA,
          yc_value = NA,
          yc_var = NA,
          
          GPP_q05 =
            NA,
          GPP_q95 = NA,
          GPP_value = NA,
          GPP_var = NA,
          
          NPP_q05 =
            NA,
          NPP_q95 = NA,
          NPP_value = NA,
          NPP_var = NA,
          
          NEE_q05 =
            NA,
          NEE_q95 = NA,
          NEE_value = NA,
          NEE_var = NA,
          
          Reco_q05 =
            NA,
          Reco_q95 = NA,
          Reco_value = NA,
          Reco_var = NA,
          
          LAI_q05 =
            NA,
          LAI_q95 = NA,
          LAI_value = NA,
          LAI_var = NA
        )
      )) 
      } else {


      site_out <-
        tryCatch({

          param_draw %>%
            dplyr::mutate(sim = mapply(SoNWaL_spat_model_run, pars, MoreArgs = list(clm,scape=scape),SIMPLIFY = F)) %>%
            dplyr::select(mcmc_id, sim) %>%
            tidyr::unnest_legacy() %>%
            group_by(Year, mcmc_id) %>%
            dplyr::summarise(
              grid_id = grid_id,
              dg = mean(dg,  na.rm = T),
              Rs = mean(Rs, na.rm = T),
              EvapTransp = mean(EvapTransp, na.rm = T),
              volSWC_rz = mean(volSWC_rz, na.rm = T),
              yc = mean(yc, na.rm = T),
              GPPsum = sum(GPP,  na.rm = T),
              NPPsum = sum(NPP, na.rm = T),
              NEEsum = sum(NEE, na.rm = T),
              GPP = mean(GPP,  na.rm = T),
              NPP = mean(NPP, na.rm = T),
              NEE = mean(NEE, na.rm = T),
              Reco = mean(Reco, na.rm = T),
              LAI = mean (LAI, na.rm = T)
              
            ) %>%
            group_by(Year) %>%
            dplyr::summarise(
              grid_id = grid_id,
              Wsbr_q05 = quantile(dg, 0.05, na.rm = T),
              Wsbr_q95 = quantile(dg, 0.95, na.rm = T),
              Wsbr_value = quantile(dg, 0.5, na.rm = T),
              Wsbr_var = var(dg, na.rm = T),
              
              Rs_q05 = quantile(Rs, 0.05, na.rm = T),
              Rs_q95 = quantile(Rs, 0.95, na.rm = T),
              Rs_value = quantile(Rs, 0.5, na.rm = T),
              Rs_var = var(Rs, na.rm = T),
              
              EvapTransp_q05 = quantile(EvapTransp, 0.05, na.rm = T),
              EvapTransp_q95 = quantile(EvapTransp, 0.95, na.rm = T),
              EvapTransp_value = quantile(EvapTransp, 0.5, na.rm = T),
              EvapTransp_var = var(EvapTransp, na.rm = T),
              
              volSWC_rz_q05 = quantile(volSWC_rz, 0.05, na.rm = T),
              volSWC_rz_q95 = quantile(volSWC_rz, 0.95, na.rm = T),
              volSWC_rz_value = quantile(volSWC_rz, 0.5, na.rm = T),
              volSWC_rz_var = var(volSWC_rz, na.rm = T),
              
              yc_q05 = quantile(yc, 0.05, na.rm = T),
              yc_q95 = quantile(yc, 0.95, na.rm = T),
              yc_value = quantile(yc, 0.5, na.rm = T),
              yc_var = quantile(yc, na.rm = T),
              
              GPP_q05 = quantile(GPP, 0.05, na.rm = T),
              GPP_q95 = quantile(GPP, 0.95, na.rm = T),
              GPP_value = quantile(GPP, 0.5, na.rm = T),
              GPP_var = quantile(GPP, na.rm = T),
              
              NPP_q05 = quantile(NPP, 0.05, na.rm = T),
              NPP_q95 = quantile(NPP, 0.95, na.rm = T),
              NPP_value = quantile(NPP, 0.5, na.rm = T),
              NPP_var = quantile(NPP, na.rm = T),
              
              NEE_q05 = quantile(NEE, 0.05, na.rm = T),
              NEE_q95 = quantile(NEE, 0.95, na.rm = T),
              NEE_value = quantile(NEE, 0.5, na.rm = T),
              NEE_var = quantile(NEE, na.rm = T),
              
              GPPsum_q05 = quantile(GPPsum, 0.05, na.rm = T),
              GPPsum_q95 = quantile(GPPsum, 0.95, na.rm = T),
              GPPsum_value = quantile(GPPsum, 0.5, na.rm = T),
              GPPsum_var = quantile(GPPsum, na.rm = T),
              
              NPPsum_q05 = quantile(NPPsum, 0.05, na.rm = T),
              NPPsum_q95 = quantile(NPPsum, 0.95, na.rm = T),
              NPPsum_value = quantile(NPPsum, 0.5, na.rm = T),
              NPPsum_var = quantile(NPPsum, na.rm = T),
              
              NEEsum_q05 = quantile(NEEsum, 0.05, na.rm = T),
              NEEsum_q95 = quantile(NEEsum, 0.95, na.rm = T),
              NEEsum_value = quantile(NEEsum, 0.5, na.rm = T),
              NEEsum_var = quantile(NEEsum, na.rm = T),
              
              Reco_q05 = quantile(Reco, 0.05, na.rm = T),
              Reco_q95 = quantile(Reco, 0.95, na.rm = T),
              Reco_value = quantile(Reco, 0.5, na.rm = T),
              Reco_var = quantile(Reco, na.rm = T),
              
              LAI_q05 = quantile(LAI, 0.05, na.rm = T),
              LAI_q95 = quantile(LAI, 0.95, na.rm = T),
              LAI_value = quantile(LAI, 0.5, na.rm = T),
              LAI_var = quantile(LAI, na.rm = T)
              
            )
        },
        #add na values where there is no data, in the sea etc.
        error = function(cond) {
          message(cond)
          
          site_out <-
            tibble::as_tibble(
              data.frame(
                Year = NA,
                grid_id = grid_id,
                Wsbr_q05 = NA,
                Wsbr_q95 = NA,
                Wsbr_value = NA,
                Wsbr_var = NA,
                
                Rs_q05 =
                  NA,
                Rs_q95 = NA,
                Rs_value = NA,
                Rs_var = NA,
                
                EvapTransp_q05 =
                  NA,
                EvapTransp_q95 = NA,
                EvapTransp_value = NA,
                EvapTransp_var = NA,
                
                volSWC_rz_q05 =
                  NA,
                volSWC_rz_q95 = NA,
                volSWC_rz_value = NA,
                volSWC_rz_var = NA,
                
                yc_q05 =
                  NA,
                yc_q95 = NA,
                yc_value = NA,
                yc_var = NA,
                
                GPP_q05 =
                  NA,
                GPP_q95 = NA,
                GPP_value = NA,
                GPP_var = NA,
                
                NPP_q05 =
                  NA,
                NPP_q95 = NA,
                NPP_value = NA,
                NPP_var = NA,
                
                NEE_q05 =
                  NA,
                NEE_q95 = NA,
                NEE_value = NA,
                NEE_var = NA,
                
                Reco_q05 =
                  NA,
                Reco_q95 = NA,
                Reco_value = NA,
                Reco_var = NA,
                
                LAI_q05 =
                  NA,
                LAI_q95 = NA,
                LAI_value = NA,
                LAI_var = NA
                
              )
            )
        })
    
    
    site_out <-
      site_out[!duplicated(site_out), ]
    return(as.data.frame(site_out))
    
  }

}



