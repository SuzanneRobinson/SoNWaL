###Morris sensitivity analysis
library(sensitivity)
library(viridis)

###########################
## Initialise Parameters ##
###########################
sitka<-getParms(weather=clm_df_full,
                waterBalanceSubMods =T, #Whether to run model using updated water balance submodels
                wiltPoint = 0.1, #Wilting point in m^3/m^3? need to convert to mm per meter with rooting depth?
                fieldCap =0.29,#Field capacity
                satPoint= 0.45, #field saturation point
                K_s=0.1, #Soil conductivity
                shared_area=4, #shared area of rooting and non-rooting zone
                V_nr=3, #Volume of non-rooting zone
                maxRootDepth=2,
                sigma_zR =0.7, #area/depth explored by 1kg of root biomass
                E_S1 =0.1, #Cumulitive evap threshold (kg^m-2) - sensitive to length of time-step, e.g. monthly time-step means wetting event only occurs at end of month
                E_S2 =0.3, #how quickly evaporation rate declines with accumulated phase 2 evaporation - based on soil structure
                K_drain=0.1,
                timeStp = 12 # time step, 52 for weekly, 12 for monthly and 365 for daily
)


#get priors and values for hydrology sub-model parameters
priors<-createPriors_sitka(sitka=sitka)
pMaxima<-priors[[1]]
pMinima<-priors[[2]]
Uprior <- createUniformPrior(lower = pMinima[1:11], upper = pMaxima[1:11])
nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","E_S1","E_S2","shared_area","maxRootDepth","K_drain")


##setup likelihood and model running functions for sense analysis
morris_setup <- createBayesianSetup(
  likelihood = NLL, 
  prior = Uprior, 
  names = nm)

set.seed(10)

#run morris sensitivity analysis
morrisOut <- morris(
  model = morris_setup$posterior$density,
  factors = nm, 
  r = 100, 
  design = list(type = "oat", levels = 20, grid.jump = 3), 
  binf = pMinima[1:11], 
  bsup = pMaxima[1:11], 
  scale = TRUE)

# summarise the sensitivity analysis
morrisOut.df <- data.frame(
  parameter = nm,
  mu.star = apply(abs(morrisOut$ee), 2, mean, na.rm = T),
  sigma = apply(morrisOut$ee, 2, sd, na.rm = T)
) %>%
  arrange( mu.star )


#plot results
morrisOut.df %>%
  gather(variable, value, -parameter) %>%
  ggplot(aes(reorder(parameter, value), value, fill = variable), color = NA)+
  geom_bar(position = position_dodge(), stat = 'identity') +
  scale_fill_viridis_d("", option = "D", labels = c('mu.star' = expression(mu * "* (high influence) "), 
                                      'sigma' = expression(sigma ~"(high interaction) "))) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5),
    axis.title = element_blank(),
    legend.position = c(0.05 ,0.95),legend.justification = c(0.05,0.95)
  )
