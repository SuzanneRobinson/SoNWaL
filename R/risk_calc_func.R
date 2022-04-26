#clm_df_full<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\NZplus\\misc_data\\AH_Harwood_CHESS_dat\\clm_hist.RDS")%>%
#  filter(site=="Harwood" & clm_var=="precip")%>%
#  mutate(value=value*86400, month=lubridate::month(date))%>%
#  select(date,month,year,value)
#
#
#monthlyPrecip <-
#  clm_df_full %>% group_by(year,month) %>% summarise(monthly_precip = sum(value))
#
##get daytime values
#
##split into seasons
#monthlyPrecip$season<-ifelse(monthlyPrecip$month==3 | monthlyPrecip$month ==4| monthlyPrecip$month==5, "spring" , "winter")
#monthlyPrecip$season<-ifelse(monthlyPrecip$month==6 | monthlyPrecip$month ==7| monthlyPrecip$month==8, "summer" , monthlyPrecip$season)
#monthlyPrecip$season<-ifelse(monthlyPrecip$month==9 | monthlyPrecip$month ==10| monthlyPrecip$month==11, "autumn" , monthlyPrecip$season)
#
#
##calc quantiles for each season
#
#hazprecip <- monthlyPrecip%>%group_by(season)%>%summarise(seasQ=quantile(monthly_precip, 0.05))%>%
#  right_join(monthlyPrecip,by="season")%>%
# arrange(year,month)%>%
#  mutate(hazMonth=ifelse(monthly_precip<=seasQ,1,0))%>%
#  group_by(year)%>%
#  summarise(hazMonths=sum(hazMonth), rain=mean(monthly_precip))%>%
#  mutate(haz=ifelse(hazMonths>=2,1,0))
#
#ggplot(data=hazprecip, aes(x=year,y=rain))+
#  geom_point(aes(colour=as.factor(haz)))+
#  geom_line(col="blue")
#
#
#
#
##Risk function
#risk_calc_func <- function(strtyr, endyr, df, hazval, modOut) {
#  fldf    <- df %>% filter(year >= strtyr & year <= endyr)
#  lowyrs  <- fldf$year[(fldf$hazMonths >= 2)]
#  highyrs <- fldf$year[!(fldf$year %in% lowyrs)]
#  yearsAv<-clm_df_full%>%filter(year>=min(modOut$Year,na.rm=T) & year<=max(modOut$Year,na.rm=T)) %>%
#    select(year)
#  
#  vulnHigh<-modOut%>%filter(Year%in% highyrs)%>%group_by(Year)%>%summarise(GPP=mean(GPP))
#  vulnLow<-modOut%>%filter(Year%in% lowyrs)%>%group_by(Year)%>%summarise(GPP=mean(GPP))
#  vuln    <-
#    mean(vulnHigh$GPP) - mean(vulnLow$GPP)
#  haz     <- length(lowyrs) / ((endyr - strtyr) + 1)
#  return(
#    tibble(
#      "startYr" = as.character(strtyr),
#      "endYr" = as.character(endyr),
#      "vulnerability" = vuln,
#      "hazard"  = haz,
#      "risk" = vuln * haz
#    )
#  )
#}#