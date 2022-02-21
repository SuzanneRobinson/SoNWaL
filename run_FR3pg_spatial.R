library("tidyr")
library("tidyverse")
library("purrr")
library("BayesianTools")
library("sensitivity")
library("dplyr")
library("future")
library("ncdf4)" )
library("raster)" )
library("ggplot2")
library("httr")
library("furrr")
library("viridis")
library("fr3PGDN")
library("tibble")
library("miscTools")
library("parallel")
library("sf")
library("rgdal")
library(raster)


data_dir="/gws/nopw/j04/hydro_jules/data/uk/driving_data/chess/chess-met/daily"#Location of input data
output_dir="/home/users/aaronm7/"#output file for bricked rasters
save_file="/home/users/aaronm7/spatialChunks/"#save file for spatial chunks of data

output_dir="C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\"
data_dir="C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\Prafor\\models\\spatial_met_data\\CHESSscape\\daily"#Location of input data
save_file="C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\Prafor\\spatialChunks\\"#save file for spatial chunks of data


spat_split_scape(data_dir=data_dir,output_dir=output_dir,save_file=save_file,dates=1960:2100,num_chunks=c(1:31))
#create chunks of approx 10,000 grid cells in parallel, saves to file
#1:35 is approx scotland
plan(multisession, workers = 8)

future_map(c(31), ~spat_dat_scape(chunk=.x,output_dir,save_file),.progress = T)


####
