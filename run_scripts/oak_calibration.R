# oak calibration
library(SoNWaL)
library(BayesianTools)
library(tidyverse)
library(dplyr)
library(coda)
library(miscTools)

#read in arguments from batch file
#args=(commandArgs(TRUE))
args<-c("weekly","weekly_1_","1","T")

print(args)
timeStep=args[1]
chainID=args[2]
chainNum=args[3]
#print arguments to log for easy reference
print(timeStep)
print(chainID)
print("oak")

#Time step to run SonWal with
timeStep<-"weekly"
args[1]<-timeStep

fileLocs<-list.files("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\NZplus\\AH_flux\\Flux data_GapFill\\Flux data_GapFill",
                     pattern = "DataSet",full.names = T)
ff<-read.csv(fileLocs[1])

