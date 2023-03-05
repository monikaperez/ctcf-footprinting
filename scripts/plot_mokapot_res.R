# plot mokapot output
library(ggplot2)
library(data.table)
library(plyr)
library(dplyr)
library(tictoc)
library(ggsci)
library(purrr)
library(scales)


#===================================================================
#========================INPUT PARAMETERS===========================
#===================================================================

args <- commandArgs()
print(args)

input_file = args[1]
output_dir = args[2]
data_dir = args[3]
file_root = args[4]

#===================================================================
#========================HELPER FUNCTIONS===========================
#===================================================================


#------------ Define dirs ------------
setwd("/mmfs1/gscratch/stergachislab/mwperez/ctcf-footprinting")

#ARGS
data_folder <- sprintf("%s/candidate_footprints", getwd())
mokapot_dir <- sprintf("%s/mokapot_res", getwd())
output_folder <- sprintf("%s/figures", mokapot_dir)
