source("try.R")
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
housename <- args[1]
streetsbase <- args[2]
landscape.crossings(housesfile=housename, streetsfile=streetsbase)

