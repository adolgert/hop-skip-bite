source("try.R")
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
house.cnt <- as.integer(args[1])
block.cnt <- as.integer(args[2])
save.name <- args[3]
landscape.generate(house.cnt, block.cnt, streetsbase=save.name)


