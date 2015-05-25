library(spatstat) # for ppp data structure
library(scales) # for alpha
source("file_format.R")

plot.all.bug.trajectory <- function(datafile, plotfile) {
  print(paste("data file is", datafile))
  print(paste("plot file is", plotfile))
  pdf(plotfile)
  plot(seq(0, 2, 2/12), seq(0, 1200, 100), type="n")
  get.trajectory <- function(locations, events) {
    with.counts <- subset(events, event!=1)
    lines(with.counts$when, with.counts$who,
      col=alpha("green", 0.1))
  }
  print("foreach.trajectory")
  foreach.trajectory(datafile, get.trajectory)
  dev.off()
}

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
datafile <- args[1]
plotfile <- args[2]
plot.all.bug.trajectory(datafile, plotfile)
