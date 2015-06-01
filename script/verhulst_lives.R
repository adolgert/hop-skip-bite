library(spatstat)
library(hopskip)
source("file_format.R")



plot.verhulst.lives <- function(datafile, plotfile) {
  trajectory.info <- single.trajectory(datafile, "dset1")
  locations <- trajectory.info$locations
  df <- trajectory.info$events

  individual.cnt <- length(unique(df$whom))
  plot.time <- max(df$when)

  # Order by Hilbert curve.
  hn <- 2^10
  hilbert.order <- hilbertXY2D(locations$x*(hn^2), locations$y* (hn^2), hn)
  locations.df <- data.frame(x=locations$x, y=locations$y,
    n=seq(0, length(locations$x)-1), hilbert=hilbert.order)
  locations.ordered<-locations.df[order(locations.df$hilbert, locations.df$n),]

  # We want birth/death events for the house, so move/death.
  move <- subset(df, event==1)
  death <- subset(df, event==6 & who==0)

  # Which was the first infected? The first to infect someone else.
  initial.infected.idx <- move$who[1]
  print(locations.ordered)
  print(move)
  print(death)
  print(initial.infected.idx)
}