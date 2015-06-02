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
  norder <- numeric(individual.cnt)
  for (locidx in 1:nrow(locations.ordered)) {
    norder[locations.ordered$n[locidx]+1]=locidx
  }
  # Or don't reorder until we figure out how to put house numbers on the plot.
  norder <- seq(1, max(df$whom)+1)

  # We want birth/death events for the house, so move/death.
  move <- subset(df, event==1)
  death <- subset(df, event==6 & who==0)
  birthdeath <- rbind(move, death)
  birthdeath <- birthdeath[order(birthdeath$when),]
  print(birthdeath)

  # Which was the first infected? The first to infect someone else.
  initial.infected.idx <- move$who[1]
  segments <- matrix(nrow=(nrow(move)+1), ncol=4)
  segment.cnt <- 0
  arrows <- matrix(nrow=nrow(move), ncol=4)
  arrow.idx <- 1
  state <- rep(NA, length=individual.cnt)
  # +1 because individuals are numbered from 0.
  state[initial.infected.idx+1] <- 0
  for (ridx in 1:nrow(birthdeath)) {
    event <- birthdeath$event[ridx]
    # +1 to go from zero based to one based
    who <- birthdeath$who[ridx]+1
    whom <- birthdeath$whom[ridx]+1
    when <- birthdeath$when[ridx]
    if (event==1) {
      arrows[arrow.idx, 1] <- when
      arrows[arrow.idx, 2] <- norder[who]
      arrows[arrow.idx, 3] <- when
      arrows[arrow.idx, 4] <- norder[whom]
      arrow.idx <- arrow.idx + 1

      if (is.na(state[whom])) {
        state[whom] <- when
      } # It's possible to have a bug move when you are already infested.
    } else if (event==6) {
      if (!is.na(state[whom])) {
        segment.cnt <- segment.cnt + 1
        segments[segment.cnt,1] <- state[whom]
        segments[segment.cnt,2] <- norder[whom]
        segments[segment.cnt,3] <- when
        segments[segment.cnt,4] <- norder[whom]
        state[whom] <- NA
      } else {
        print(paste("The state should be nonzero before it dies out.",
          "who", who, "whom", whom, "when", when, "state", state))
      }
    } else {
      print("no way you get here")
    }
  }
  for (cleanidx in 1:length(state)) {
    if (!is.na(state[cleanidx])) {
      segment.cnt <- segment.cnt + 1
      segments[segment.cnt,1] <- state[cleanidx]
      segments[segment.cnt,2] <- norder[cleanidx]
      segments[segment.cnt,3] <- plot.time
      segments[segment.cnt,4] <- norder[cleanidx]
    }
  }
  segments <- segments[1:segment.cnt,]

  pdf(plotfile)
  print("plot")
  plot(seq(0, plot.time, plot.time/individual.cnt), seq(0, individual.cnt),
    type="n", xlab="Time [years]", ylab="House", main="House Lifelines")
  print("segments")
  segments(segments[,1], segments[,2], segments[,3], segments[,4])
  print("arrows")
  arrows(arrows[,1], arrows[,2], arrows[,3], arrows[,4], length=0.05)
  dev.off()
}

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
datafile <- args[1]
plotfile <- args[2]
plot.verhulst.lives(datafile, plotfile)
