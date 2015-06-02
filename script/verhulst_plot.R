library(spatstat)
library(ggplot2)
source("file_format.R")


plot.verhulst.set <- function(datafile, plotfile) {
  trajectory.info <- single.trajectory(datafile, "dset1")
  locations <- trajectory.info$locations
  df <- trajectory.info$events
  out.df <- data.frame(when=as.numeric(character()), who=character(),
    bugs=as.numeric(character()))

  individual.cnt <- length(unique(df$whom))
  plot.time <- max(df$when)
  division.cnt <- 201
  time.points <- seq(0, plot.time, plot.time/division.cnt)
  for (individual.idx in unique(df$whom)) {
    ind.df <- subset(df, whom==individual.idx)

    bug.cnt <- 0
    bug.cnt.prev <- bug.cnt
    event.idx <- 0
    event.time <- 0
    event.cnt <- nrow(ind.df)
    cumulant <- numeric(length(time.points))
    for (tidx in 1:length(time.points)) {
      time <- time.points[tidx]
      while (event.idx<event.cnt & event.time<time) {
        event.idx <- event.idx + 1
        event.time <- ind.df$when[event.idx]
        if (ind.df$event[event.idx] == 5 | ind.df$event[event.idx] == 6) {
          bug.cnt.prev <- bug.cnt
          bug.cnt <- ind.df$who[event.idx]
        }
      }
      cumulant[tidx] <- bug.cnt.prev
    }
    print(paste("read", individual.idx))
    add.df <- data.frame(when=time.points, bugs=cumulant,
      who=rep(as.character(individual.idx), length(time.points)))
    out.df <- rbind(out.df, add.df)
  }
  print(paste("output is", plotfile))
  qplot(data=out.df, when, bugs, geom="line", colour=who) +
    labs(x="Time [years]", y="Individuals", colour="",
      title="Vector Trajectory")
  ggsave(file=plotfile)
}


options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
datafile <- args[1]
plotfile <- args[2]
plot.verhulst.set(datafile, plotfile)
