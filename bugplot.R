library(spatstat) # for ppp data structure
library(scales) # for alpha
source("file_format.R")

plot.all.bug.trajectory <- function(datafile, plotfile) {
  print(paste("data file is", datafile))
  ending<-substring(plotfile, nchar(plotfile)-3)
  if (ending==".pdf") {
    print(paste("plot file is pdf", plotfile))
    pdf(plotfile)
  } else if (ending==".png") {
    print(paste("plot file is png", plotfile))
    png(plotfile, width=6000, height=6000, units="px",
      pointsize=200)
  } else {
    print(paste("Unknown type of plot file, ending", ending))
  }
  plot.time <- 2
  plot(seq(0, plot.time, 2/12), seq(0, 1200, 100), type="n",
    xlab="Time [years]", ylab="Individuals",
    main="Verhulst Model Bug Counts")

  callback.env <- new.env(parent=emptyenv())
  # resolution of average line.
  time.points <- seq(0, plot.time, plot.time/201)
  callback.env$cumulant <- numeric(length(time.points))
  callback.env$dieoff <- numeric(length(time.points))
  callback.env$trajectory.cnt <- 0
  callback.env$die.cnt <- 0
  get.trajectory <- function(locations, events) {
    with.counts <- subset(events, event!=1)
    lines(with.counts$when, with.counts$who,
      col=alpha("green", 0.1))

    bug.cnt <- 1
    bug.cnt.prev <- bug.cnt
    event.idx <- 0
    event.time <- 0
    event.cnt <- nrow(events)
    print(paste("events", event.cnt, events$what))
    cumulant <- numeric(length(time.points))
    for (tidx in 1:length(time.points)) {
      time <- time.points[tidx]
      while (event.idx<event.cnt & event.time<time) {
        event.idx <- event.idx + 1
        event.time <- events$when[event.idx]
        if (events$event[event.idx] == 5 | events$event[event.idx] == 6) {
          bug.cnt.prev <- bug.cnt
          bug.cnt <- events$who[event.idx]
        }
      }
      cumulant[tidx] <- bug.cnt.prev
    }
    if (bug.cnt.prev>100) {
      callback.env$cumulant <- callback.env$cumulant + cumulant
      callback.env$trajectory.cnt <- callback.env$trajectory.cnt + 1
    } else {
      callback.env$dieoff <- callback.env$dieoff + cumulant
      callback.env$die.cnt <- callback.env$trajectory.cnt + 1
    }
  }
  print("foreach.trajectory")
  foreach.trajectory(datafile, get.trajectory)

  callback.env$cumulant <- callback.env$cumulant/callback.env$trajectory.cnt
  callback.env$dieoff <- callback.env$dieoff/callback.env$die.cnt
  write.csv(callback.env$cumulant, file="verhulst_single_cumulant.csv")
  lines(time.points, callback.env$cumulant, col="black", lwd=15)
  lines(time.points, callback.env$dieoff, col="blue", lwd=10)
  print(paste("trajectories survived", callback.env$trajectory.cnt,
    "trajectories died", callback.env$die.cnt))
  dev.off()
}


options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
datafile <- args[1]
plotfile <- args[2]
plot.all.bug.trajectory(datafile, plotfile)

