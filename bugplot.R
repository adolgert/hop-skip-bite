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
      pointsize=144)
  } else {
    print(paste("Unknown type of plot file, ending", ending))
  }
  plot(seq(0, 2, 2/12), seq(0, 1200, 100), type="n",
    xlab="Time [years]", ylab="Individuals",
    main="Verhulst Model Bug Counts")
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
