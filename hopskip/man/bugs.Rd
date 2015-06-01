\name{bugs}
\alias{bugs}
\docType{package}
\title{
Discrete Vector Movement in Continuous Time
}
\description{
Simulates movement of vectors on a landscape
}
\usage{
bugs(distance.matrix, crossings.matrix, parameters, output.function)
}
\examples{
\dontrun{
verhulst.generate <- function(run_cnt=1, outfile="verhulst.h5",
    housesfile="houses.txt", crossingsfile="crossings.txt", max_time=5.0) {
  houses.table <- read.csv(housesfile)
  print(houses.table)
  houses <- ppp(houses.table$x, houses.table$y)
  print(houses)
  crossings <- read.csv(crossingsfile)$x

  dx<-pairdist(houses)
  start<-center.point(houses)
  print(paste("start is", start))
  # Translate from deterministic to parameters for Verhulst.
  K1=1000 # carrying. We kind of know carrying.
  R0=1.5 # lambda/mu. Um, R0>1.15 is considered large.
  r=4.5 # logistic growth rate that puts it in a 2-years to carrying.
  mu=r/(R0-1) # so find death rate
  lambda=mu+r # find birth rate
  N=K1*R0/(R0-1) # Find N, a parameter for Verhulst, not the carrying.
  # assume alpha1=1 and alpha2=0.


  p<-c(individual_cnt=dim(dx)[[1]], seed=33333, birth=lambda, death=mu,
    carrying=N, move0=0.3, move1=0.001, gamma=0.1, initial_bug_cnt=1,
    cutoff=0.1, alpha1=1.0, alpha2=0.0, streetfactor=0.4, runs=run_cnt,
    initial=start, max_time=max.time)

  create_file(outfile)
  write_locations(outfile, houses)
  callback_env<-new.env(parent=emptyenv())
  callback_env$results=list()
  save.to.file <- function(arg) {
    write_events(outfile, arg, p)
  }
  add.to.results<-function(arg) {
    callback_env$results[[length(callback_env$results)+1]]<- arg
  }
  
  print("running")
  bugs(dx, crossings, p, save.to.file)
  print("ran")
}
}
}
