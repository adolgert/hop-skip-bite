library("spatstat")
library("survival")
library("animation")
library("hopskip")

source("file_format.R")

repulsion_matrix <- function(n) {
  X<-rHardcore(n, 0.02, square(1))
  distance_matrix<-pairdist(X)
}


center.point <- function(spatdata) {
  center <- ppp(0.5, 0.5)
  which.min(crossdist(spatdata, center))
}


sirgenerate <- function(house_cnt, run_cnt) {
  X<-rHardcore(house_cnt, 0.02, square(1))
  dx<-pairdist(X)
  start<-center.point(X)
  p<-c(individual_cnt=dim(dx)[[1]], seed=33333, N0=1.0, beta0=0.3, beta1=0.01,
    beta2=0.1, gamma=0.1, cutoff=0.1, growthrate=4.5, carrying=1000.0,
    runs=run_cnt, initial=start
    )

  outfile="z.h5"
  create_file(outfile)
  write_locations(outfile, X)
  callback_env<-new.env(parent=emptyenv())
  callback_env$results=list()
  save.to.file <- function(arg) {
    write_events(outfile, arg)
  }
  add.to.results<-function(arg) {
    callback_env$results[[length(callback_env$results)+1]]<- arg
  }
  
  simple_hazard(dx, p, save.to.file)
}

#  list("locations"=X, "events"=callback_env$results[[1]])



bugtest <- function(n) {
  X<-rHardcore(n, 0.02, square(1))
  dx<-pairdist(X)
  p<-c(individual_cnt=dim(dx)[[1]], seed=33333, birth=0.1, death=0.005,
    carrying=1000.0, move0=0.3, move1=0.001, gamma=0.1, initial_bug_cnt=50,
    cutoff=0.1)
  list("locations"=X, "events"=bugs(dx, p))
}


create.web.movie <- function(single_result) {
  initial_marking<-factor(rep('s', res$locations$n), levels=c('s', 'i', 'n'))
  locations<-ppp(res$locations$x, res$locations$y, marks=initial_marking)
  event_cnt<-dim(res$events)[[1]]
  cat("There were", event_cnt, "events\n")

  saveHTML({
     for (i in 1:event_cnt) {
        event_tag<-res$events$event[[i]]
        if (event_tag==1 || event_tag==2) {
          locations$marks[[res$events$whom[[i]]+1]]='i'
        } else if (event_tag==3) {
          locations$marks[[res$events$whom[[i]]+1]]='n'
        } else if (event_tag==4) {
          # Returning to susceptible from infected.
          locations$marks[[res$events$whom[[i]]+1]]='s'
        }
        plot(locations, cols=c("blue", "red", "black"))
     }
    }, img.name = "infection", title="spatial infection",
    interval = 0.1, nmax = event_cnt,
    ani.width = 600, ani.height = 600, loop = TRUE)
}


end.times <- function(locations, events) {
  cnt=nrow(events)
  print(events$when[[cnt]])
}


single.trajectory <- function(filename, which=-1) {
  env <- new.env(parent=emptyenv())
  get.trajectory <- function(locations, arg) {
    env$locations <- locations
    env$trajectory <- arg
  }
  foreach.trajectory(filename, get.trajectory, 1)
  list(locations=env$locations, trajectory=env$trajectory)
}


infection.times.hazard <- function(filename) {
  trajectory_cnt <- trajectory.count(filename)
  env <- new.env(parent=emptyenv())
  env$times<-vector(mode="list", length=trajectory_cnt)
  env$idx<-1
  get.times<-function(locations, arg) {
    cnt=nrow(locations)-1 # exclude first infection
    print(paste("number of events", nrow(arg)))

    time_to_infect <- subset(arg, event==1 | event==2)$when
    print(paste("number infections", length(time_to_infect)))
    print(env$idx)
    env$times[[env$idx]] <- time_to_infect
    env$idx <- env$idx + 1
  }
  foreach.trajectory(filename, get.times, 3)
}

# res<-bugtest(1000)
# res<-sirgenerate(1000, 10)

# foreach.trajectory("1000.h5", end.times, 3)
infection.times.hazard("z.h5")
