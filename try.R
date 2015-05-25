library("spatstat")
library("survival")
library("animation")
library("deldir")
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


houses.plot <- function(houses, streets, filename="points.pdf", crossings=NA) {
  house_cnt<-length(houses$x)
  if (!is.na(crossings)) {
    stopifnot(length(crossings)==house_cnt*house_cnt)
    print(crossings)
    p1x=vector(mode="numeric", house_cnt*house_cnt)
    p1y=vector(mode="numeric", house_cnt*house_cnt)
    p2x=vector(mode="numeric", house_cnt*house_cnt)
    p2y=vector(mode="numeric", house_cnt*house_cnt)
    conn_idx <- 1
    for (i in 1:(house_cnt-1)) {
      for (j in (i+1):house_cnt) {
        if (crossings[(i-1)*house_cnt+j]==0) {
          p1x[conn_idx] <- houses$x[i]
          p1y[conn_idx] <- houses$y[i]
          p2x[conn_idx] <- houses$x[j]
          p2y[conn_idx] <- houses$y[j]
          conn_idx <- conn_idx + 1
        } # else no crossing
      }
    }
    conn_cnt <- conn_idx - 1
    print(paste("There are", conn_cnt, "crossings"))
  }

  pdf(filename)
  plot(0:1, 0:1, type="n")
  for (hi in 1:house_cnt) {
    text(houses$x[hi], houses$y[hi], hi-1)
  }
  for (si in 1:length(streets$x)) {
    text(streets$x[si], streets$y[si], house_cnt+si-1)
  }
  #points(houses, col="black")
  dirsgs=streets$dirsgs
  segments(dirsgs$x1, dirsgs$y1, dirsgs$x2, dirsgs$y2, col="blue")
  if (!is.na(crossings)) {
    segments(p1x[1:conn_cnt], p1y[1:conn_cnt], p2x[1:conn_cnt], p2y[1:conn_cnt],
      col="red")
  }
  dev.off()
}


generate.streets <- function(block_cnt) {
  X <- rHardcore(block_cnt, 0.05, square(1))
  voronoi<-deldir(X, rw=c(0, 1, 0, 1))
  dsgs<-voronoi$dirsgs
  # They give points as x, y, not as a complex.
  # Find the ones that are distinct. We have an unusual advantage in
  # that we can calculate the number of distinct points there should
  # be, so a stable method is to find every point distance, sort them,
  # and pick a cutoff tolerance.
  print(dsgs)
  boundary_cnt <- 0
  edge_cnt <- nrow(dsgs)
  for (vcp_idx in 1:edge_cnt) {
    if (dsgs$bp1[vcp_idx]) {
      boundary_cnt <- boundary_cnt + 1
    }
    if (dsgs$bp2[vcp_idx]) {
      boundary_cnt <- boundary_cnt + 1
    }
  }
  pt_cnt <- 2*edge_cnt
  internal_cnt <- (pt_cnt-boundary_cnt)/3
  unique_cnt <- internal_cnt + boundary_cnt
  print(paste("unique", unique_cnt, "internal", internal_cnt,
    "pt_cnt", pt_cnt, "edge_cnt", edge_cnt, "boundary", boundary_cnt))
  x <- vector(mode="numeric", unique_cnt)
  y <- vector(mode="numeric", unique_cnt)
  p0 <- vector(mode="numeric", edge_cnt)
  p1 <- vector(mode="numeric", edge_cnt)

  all_pts <- ppp(c(dsgs$x1, dsgs$x2), c(dsgs$y1, dsgs$y2))
  dx <- pairdist(all_pts)
  dx <- as.vector(dx)
  dx <- dx[order(dx)]
  total <- pt_cnt*pt_cnt
  stopifnot(total==length(dx))
  zero_cnt <- boundary_cnt + 9*internal_cnt
  print(paste("expect zeros", zero_cnt, "find", sum(dx<0.000001)))
  #print(dx)
  tolerance=0
  if (dx[zero_cnt]>0) {
    # geometric mean, then squared to compare with squared distance.
    tolerance=dx[zero_cnt]*dx[zero_cnt+1]
  } else {
    # minimum size, squared to compare with squared distance.
    tolerance=0.1*(dx[zero_cnt+1])^2
  }
  print(paste("tolerance", tolerance))

  pt_idx <- 1
  find_or_make <- function(ax, ay, env=parent.frame()) {
    for (idx in 1:env$pt_idx) {
      if ((ax-x[idx])^2+(ay-y[idx])^2 < tolerance) {
        return(idx)
      }
    }
    add_idx<-env$pt_idx+1
    env$pt_idx<-add_idx
    env$x[add_idx]<-ax
    env$y[add_idx]<-ay
    # print(paste("inx", x, "iny", y))
    add_idx
  }
  x[pt_idx]<-dsgs$x1[1]
  y[pt_idx]<-dsgs$y1[1]
  for (e_idx in 1:edge_cnt) {
    p0p<-find_or_make(dsgs$x1[e_idx], dsgs$y1[e_idx])
    pt_idx<-max(p0p, pt_idx)
    p1p<-find_or_make(dsgs$x2[e_idx], dsgs$y2[e_idx])
    p0[e_idx]=p0p
    p1[e_idx]=p1p
    # print(paste("p0", p0p, "p1", p1p, "idx", pt_idx))
    # print(paste("x", x, "y", y))
  }
  stopifnot(pt_idx==unique_cnt)
  list(x=x, y=y, p0=p0, p1=p1, dirsgs=dsgs)
}

streets.write <- function(streets, basename="streets") {
  streetsxy<-list(x=streets$x, y=streets$y)
  write.csv(streetsxy, file=paste(basename, "xy.txt", sep=""))
  streetspt<-list(p0=streets$p0, p1=streets$p1)
  write.csv(streetspt, file=paste(basename, "pts.txt", sep=""))
  write.csv(streets$dirsgs, file=paste(basename, "dirsgs.txt", sep=""))
}

streets.read <- function(basename="streets") {
  xy <- read.csv(paste(basename, "xy.txt", sep=""))
  pts <- read.csv(paste(basename, "pts.txt", sep=""))
  dirsgs <- read.csv(paste(basename, "dirsgs.txt", sep=""))
  list(x=xy$x, y=xy$y, p0=pts$p0, p1=pts$p1, dirsgs=dirsgs)
}

# Save houses and street blocks to files.
landscape.generate <- function(house_cnt, block_cnt, streetsbase="streets") {
  houses<-rHardcore(house_cnt, 0.01, square(1))

  streets<-generate.streets(block_cnt)
  houses.plot(houses, streets, filename="nocross.pdf")

  write.csv(houses, file="houses.txt")
  streets.write(streets, basename=streetsbase)

  list(houses=houses, streets=streets)
}

landscape.crossings <- function( housesfile="houses.txt",
    streetsfile="streets" ) {
  houses <- read.csv(housesfile)
  streets <- streets.read(streetsfile)
  crossings<-intersections(houses$x, houses$y, streets$x, streets$y,
    streets$p0, streets$p1)

  houses.plot(houses, streets, crossings=crossings, filename="landscape.pdf")
  write.csv(crossings, file="crossings.txt")
  list(houses=houses, crossings=crossings)
}


sir.generate <- function(run_cnt=1, outfile="sir.h5",
    housesfile="houses.txt", crossingsfile="crossings.txt") {
  houses.table <- read.csv(housesfile)
  print(houses.table)
  houses <- ppp(houses.table$x, houses.table$y)
  print(houses)
  crossings <- read.csv(crossingsfile)$x

  dx<-pairdist(houses)
  start<-center.point(houses)
  print(paste("start is", start))
  p<-c(individual_cnt=dim(dx)[[1]], seed=33333, N0=1.0, beta0=0.3, beta1=0.01,
    beta2=0.1, gamma=0.1, cutoff=0.1, growthrate=4.5, carrying=1000.0,
    runs=run_cnt, initial=start, streetfactor=0.4
    )
  # fully infected takes n days to infect neighbor
  ndays <- 30
  max_hazard <- 365/ndays
  p["beta0"] <- max_hazard/(p["growthrate"]*p["carrying"])
  p["beta1"] <- p["beta0"]/20
  print(p)

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
  simple_hazard(dx, crossings, p, save.to.file)
  print("ran")
}

#  list("locations"=X, "events"=callback_env$results[[1]])



verhulst.generate <- function(run_cnt=1, outfile="verhulst.h5",
  housesfile="houses.txt", crossingsfile="crossings.txt") {
  houses.table <- read.csv(housesfile)
  print(houses.table)
  houses <- ppp(houses.table$x, houses.table$y)
  print(houses)
  crossings <- read.csv(crossingsfile)$x

  dx<-pairdist(houses)
  start<-center.point(houses)
  print(paste("start is", start))

  p<-c(individual_cnt=dim(dx)[[1]], seed=33333, birth=0.1, death=0.005,
    carrying=1000.0, move0=0.3, move1=0.001, gamma=0.1, initial_bug_cnt=50,
    cutoff=0.1, alpha1=1.0, alpha2=0.0, streetfactor=0.4, runs=run_cnt,
    initial=start)

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


create.web.movie <- function(locations, events) {
  initial_marking<-factor(rep('s', locations$n), levels=c('s', 'i', 'n'))
  locations<-ppp(locations$x, locations$y, marks=initial_marking)
  event_cnt<-dim(events)[[1]]
  cat("There were", event_cnt, "events\n")

  saveHTML({
     for (i in 1:event_cnt) {
        event_tag<-events$event[[i]]
        if (event_tag==1 || event_tag==2) {
          locations$marks[[events$whom[[i]]+1]]='i'
        } else if (event_tag==3) {
          locations$marks[[events$whom[[i]]+1]]='n'
        } else if (event_tag==4) {
          # Returning to susceptible from infected.
          locations$marks[[events$whom[[i]]+1]]='s'
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


get.single.trajectory <- function(filename, which=1) {
  env <- new.env(parent=emptyenv())
  env$idx <- 1
  get.trajectory <- function(locations, arg) {
    if (env$idx==which) {
      env$locations <- locations
      env$trajectory <- arg
    }
    env$idx <- env$idx + 1
  }
  foreach.trajectory(filename, get.trajectory, 1)
  list(locations=env$locations, trajectory=env$trajectory)
}

to.filename <- function(human_name) {
  gsub("[ %%]", "_", human_name)
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
  foreach.trajectory(filename, get.times)
  all_times <- unlist(env$times)
  surv.object <- Surv(all_times, rep(1, length(all_times)))
  x.fit<-survfit(surv.object ~ 1)
  name <- "Infection Times"
  pdf(paste(to.filename(name), ".pdf", sep=""))
  plot(x.fit)
  title(paste("Survival for", name),
    xlab = "Time [years]", ylab = "Survival Fraction")
  dev.off()
}

# res<-bugtest(1000)
# res<-sirgenerate(10, 5, 1)
#landscape.generate(50, 10)

# warnings()
# foreach.trajectory("1000.h5", end.times, 3)
# infection.times.hazard("z.h5")
#t<-get.single.trajectory("z.h5", 1)
#create.web.movie(t$locations, t$trajectory)
