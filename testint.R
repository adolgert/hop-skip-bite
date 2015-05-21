# Test line intersection

library("spatstat")
library("deldir")
library("hopskip")


center.point <- function(spatdata) {
  center <- ppp(0.5, 0.5)
  which.min(crossdist(spatdata, center))
}

gen_points <- function(house_cnt) {
  X<-rHardcore(house_cnt, 0.02, square(1))
  dx<-pairdist(X)
  X
}

gen_streets <- function(block_cnt) {
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


ragged.read <- function(name, event_cnt) {
  lscan <- as.numeric(scan(file=name))
  print(lscan)
  llist<-list()
  lidx <- 1
  for (li in 1:event_cnt) {
    cnt <- lscan[lidx]
    print(paste("reading", cnt, "lidx", lidx))
    if (cnt>0) {
      print(lscan[(lidx+1):(lidx+cnt)]+1)
      llist[[li]] <- lscan[(lidx+1):(lidx+cnt)]
    }
    lidx <- lidx + cnt + 1
  }
  llist
}

plot.graphical <- function() {
  points <- read.table("points.txt")
  print(points)
  intpoints <- read.table("intpoints.txt")
  print(intpoints)
  lines <- read.table("lines.txt")
  print(lines)
  curp <- read.table("current.txt")
  print(curp)
  event_cnt <- length(curp$V1)
  llist <- ragged.read("l.txt", event_cnt)
  ulist <- ragged.read("u.txt", event_cnt)
  clist <- ragged.read("c.txt", event_cnt)
  tlist <- ragged.read("tstatus.txt", event_cnt)

  for (gidx in 1:event_cnt) {
    png(paste("plot", sprintf("%03d", gidx), ".png", sep=""),
      width=800, height=800, res=120)
    plot(0:1, 0:1, type="n")
    for (hi in 1:length(points$V1)) {
      text(points$V1[hi], points$V2[hi], hi-1)
    }
    px<-points$V1[curp$V1[gidx]+1]
    py<-points$V2[curp$V1[gidx]+1]

    l0x<-points$V1[lines$V1+1]
    l0y<-points$V2[lines$V1+1]
    l1x<-points$V1[lines$V2+1]
    l1y<-points$V2[lines$V2+1]

    segments(l0x, l0y, l1x, l1y, col="grey")
    points(px, py, col="red")

    segments(0, py, 1, py, col="green")
    if (!is.na(tlist[gidx])) {
      tx<-as.numeric(tlist[gidx][[1]])
      print(paste("tx is", tx))
      segments(tx, rep(py-0.05, length(tx)), tx, rep(py+0.05, length(tx)),
        col="green")
    }
    intx<-intpoints$V1
    inty<-intpoints$V2
    abovex<-intx[inty>py-1e-7]
    abovey<-inty[inty>py-1e-7]
    print(abovex)
    print(abovey)
    points(abovex, abovey, col="blue")

    dev.off()
  }
}


make.plot <- function(houses, streets, crossings) {
  house_cnt<-length(houses$x)
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

	pdf("points.pdf")
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
  segments(p1x[1:conn_cnt], p1y[1:conn_cnt], p2x[1:conn_cnt], p2y[1:conn_cnt],
    col="red")
	dev.off()
}

gen.world <- function() {
  houses<-gen_points(100)
  streets<-gen_streets(5)
  print(paste("houses", houses$x, houses$y))
  crossings<-intersections(houses$x, houses$y, streets$x, streets$y,
    streets$p0, streets$p1)
  make.plot(houses, streets, crossings)
}

gen.sample <- function() {
  houses<-list(x=c(0.2, 0.8), y=c(0.2, 0.8))
  streets<-list(x=c(0.3, 0.7), y=c(0.85, 0.4), p0=c(1), p1=c(2))
  crossings<-intersections(houses$x, houses$y, streets$x, streets$y,
    streets$p0, streets$p1)
  make.plot(houses, streets, crossings)
}

gen.world()
#plot.graphical()

