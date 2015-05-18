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
  pt_cnt=max(dsgs$ind1)
  pt_cnt=max(pt_cnt, max(dsgs$ind2))
  print(pt_cnt)
  x <- vector(mode="numeric", pt_cnt)
  y <- vector(mode="numeric", pt_cnt)
  for (vcp_idx in 1:nrow(dsgs)) {
    x[dsgs$ind1[vcp_idx]] <- dsgs$x1[vcp_idx]
    y[dsgs$ind1[vcp_idx]] <- dsgs$y1[vcp_idx]
    x[dsgs$ind2[vcp_idx]] <- dsgs$x2[vcp_idx]
    y[dsgs$ind2[vcp_idx]] <- dsgs$y2[vcp_idx]
  }
  list(x=x, y=y, p0=dsgs$ind1, p1=dsgs$ind2, delsgs=dsgs)
}




make.plot <- function(houses, streets, crossings) {
  house_cnt<-length(houses$x)
  stopifnot(length(crossings)==house_cnt*house_cnt)
  p1x=vector(mode="numeric", house_cnt*house_cnt)
  p1y=vector(mode="numeric", house_cnt*house_cnt)
  p2x=vector(mode="numeric", house_cnt*house_cnt)
  p2y=vector(mode="numeric", house_cnt*house_cnt)
  conn_idx<-1
  for (i in 1:(house_cnt-1)) {
    for (j in (i+1):house_cnt) {
      if (crossings[i*house_cnt+j]>0) {
        p1x[conn_idx] <- houses$x[i]
        p1y[conn_idx] <- houses$y[i]
        p2x[conn_idx] <- houses$x[j]
        p2y[conn_idx] <- houses$y[j]
        conn_idx <- conn_idx + 1
      }
    }
  }
  conn_cnt <- conn_idx - 1
  print(paste("There are", conn_cnt, "crossings"))

	pdf("points.pdf")
	plot(0:1, 0:1, type="n")
	points(houses, col="black")
  delsgs=streets$delsgs
  segments(delsgs$x1, delsgs$y1, delsgs$x2, delsgs$y2, col="blue")
  segments(p1x, p1y, p2x, p2y, col="red")
	dev.off()
}

houses<-gen_points(50)
streets<-gen_streets(5)
crossings<-intersections(houses$x, houses$y, streets$x, streets$y,
  streets$p0, streets$p1)
make.plot(houses, streets, crossings)
