# Test line intersection

library("spatstat")


center.point <- function(spatdata) {
  center <- ppp(0.5, 0.5)
  which.min(crossdist(spatdata, center))
}

gen_points <- function(house_cnt) {
  X<-rHardcore(house_cnt, 0.02, square(1))
  dx<-pairdist(X)
  start<-center.point(X)
  write.table(X, "points.txt", sep=" ", row.names=FALSE, col.names=FALSE)
  cnt=length(X)
  print(cnt)
  segments<-matrix(nrow=cnt*(cnt-1)/2, ncol=2)
  print(dim(segments))
  running_idx <- 1
  for (i in 1:(cnt-1)) {
  	for (j in (i+1):cnt) {
  		print(paste(running_idx, i, j))
  	  segments[running_idx, 1] <- i-1
  	  segments[running_idx, 2] <- j-1
  	  running_idx <- running_idx + 1
  	}
  }
  write.table(segments, "lines.txt", sep=" ", row.names=FALSE, col.names=FALSE)
  list(points=X, segments=segments)
}

read.intersections <- function() {
  # +1 to convert from 0-based indexing of vertices.
  ll<-read.table("intlines.txt") + 1
  ip<-read.table("intpoints.txt")
  pxy<-ppp(ip$V1, ip$V2)
  list(points=pxy, segments=ll)
}


make.plot <- function(ps, is) {
	pdf("points.pdf")
	plot(0:1, 0:1, type="n")
	points(ps$points, col="black")
	points(is$points, col="red")
	dev.off()
}

ps<-gen_points(10)
system("./segint")
is<-read.intersections()
make.plot(ps, is)
