
repulsion_matrix <- function(n) {
  library("spatstat")
  X<-rHardcore(n, 0.02, square(1))
  distance_matrix<-pairdist(X)
}


sirtest <- function(n) {
  dx<-repulsion_matrix(n)
  p<-c(individual_cnt=dim(dx)[[1]], seed=33333, beta0=0.1, beta1=0.1,
    beta2=0.1, gamma=0.1)
  c(points=dx, events=rcpp_hello_world(dx, p))
}

