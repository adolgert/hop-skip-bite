library("spatstat")

repulsion_matrix <- function() {
  X<-rHardcore(10, 0.02, square(1))
  distance_matrix<-pairdist(X)
}


sirtest <- function() {
  dx<-repulsion_matrix()
  p<-c(individual_cnt=dim(dx)[[1]], seed=33333, beta0=0.1, beta1=0.1,
    beta2=0.1, gamma=0.1)
  rcpp_hello_world(dx, p)
}
