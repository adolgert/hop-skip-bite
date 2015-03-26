library("spatstat")
library("animation")
library("hopskip")

repulsion_matrix <- function(n) {
  X<-rHardcore(n, 0.02, square(1))
  distance_matrix<-pairdist(X)
}


sirtest <- function(n) {
  X<-rHardcore(n, 0.02, square(1))
  dx<-pairdist(X)
  p<-c(individual_cnt=dim(dx)[[1]], seed=33333, beta0=0.3, beta1=0.01,
    beta2=0.1, gamma=0.1, cutoff=0.1)
  list("locations"=X, "events"=simple_hazard(dx, p))
}

res<-sirtest(400)
initial_marking<-factor(rep('s', res$locations$n), levels=c('s', 'i', 'n'))
locations<-ppp(res$locations$x, res$locations$y, marks=initial_marking)
event_cnt<-dim(res$events)[[1]]


saveHTML({
   for (i in 1:event_cnt) {
      event_tag<-res$events$event[[i]]
      if (event_tag==1 || event_tag==2) {
        locations$marks[[res$events$who[[i]]+1]]='i'
      } else if (event_tag==3) {
        locations$marks[[res$events$who[[i]]+1]]='n'
      }
      plot(locations, cols=c("blue", "red", "black"))
   }
  }, img.name = "infection", title="spatial infection",
  interval = 0.1, nmax = event_cnt,
  ani.width = 600, ani.height = 600, loop = TRUE)
