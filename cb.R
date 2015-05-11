# Testing how to do a callback
library(hopskip)

# First do a regular callback without C++.

list_env<-new.env(parent=emptyenv())
list_env$alist=list()

findlist<-function(arg) {
  print(environment(alist))
  print(paste("exists", exists("alist")))
  print(length(alist))
  list_env$alist[[length(list_env$alist)+1]]<-arg
  print(alist)
}

findlist(7)
findlist(8)
findlist(9)
print(list_env$alist)


callback_maker <- function(env=parent.frame()) {
  get<-function() env$alist
  set<-function(arg) env$alist[[length(env$alist)+1]]<- arg
  list(set=set, get=get)
}

a=callback_maker()
print(typeof(a$set))
a$set(7)
a$set(5)
a$set(3)
print(a$get())

TestCallback(a$set)
print(alist)
