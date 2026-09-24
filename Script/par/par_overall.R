TM <- 1 # maximal time
library(simcausal); options(simcausal.verbose = FALSE)
library(tibble)
library(dplyr)
library(purrr)
library(ggplot2)
library(randomForest)

deltaA = function(L) {
  tanh(0.7 + 0.4 * sin(0.25 * L))
}
OPA = function(L, U) {
  exp(-0.25 * U)
}
pi = function(OPA, deltaA) {
  numerator1 <- OPA * (2 - deltaA) + deltaA
  numerator2 <- (OPA * (deltaA - 2) - deltaA)^2 + 4 * OPA * (1 - OPA) * (1 - deltaA)
  denumerator <- 2 * (OPA - 1)
  return( (numerator1 - sqrt(numerator2)) / denumerator )
}

D <- DAG.empty() +
  node("L",t=0, distr = "rnorm", mean = 0, sd = 1.5) +
  node("U",t=0, distr = "rnorm", mean = 0, sd = 1.5) +
  node("Z",t=0, distr = "rbern", prob = plogis(-0.5 * L[t])) +
  node("A",t=0, distr = "rbern", prob = pi(OPA(L[t],U[t]),deltaA(L[t]))+Z[t]*deltaA(L[t]))+
  
  node("L",t=1:TM, distr = "rnorm", mean = (A[t-1]-0.5)+0.5*L[t-1]+0.3*U[t-1], sd = 0.5) +
  node("U",t=1:TM, distr = "rnorm", mean = (A[t-1]-0.5)+0.5*U[t-1]+0.3*L[t], sd = 0.5) +
  node("Z",t=1:TM, distr = "rbern", prob = plogis(-0.5*L[t]-0.3*(A[t-1]-0.5)-0.1*(Z[t-1]-0.5))) +
  node("A",t=1:TM, distr = "rbern", prob = pi(OPA(L[t],U[t]),deltaA(L[t]))+Z[t]*deltaA(L[t]))+
  node("Y",t=TM+1,distr = "rnorm", mean = (A[t-1]-0.5)+2*L[t-1]+U[t-1], sd = 0.5)

D <- set.DAG(D, vecfun=c("OPA","pi","deltaA"))
# plotDAG(D, xjitter = 1)
