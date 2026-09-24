TM <- 1 # maximal time
library(simcausal); options(simcausal.verbose = FALSE)
library(tibble)
library(dplyr)
library(mgcv)
library(purrr)
library(ggplot2)

# install.packages("remotes")
# remotes::install_github("tedwestling/survSuperLearner")
# library(survSuperLearner)

compliance_map <- function(Z, R) {
  if (length(Z) != length(R)) {
    stop("Z and R must have the same length.")
  }
  
  # 统一处理大小写
  R <- tolower(R)
  
  A <- rep(NA_integer_, length(Z))
  
  A[R == "nt"] <- 0
  A[R == "at"] <- 1
  A[R == "co"] <- Z[R == "co"]
  A[R == "de"] <- 1 - Z[R == "de"]
  
  return(A)
}



D <- DAG.empty() +
  node("L",t=0, distr = "rnorm", mean = 0, sd = 1.5) +
  node("U",t=0, distr = "rnorm", mean = 0, sd = 1.5) +
  node("Z",t=0, distr = "rbern", prob = plogis(-0.5 * L[t])) +
  node("scoreCo",t=0, distr = "rconst", const = exp(2+0.5 * L[t] - 0.5 * U[t])) +
  node("scoreAt",t=0, distr = "rconst", const = exp(-0.5 * L[t] + 0.5 * U[t])) +
  node("scoreNt",t=0, distr = "rconst", const = 1) +
  node("scoreSum",t=0, distr = "rconst", const = scoreCo[t] + scoreAt[t] + scoreNt[t])+
  node("Rtype",t=0, distr = "rcat.b1",
       probs = cbind(scoreCo[t] / scoreSum[t],
                     scoreAt[t] / scoreSum[t],
                     scoreNt[t] / scoreSum[t]))+
  node("R",t=0, distr = "rconst",
       const = ifelse(Rtype[t] == 1, "co",
                      ifelse(Rtype[t] == 2, "at", "nt")))+
  node("A",t=0, distr = "rconst", const = compliance_map(Z[t], R[t]))+
  
  node("L",t=1:TM, distr = "rnorm", mean = (A[t-1]-0.5)+0.5*L[t-1]+0.3*U[t-1], sd = 0.5) +
  node("U",t=1:TM, distr = "rnorm", mean = (A[t-1]-0.5)+0.5*U[t-1]+0.3*L[t], sd = 0.5) +
  node("Z",t=1:TM, distr = "rbern", prob = plogis(-0.5*L[t]-0.3*(A[t-1]-0.5)-0.1*(Z[t-1]-0.5))) +
  node("scoreCo",t=1:TM, distr = "rconst", const = exp(2+0.5 * L[t] - 0.5 * U[t]-(A[t-1]-0.5))) +
  node("scoreAt",t=1:TM, distr = "rconst", const = exp(-0.5 * L[t] + 0.5 * U[t]-(A[t-1]-0.5))) +
  node("scoreNt",t=1:TM, distr = "rconst", const = 1) +
  node("scoreSum",t=1:TM, distr = "rconst", const = scoreCo[t] + scoreAt[t] + scoreNt[t])+
  node("Rtype",t=1:TM, distr = "rcat.b1",
       probs = cbind(scoreCo[t] / scoreSum[t],
                     scoreAt[t] / scoreSum[t],
                     scoreNt[t] / scoreSum[t]))+
  node("R",t=1:TM, distr = "rconst",
       const = ifelse(Rtype[t] == 1, "co",
                      ifelse(Rtype[t] == 2, "at", "nt")))+
  node("A",t=1:TM, distr = "rconst", const = compliance_map(Z[t], R[t]))+
  node("Y",t=TM+1,distr = "rnorm", mean = (A[t-1]-0.5)+2*L[t-1]+U[t-1], sd = 0.5)

D <- set.DAG(D,vecfun = c("cbind","compliance_map"))

# sim(D, n = 100)
# plotDAG(D, xjitter = 1)
