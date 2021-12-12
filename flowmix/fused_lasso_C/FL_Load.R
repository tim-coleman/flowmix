setwd("C:/Users/drain/Box Sync/USC_Stuff/Featureless_Flowmix/fused_lasso_C")
dyn.load("tf_dp_ryan_orig.dll")

getLoadedDLLs()



prox = function(z, lam) {
  if (!is.loaded("prox_dp_R")) {
    dyn.load("prox_R.so")
  }
  
  o = .C("prox_dp_R",
         n=as.integer(length(z)),
         z=as.double(z),
         lam=as.double(lam),
         beta=as.double(numeric(length(z))),
         dup=FALSE)
  
  return(o$beta)
}


set.seed(0)
z <- rnorm(100)
lam <- 0.5


out <- prox(z, lam)
plot(z)
lines(out, col = 'blue')


out <- prox(z, 10*lam)
plot(z)
lines(out, col = 'blue')


run_prox <- function(n = 1e10){
  z <- rnorm(n)
  lam <- 0.5
  out <- prox(z, lam)
}

library(microbenchmark)
microbenchmark(run_prox)
