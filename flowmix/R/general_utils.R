
# packages ----------------------------------------------------------------


library(scales)
library(dplyr)
library(ggplot2)
library(tidyr)
library(Matrix)
library(plotly)
library(CVXR)

# difference matrix generation --------------------------------------------


gen_diff_mat <- function(n, l){
  get_D1 <- function(t) {do.call(rbind, lapply(1:(t-1), FUN = function(x){
    v <- rep(0, t)
    v[x] <- -1
    v[x+1] <- 1
    v
  }))}
  if(l == 0){
    return(diag(rep(1,n)))
  }
  if(l == 1){
    return(get_D1(n))
  }
  if(l > 1){
    D <- get_D1(n)
    for(k in 1:(l-1)){
      D <- get_D1(n-k) %*% D
    }
    return(D)
  }
}




# data generation ---------------------------------------------------------


y_sine_data <- function(TT, nt){
  TT.scaled <- seq(-1, 1, length.out = TT)
  y <- lapply(TT.scaled, function(t){
    y1 <- replicate(nt, rnorm(n = 1, mean = sin(4*pi*t + pi/3), sd = 0.25))
    y2 <- replicate(nt, rnorm(n = 1, mean = 5 + sin(2*pi*t), sd = 0.25))
    return(cbind(y1, y2))
  })
  return(y)
}



# Centering a ylist -------------------------------------------------------

weight_ylist <- function(iclust, resp, resp.sum, ylist){
  
  ## Setup
  dimdat = ncol(ylist[[1]])
  TT = length(ylist)
  
  ## All weighted data
  weighted_ylist = Map(function(myresp, y){
    myresp[,iclust, drop = TRUE] * y
  }, resp, ylist)
  
  ## All weighted data SUMs
  weighted_ysum = lapply(weighted_ylist, colSums)
  
  ## Grand mean of data
  resp.sum.thisclust = sum(resp.sum[,iclust])
  ybar = Reduce("+", weighted_ysum) / resp.sum.thisclust
  
  ## Centered weighted ylist
  weighted_ybar = resp.sum[,iclust] * matrix(ybar, TT, dimdat, byrow=TRUE)
  ycentered = do.call(rbind, weighted_ysum) - weighted_ybar##sweep(do.call(cbind, centered_y), 1, ybar)
  return(list(weighted_ylist = weighted_ylist,
              weighted_ysum = weighted_ysum,
              ybar = ybar,
              ycentered = t(ycentered)))
}



# Outer convergence check -------------------------------------------------

outer_converge <- function(objectives){
  consec = 4
  if(length(objectives) < consec){
    return(FALSE)
  } else {
    mytail = utils::tail(objectives, consec)
    rel_diffs = mytail[1:(consec-1)]/mytail[2:consec]
    
    
    return(all(abs(rel_diffs) - 1 < 1E-5))
  }
}

