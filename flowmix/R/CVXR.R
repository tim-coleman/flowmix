

# main function -----------------------------------------------------------


cvxr_tf <- function(y, resp, lambda, l, 
                    Sigma_inv, 
                    thresh = 1E-8,
                    maxdev = NULL,
                    dimdat,
                    N,
                    ecos_thresh = 1E-8,
                    scs_eps = 1E-5){
  
  ## Define dimensions
  TT = length(y)
  
  ## Responsibility Weighted Data
  ytildes <- lapply(1:TT, FUN = function(tt){
    yy <- y[[tt]]
    g <- resp[[tt]]
    
    yy <- apply(yy, MARGIN = 2, FUN = function(x) x * g)
    colSums(yy)
    
  })  %>% bind_rows() %>% as.matrix()
  
  
  ## Auxiliary term, needed to make the objective interpretable
  aux.y <- Reduce("+", lapply(1:TT, FUN = function(tt){
    yy <- y[[tt]]
    g <- sqrt(resp[[tt]])
    
    yy <- apply(yy, MARGIN = 2, FUN = function(x) x * g)
    
    
    sum(diag(yy %*% Sigma_inv %*% t(yy)))
    
  }))
  
  ## Mu, d x T matrix
  mumat <- CVXR::Variable(rows=dimdat, cols=TT)
  
  ## Summed sqrt responsibilities - needed in the objective.
  resp.sum.sqrt <- lapply(resp, FUN = function(x) sqrt(sum(x)))
  
  
  ## Differencing Matrix, TT-l + 1 x TT 
  Dl <- gen_diff_mat(n = TT, l = l)
  
  ## Forming the objective
  
  obj = 1/(2*N) *( Reduce("+", lapply(1:TT, FUN = function(tt) quad_form(resp.sum.sqrt[[tt]]*mumat[,tt], Sigma_inv))) -2 * Reduce("+", lapply(1:TT, FUN = function(tt) t(ytildes[tt,]) %*% Sigma_inv %*% mumat[,tt])) + aux.y) + lambda * sum(sum_entries(abs(Dl %*% t(mumat)), axis = 1))
  
  ## Putting together the ball constraint
  
  rowmns <- matrix(rep(1, TT^2), ncol = TT)/TT
  
  mu_dotdot <- mumat %*% rowmns
  
  constraints = list()
  if(!is.null(maxdev)){
    constraints = list(CVXR::sum_entries(CVXR::square(mumat - mu_dotdot), axis = 2) <= rep(maxdev^2, TT) )
    
  }
  
  ## Try all two CVXR solvers.
  prob <- CVXR::Problem(CVXR::Minimize(obj), constraints)
  result = NULL
  result <- tryCatch({
    solve(prob, solver="ECOS",
          FEASTOL = ecos_thresh, RELTOL = ecos_thresh, ABSTOL = ecos_thresh)
  }, error=function(err){
    err$message = paste(err$message,
                        "\n", "Lasso solver using ECOS has failed." ,sep="")
    cat(err$message, fill=TRUE)
    return(NULL)
  })
  
  
  
  ## If anything is wrong, flag to use SCS solver
  scs = FALSE
  if(is.null(result)){
    scs = TRUE
  } else {
    if(result$status != "optimal") scs = TRUE
  }
  
  ## Use the SCS solver
  if(scs){
    result = solve(prob, solver="SCS", eps = scs_eps)
    if(any(is.na(result$getValue(mumat)))){ ## A clumsy way to check.
      stop("Lasso solver using both ECOS and SCS has failed.", sep="")
    }
  }
  
  ## Record Interesting Parameters
  num_iters <- result$num_iters
  status <- result$status
  mumat <- result$getValue(mumat)
  val <- result$value
  
  return(list(mu = mumat, value = val, status = status, num_iters = num_iters))
}


# plotting utility --------------------------------------------------------


plot_cvxr_out <- function(cvxr.mn, ylist, main = "CVXR"){
  
  ydf <- do.call(rbind, lapply(1:length(ylist), FUN = function(t){
    cbind(ylist[[t]], t)
  }))
  plot(c(0,0), xlim = c(0, length(ylist)), ylim = c(min(ydf[,-3])-1.5, max(ydf[,-3]))+0.5, type = 'n',
       ylab = "", xlab = "Time", main = main)
  grid()
  points(ydf[,3], ydf[,1], col = alpha('rosybrown1', 0.25), cex = 0.75)
  points(ydf[,3], ydf[,2], col = alpha('lightskyblue3', 0.25) , cex = 0.75)
  
  lines(x = 1:length(ylist), cvxr.mn[1,], col = 'red', lwd = 3.5)
  lines(x = 1:length(ylist), cvxr.mn[2,], col = 'blue', lwd = 3.5)
}