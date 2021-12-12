
# get matrices ------------------------------------------------------------

etilde_mat <- function(TT){
  
  mats <- lapply(1:TT, FUN = function(t){
    e_vec <- rep(0, TT)
    e_vec[t] <- 1
    # browser()
    (e_vec - 1/TT) %*% t(e_vec - 1/TT)
  })
  
  Reduce('+', mats)
  
}

get_AB_mats <- function(y, resp, Sigma_inv, e_mat, N, Dl, Dlm1, rho, z, w, uz, uw){
  
  # A matrix
  A <- 1/N * Sigma_inv
  
  # B matrix
  sum_resp <- lapply(resp, sum)
  B <- rho*(t(Dlm1)%*%Dlm1 + e_mat)%*% diag(1/unlist(sum_resp))
  
  return(list(A = A, B = B))
}

get_C_mat <- function(y, resp, Sigma_inv, e_mat, N, Dl, Dlm1, rho, z, w, uz, uw){
  sum_resp <- lapply(resp, sum)
  
  
  TT <- length(y)
  
  d <- ncol(y[[1]])
  
  
  # first component
  C1 <- do.call(cbind, lapply(1:TT, FUN = function(t){
    multmat <- apply(y[[t]], FUN = function(yy) yy*resp[[t]], MARGIN = 2)
    Sigma_inv %*% colSums(multmat)
  }))
  
  
  # second component 
  C2 <- do.call(cbind, lapply(1:TT, FUN = function(t){
    uz[t,] - rho*z[t,]
  }))
  
  # averaging
  C2 <- C2 - rowMeans(C2)
  
  # third component
  
  C3 <- do.call(rbind, lapply(1:d, FUN = function(j){
    (uw[j,] - rho*w[j,]) %*% Dlm1
  }))
  
  
  # combining
  C <- (-1/N*C1 + C2 + C3) %*% diag(1/unlist(sum_resp))
  
  return(C)
}

myschur <- function(mat){
  stopifnot(nrow(mat) == ncol(mat))
  if(is.numeric(mat) & length(mat)==1) mat = mat %>% as.matrix()
  obj = Matrix::Schur(mat)
  obj$tQ = t(obj$Q)
  obj$orig = mat
  return(obj)
}


# convergence/augmented lagrangian/objective ------------------------------


converge <- function(mu, rho, w, z, w_prev, z_prev, uw, uz, Dlm1,
                     err_rel = 1E-3,
                     err_abs = 0
){
  
  
  prim1 = rbind(t(mu - rowMeans(mu)), Dlm1 %*% t(mu))
  prim2 = rbind(z, t(w))
  primal_resid = prim1 - prim2
  
  
  dual_resid = -rho * rbind((z-z_prev - colMeans(z-z_prev)), t(Dlm1) %*% t(w - w_prev))
  tAU = rbind(uz, -t(uw)) 
  
  ## Form primal and dual tolerances.
  primal_err = sqrt(length(primal_resid)) * err_abs +
    err_rel * max(norm(prim1, "F"), norm(prim2, "F"))
  dual_err = sqrt(length(dual_resid)) * err_abs +
    err_rel * norm(tAU, "F")
  
  ## Check convergence.
  primal_resid_size = norm(primal_resid, "F")
  dual_resid_size = norm(dual_resid, "F")
  primal_converge = ( primal_resid_size  <= primal_err )
  dual_converge = ( dual_resid_size <= dual_err )
  
  ## Some checks (trying to address problems with |converge|).
  assertthat::assert_that(is.numeric(primal_resid_size))
  assertthat::assert_that(is.numeric(primal_err))
  assertthat::assert_that(is.numeric(dual_resid_size))
  assertthat::assert_that(is.numeric(dual_err))
  
  ## return(primal_converge & dual_converge)
  converge = primal_converge & dual_converge
  
  return(list(primal_resid = primal_resid,
              primal_err = primal_err,
              dual_resid = dual_resid,
              dual_err = dual_err,
              converge = converge))
}


aug_lagr <- function(y, TT, d, z, w, l, uz, uw, mu, resp, Sigma_inv, Dl, Dlm1, maxdev, lambda, rho, N){
  mu_dd <- rowMeans(mu)
  
  # Check the Z's for ball constraint, up to a tolerance of 1e-4
  znorms <- apply(z, FUN = function(zz) sqrt(sum(zz^2)), MARGIN = 1)
  if(any(is.na(znorms))) browser()
  if(any(znorms > (maxdev + 1e-4))){
    warning("||z||_2 > maxdev, current iterate not feasible.")
    return(Inf)
  }
  
  aug1 <- sum(do.call(cbind, lapply(1:TT, FUN = function(t){
    multmat <- apply(y[[t]], FUN = function(yy){
      t(yy - mu[,t]) %*% Sigma_inv %*% (yy - mu[,t])}, MARGIN = 1)
    sum(1/(2*N) * resp[[t]]*multmat)
  })))
  
  aug2 <- lambda*sum(do.call(rbind, lapply(1:d, FUN = function(j)
    sum(abs(diff(w[j,], differences = 1))))))
  
  aug3 <- sum(do.call(cbind, lapply(1:TT, FUN = function(t){
    uz[t,]%*%(mu[,t] - mu_dd - z[t,]) + rho/2 * sum((mu[,t] - mu_dd - z[t,])^2)
  })))
  
  
  aug4 <- sum(do.call(rbind, lapply(1:d, FUN = function(j){
    uw[j,] %*% (Dlm1 %*% mu[j,]  - w[j,]) + rho/2 * sum((Dlm1 %*% mu[j,] - w[j,])^2)
  })))
  
  total <- aug1 + aug2 + aug3 + aug4
  
  
  return(total)
  
}

obj <- function(y, TT, d, l, mu, resp, Sigma_inv, Dl, Dlm1, maxdev, lambda, rho, N){
  mu_dd <- rowMeans(mu)
  
  aug1 <- sum(do.call(cbind, lapply(1:TT, FUN = function(t){
    multmat <- apply(y[[t]], FUN = function(yy){
      t(yy - mu[,t]) %*% Sigma_inv %*% (yy - mu[,t])}, MARGIN = 1)
    sum(1/(2*N) * resp[[t]]*multmat)
  })))
  
  aug2 <- lambda*sum(do.call(rbind, lapply(1:d, FUN = function(j)
    sum(abs(diff(mu[j,], differences = l))))))
  
  total <- aug1 + aug2 
  
  return(total)
}


# W, Z, U updates ---------------------------------------------------------


W_update_fused <- function(Dlm1, TT, mu, uw, rho, lambda){
  
  # modified lambda for fused lasso routine
  mod_lam <- lambda/rho
  #print(mod_lam)
  
  # generating pseudo response xi
  xi <- Dlm1 %*% mu + 1/rho * uw
  
  #browser()
 # print(c(xi = norm(xi, "F")))
  
  # running the fused LASSO
  fit <- fused_lasso(z = xi, lam = mod_lam)
  
  return(fit)
}

Z_update  <- function(m, Uz, C, rho){
  mat = m + Uz/rho
  Z = projCmat(mat, C)
  return(Z)
}

projCmat <- function(mat, C){
  if(!is.null(C)){
    vlens = sqrt(rowSums(mat * mat))
    inds = which(vlens > C)
    if(length(inds) > 0){
      mat[inds,] = mat[inds,] * C / vlens[inds]
    }
  }
  return(mat)
}

U_update_Z <- function(U, rho, mu, Z){
  return(U + rho * (scale(mu, scale = F) - Z))
}

U_update_W <- function(U, rho, mu, W, Dlm1){
  return(U +  rho * ( t(Dlm1%*%mu) - W))
}


# main function -----------------------------------------------------------


admm_oneclust <- function(iclust = 1, niter, y,
                          Dlm1, Dl, l = NULL,
                          TT, N, dimdat, maxdev,
                          rho,
                          rhoinit = rho,
                          Xinv,
                          schurA,
                          schurB,
                          sigmainv,
                          lambda,
                          resp,
                          ylist, err_rel, err_abs,
                          zerothresh,
                          z,
                          w,
                          uw,
                          uz,
                          first_iter,## Not used
                          outer_iter,
                          local_adapt,
                          sigma,
                          sigma_eig_by_clust,
                          space = 20, norms = F, auglag = F, objective  = F){
  
  ## Initialize the variables ###
  resid_mat = matrix(NA, nrow = ceiling(niter/5), ncol = 4)
  colnames(resid_mat) = c("primresid", "primerr", "dualresid", "dualerr")
  
  rhofac = rho / rhoinit
  
  ## Main inner LA-ADMM loop
  fits = rep(NA, ceiling(niter/space))
  converge = FALSE
  start.time = Sys.time()
  Zlist = list()
  
  ## This doesn't change over iterations
  schurA = myschur(schurA$orig * rhofac)
  TA = schurA$T ##* rhofac
  TB = schurB$T
  UA = schurA$Q
  UB = schurB$Q
  tUA = schurA$tQ
  tUB = schurB$tQ
  
  
  auglags <- NULL
  
  if(auglag){
    auglags <- data.frame(iter = NULL, init_aug = NULL, mu_update = NULL, z_update = NULL, w_update = NULL, final_aug = NULL)
    l <- sum(Dlm1[1,] != 0)
    init_aug = Inf
  }
  
  for(iter in 1:niter){
    
    if(auglag & iter > 1){
      init_aug <- aug_lagr(y = y, mu = mu, resp = resp, Sigma_inv = Sigma_inv, TT = TT, d = d, z = z, w = w, uw = uw, uz = uz, Dl = Dl, Dlm1 = Dlm1, l = l, maxdev = maxdev, lambda = lambda, rho = rho, N = N)
    }
    
    syl_C <- get_C_mat(y = y, resp = resp, Sigma_inv = sigmainv, N = N, Dlm1 = Dlm1, rho = rho, z = z, w = w, uz = uz, uw = uw)
    FF =  (-1) * tUA %*% syl_C %*% UB
    mu = UA %*% matrix_function_solve_triangular_sylvester_barebones(TA, TB, FF) %*% tUB
    
    if(auglag){
      mu_aug <- aug_lagr(y = y, mu = mu, resp = resp, Sigma_inv = Sigma_inv, TT = TT, d = d, z = z, w = w, uw = uw, uz = uz, Dl = Dl, Dlm1 = Dlm1, l = l, maxdev = maxdev, lambda = lambda, rho = rho, N = N)
      if(mu_aug > init_aug) {cat(paste("Mu update at iteration:", iter, "increased AUGL, difference of ", mu_aug - init_aug, "\n"))
        # browser()
      }
      
    }
    #  browser()
    z <- Z_update( scale(t(mu), scale = F), Uz = uz, C = maxdev, rho = rho)
    
    
    if(auglag){
      z_aug <- aug_lagr(y = y, mu = mu, resp = resp, Sigma_inv = Sigma_inv, TT = TT, d = d, z = z, w = w, uw = uw, uz = uz, Dl = Dl, Dlm1 = Dlm1, l = l, maxdev = maxdev, lambda = lambda, rho = rho, N = N)
      
      if(z_aug > mu_aug){cat(paste("Z update at iteration:", iter, "increased AUGL, difference of ", z_aug - mu_aug, "\n"))
        #browser()
      }
    }
    
    w <- do.call(rbind, lapply(1:dimdat, 
                               FUN = function(j) W_update_fused(Dlm1 = Dlm1, TT = TT, mu = mu[j,], rho = rho, lambda = lambda, uw = uw[j,])))
    
    if(auglag){
      w_aug <- aug_lagr(y = y, mu = mu, resp = resp, Sigma_inv = Sigma_inv, TT = TT, d = d, z = z, w = w, uw = uw, uz = uz, Dl = Dl, Dlm1 = Dlm1, l = l, maxdev = maxdev, lambda = lambda, rho = rho, N = N)
      final_aug = w_aug
      if(w_aug > z_aug) {cat(paste("W update at iteration:", iter, "increased AUGL, difference of ", w_aug - z_aug, "\n"))
      }
      
    }
    
    uz = U_update_Z(uz, rho, t(mu), z)
    uw = U_update_W(uw, rho, t(mu), w, Dlm1)
    
    
    if(norms & iter %% 1 == 0){
      print(round(c("mu" = norm(mu, "F"), "C" = norm(syl_C, "F"), "FF" = norm(FF, "F"),
                    "z" = norm(z, "F"), "w" = norm(w, "F"), "uz" = norm(uz, "F"), "uw" = norm(uw, "F")),3))
    }
    
    ## Check convergence
    if( iter > 1  & iter %% 5 == 0){## & !local_adapt){
      
      ## Calculate convergence criterion
      obj = converge(mu, rho,
                     w, z,
                     w_prev, z_prev,
                     uw, uz,
                     Dlm1)
      
      jj = (iter/ 5)
      resid_mat[jj,] = c(norm(obj$primal_resid, "F"),
                         obj$primal_err,
                         norm(obj$dual_resid,"F"),
                         obj$dual_err)


      if(auglag){
        
        auglags <- rbind(auglags, 
                         data.frame(iter = iter, init_aug = init_aug, mu_update = mu_aug, z_update = z_aug, w_update = w_aug, final_aug = final_aug))
      }
      if(is.na(obj$converge)){
        obj$converge <- FALSE
        warning("Convergence was NA")
      }
      if(obj$converge){
        converge = TRUE
        break
      }
    }
    
    ## ## 3. Calculate objective values for this cluster.
    w_prev = w
    z_prev = z
    
  }
  
  
  ## Gather results.
  mu = mu
  
  yhat = mu
  
  fit <- NULL
  fits <- NULL
  obj.value <- NULL 
  
  if(objective){
    obj.value <- obj(y = y, mu = mu, resp = resp, Sigma_inv = sigmainv, TT = TT, d = dimdat, Dl = Dl, Dlm1 = Dlm1, l = l, maxdev = maxdev, lambda = lambda, rho = rho, N = N)
  }
  
  return(list(mu = mu,
              yhat = yhat,
              resid_mat = resid_mat,
              fits = fits,
              converge = converge,
              fit = obj.value,
              ## Other variables to return.
              Z = z,
              W = w,
              uz = uz,
              uw = uw,
              auglags = auglags,
              objective = obj.value
  ))
}

# plotting utility -----------------------------------------------------------


plot_admm_out <- function(admm, ylist, main = ""){
  
  ydf <- do.call(rbind, lapply(1:length(ylist), FUN = function(t){
    cbind(ylist[[t]], t)
  }))
  plot(c(0,0), xlim = c(0, ncol(admm$yhat)), ylim = c(min(ydf[,-3])-1.5, max(ydf[,-3]))+0.5, type = 'n',
       ylab = "", xlab = "Time", main = main)
  grid()
  points(ydf[,3], ydf[,1], col = alpha('rosybrown1', 0.25), cex = 0.75)
  points(ydf[,3], ydf[,2], col = alpha('lightskyblue3', 0.25) , cex = 0.75)
  
  lines(x = 1:length(ylist), admm$mu[1,], col = 'red', lwd = 3.5)
  lines(x = 1:length(ylist), admm$mu[2,], col = 'blue', lwd = 3.5)
}


