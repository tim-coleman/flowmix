flowmix_tf_once <- function(ylist,
                         countslist = NULL,
                         numclust, niter = 1000, l,
                         mn = NULL, lambda, verbose = FALSE,
                         sigma_fac = 1, tol_em = 1E-4,
                         maxdev = NULL,
                         countslist_overwrite = NULL,
                         zero_stabilize  = FALSE,
                         init_mn_flatten = FALSE,
                         ## beta Mstep (CVXR) settings
                         mstep_cvxr_ecos_thresh = 1E-8,
                         mstep_cvxr_scs_eps = 1E-5,
                         zerothresh = 1E-6,
                         ## beta Mstep (ADMM) settings
                         admm = TRUE,
                         admm_rho = 0.01,
                         admm_err_rel = 1E-3,
                         admm_err_abs = 1E-4,
                         ## beta M step (Locally Adaptive ADMM) settings
                         admm_local_adapt = TRUE,
                         admm_local_adapt_niter = 10,
                         admm_niter = (if(admm_local_adapt)1E3 else 1E4),
                         CVXR =FALSE, ## temporary
                         plot_vals = F
){
  
  ## Basic checks
  if(!is.null(maxdev)){
    assertthat::assert_that(maxdev!=0)
  } else {
    maxdev = 1E10 ## Some large number
  }
  ## assert_that(!(is.data.frame(ylist[[1]])))
  #assertthat::assert_that(!(is.data.frame(X)))
  #assertthat::assert_that(sum(is.na(X)) == 0)
  #assertthat::assert_that(length(ylist) == nrow(X))
  ## assertthat::assert_that(prob_lambda > 0)
  assertthat::assert_that(numclust > 1)
  assertthat::assert_that(niter > 1)
  
  ## assertthat::assert_that(all(sapply(ylist, nrow) == sapply(countslist, length)))
  
  ## Setup
  TT = length(ylist)
  dimdat = ncol(ylist[[1]])
  #p = ncol(X)
  
  if(is.null(mn)) mn = init_mn(ylist, numclust, TT, dimdat)
  ntlist = sapply(ylist, nrow)
  N = sum(ntlist)
  
  ## Initialize some objects
  prob = matrix(1/numclust, nrow = TT, ncol = numclust) ## Initialize to all 1/K.
  denslist_by_clust <- NULL
  objectives = c(+1E20, rep(NA, niter-1))
  sigma = init_sigma(ylist, numclust, sigma_fac) ## (T x numclust x dimdat x dimdat)
  sigma_eig_by_clust = NULL
  zero.betas = zero.alphas = list()
  admm_niters = list()
  
  ## Warm startable variables
  mus = NULL
  Zs = NULL
  Ws = NULL
  uws = NULL
  uzs = NULL
  


  ## The least elegant solution I can think of.. used only for blocked cv
  if(!is.null(countslist_overwrite)) countslist = countslist_overwrite
  #if(!is.null(countslist)) check_trim(ylist, countslist)
  
  
  vals <- vector(length = niter)
  start.time = Sys.time()
  for(iter in 2:niter){
    if(verbose){
      print_progress(iter-1, niter-1, "EM iterations.", start.time = start.time)
    }
    resp <- Estep_nocpp(mn, sigma, prob, ylist = ylist, numclust = numclust,
                  denslist_by_clust = denslist_by_clust,
                  first_iter = (iter == 2), countslist = countslist)
    ## M step (three parts)
    ## 1. Alpha
    # res.alpha = Mstep_alpha(resp, X, numclust, lambda = prob_lambda,
    #                         zerothresh = zerothresh)
    # prob = res.alpha$prob
    # alpha = res.alpha$alpha
    # rm(res.alpha)
    
    ## 2. Beta
    

    ## temporary
    if(CVXR){
      res.mu = Mstep_beta(resp, ylist, l = l,
                            lambda = lambda,
                            first_iter = (iter == 2),
                            sigma_eig_by_clust = sigma_eig_by_clust,
                            sigma = sigma, maxdev = maxdev)}
    } else {
      res.mu = Mstep_admm_tf(resp, ylist,
                                 lambda = lambda,
                                 first_iter = (iter == 2),
                                 l = l, 
                                 sigma_eig_by_clust = sigma_eig_by_clust,
                                 sigma = sigma, maxdev = maxdev, 
                                 betas = betas,
                                 Zs = Zs,
                                 Ws = Ws,
                                 uws = NULL,
                                 uzs =  NULL,
                                 err_rel = admm_err_rel,
                                 err_abs = admm_err_abs,
                                 local_adapt = admm_local_adapt,
                                 local_adapt_niter = admm_local_adapt_niter)
    }
    
   # admm_niters[[iter]] = unlist(res.beta$admm_niters)
    
    ## Harvest means
    mn = res.mu$mns
    betas = beta = res.mu$beta
    
    ## Harvest other things for next iteration's ADMM.
    # Zs = res.beta$Zs
    # Ws = res.beta$Ws
    # Us = res.beta$Us
    # ## rm(res.beta)
    
    ## Check if the number of zeros in the alphas and betas have stabilized.
    zero.betas[[iter]] = lapply(beta, function(mybeta) which(mybeta==0))
    #zero.alphas[[iter]] = which(alpha == 0)
    if(zero_stabilize & iter >= 30){ ## If 5 is to low, try 10 instead of 5.
      if(check_zero_stabilize(zero.betas, iter)) break
    }
    
    ## 3. Sigma
    sigma = Mstep_sigma(resp, ylist, mn, numclust)
    
    ## 3. (Continue) Decompose the sigmas.
    sigma_eig_by_clust <- eigendecomp_sigma_array(sigma)
    denslist_by_clust <- make_denslist_eigen(ylist, mn, TT, dimdat, numclust,
                                             sigma_eig_by_clust,
                                             countslist)
    
    
    ## 4. Probabilities
    prob = do.call(rbind, mstep_pi(resp, countslist = countslist))

    
    ## Calculate the objectives
    objectives[iter] = objective(mn, prob, sigma, ylist,
                                 lambda = lambda, l = l,
                                 denslist_by_clust = NULL,
                                 countslist = countslist,
                                 sep = F,
                                 each = F)
    
    vals[iter] <- res.beta$value

    ## Check convergence
    ## if(iter > 10){ ## don't stop super early. ## We might not need this.
    # if(check_converge_rel(objectives[iter-1],
    #                       objectives[iter],
    #                       tol = tol_em)) break
    # ## }
    ## if(objectives[iter] > objectives[iter-1] * 1.01 ) break # Additional stopping
    ## of the likelihood
    ## increasing more
    ## than 1%.
  }
  
  ## Measure time
  lapsetime = difftime(Sys.time(), start.time, units = "secs")
  time_per_iter = lapsetime / (iter-1)
  
  ## Also calculate per-cytogram likelihoods (NOT divided by nt)
  # loglikelihoods = objective(mn, prob, sigma, ylist,
  #                            prob_lambda = prob_lambda,
  #                            mean_lambda = mean_lambda,
  #                            alpha = alpha, beta = beta,
  #                            denslist_by_clust = denslist_by_clust,
  #                            countslist = countslist,
                           #  each = TRUE)
  ## loglikelihoods_particle = objective(mn, prob, sigma, ylist,
  ##                            prob_lambda = prob_lambda,
  ##                            mean_lambda = mean_lambda,
  ##                            alpha = alpha, beta = beta,
  ##                            denslist_by_clust = denslist_by_clust,
  ##                            countslist = countslist,
  ##                            each = FALSE,
  ##                            sep=TRUE)
  
  ## Also reformat the coefficients
  #obj <- reformat_coef(alpha, beta, numclust, dimdat, X)
 # alpha = obj$alpha
 # beta = obj$beta
  
  if(plot_vals){
  plot(vals, type = 'o', main = "CVX Objective", xlab = "iteration")
  }
  return(structure(list(beta = beta,
                        mn = mn,
                        prob = prob,
                        sigma = sigma,
                        ## denslist_by_clust = denslist_by_clust,
                        objectives = objectives[2:iter],
                        final.iter = iter,
                        time_per_iter = time_per_iter,
                        total_time = lapsetime,
                      #  loglikelihoods = loglikelihoods,
                        ## loglikelihoods_particle = loglikelihoods_particle,
                        ## Above is output, below are data/algorithm settings.
                        dimdat = dimdat,
                        TT = TT,
                        N = N,
                        l = l,
                        numclust = numclust,
                        lambda = lambda,
                        maxdev=maxdev,
                        niter = niter,
                        admm_niters = admm_niters
  ), class = "flowmix"))
}

## ## Some tests to add
## ## object is the result of having run flowmix() or flowmix_once().
## check_size <- function(obj){
##   assert_that(check_beta_size(res$beta, p, dimdat, numclust))
##   assert_that(check_alpha_size(res$alpha, p, dimdat))
## }

## check_beta_size <- function(beta, p, dimdat, numclust){
##   all.equal(dim(beta), c(p+1, dimdat, numclust))
## }
## check_alpha_size <- function(alpha, p, dimdat){
##   all.equal(dim(alpha), c(dimdat, p+1))
## }



##' Prediction: Given new X's, generate a set of means and probs (and return the
##' same Sigma).
##'
##' @param res Object returned from covariate EM flowmix().
##' @param newx New covariate.
##'
##' @return List containing mean, prob, and sigma.
##'
##' @export
##'
predict.flowmix <- function(res, newx = NULL){
  
  ## ## Check the dimensions
  ## stopifnot(ncol(new.x) == ncol(res$X))
  ## newx = X[1,,drop=FALSE]
  if(is.null(newx)){
    newx = res$X
  }
  
  ## Check if the variable names are the same.
  cnames = res$X %>% colnames()
  cnames_new = newx %>% colnames()
  stopifnot(all(cnames == cnames_new))
  
  ## Augment it with a dummy variable 1
  if(nrow(newx)>1){
    newx.a = cbind(rep(1, nrow(newx)), newx)
  } else {
    newx.a = c(1, newx)
  }
  
  TT = nrow(newx) ## This used to be nrow(X)..
  numclust = res$numclust
  dimdat = res$dimdat
  if(is.null(dimdat)) dimdat = res %>%.$mn %>% dim() %>% .[2] ## for back=compatibility
  
  ## Predict the means (manually).
  newmn = lapply(1:numclust, function(iclust){
    newx.a %*% res$beta[[iclust]]
  })
  newmn_array = array(NA, dim=c(TT, dimdat, numclust))
  for(iclust in 1:numclust){ newmn_array[,,iclust] = newmn[[iclust]] }
  
  ## Predict the probs.
  ## newprob = predict(res$alpha.fit, newx=newx, type='response')[,,1]
  
  probhatmat = as.matrix(exp(cbind(1,newx) %*% t(res$alpha)))
  newprob = probhatmat / rowSums(probhatmat)
  ## predict(fit, newx=X, type="response")[,,1]
  stopifnot(all(dim(newprob) == c(TT,numclust)))
  stopifnot(all(newprob >= 0))
  
  ## Return all three things
  return(list(mn = newmn_array,
              prob = newprob,
              pie = newprob, ## Just a copy of prob, for back-compatibility
              alpha = res$alpha,
              beta = res$beta,
              sigma = res$sigma,
              TT = res$TT,
              N = res$N,
              numclust = res$numclust,
              X = newx))
}



##' Helper for making list of densities. Returns list by cluster then time
##' e.g. access by \code{denslist_by_clust[[iclust]][[tt]]}
##'
##' @param ylist T-length list each containing response matrices of size (nt x
##'   3), which contains coordinates of the 3-variate particles, organized over
##'   time (T) and with (nt) particles at every time.
##' @param mu (T x dimdat x numclust) array.
##' @param dimdat dimension of data.
##' @param numclust number of clusters.
##' @param TT number of time points
##' @param sigma_eig_by_clust Result of running
##'   \code{eigendecomp_sigma_array(sigma.list[[iter]])}.
##'
##' @return numclust-lengthed list of TT-lengthed.
##'
make_denslist_eigen <- function(ylist, mu,
                                TT, dimdat, numclust,
                                sigma_eig_by_clust,
                                countslist){ ## Temporary
  
  ## Basic checks
  assertthat::assert_that(!is.null(sigma_eig_by_clust))
  
  ## Calculate densities (note to self: nested for loop poses no problems)
  lapply(1:numclust, function(iclust){
    mysigma_eig <- sigma_eig_by_clust[[iclust]]
    lapply(1:TT, function(tt){
      ## return(dmvnorm_fast(ylist[[tt]],
      ##                     mu[tt,,iclust],
      ##                     sigma_eig=mysigma_eig))
      mn = mu[tt,,iclust]
      sgm = mysigma_eig$sigma
      if(dimdat == 1){
        mn = as.matrix(mn)
        sgm = sgm %>% as.matrix()
      }
      return(mvtnorm::dmvnorm(ylist[[tt]],
                               mn,
                               sgm))
    })
  })
}

