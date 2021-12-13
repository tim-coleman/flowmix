source("C:/Users/drain/Box Sync/USC_Stuff/Featureless_Flowmix/flowmix_CVXR/Cpp_loading.R")
source("C:/Users/drain/Box Sync/USC_Stuff/Featureless_Flowmix/flowmix_CVXR/FL_loading_C.R")
source("C:/Users/drain/Box Sync/USC_Stuff/Featureless_Flowmix/flowmix_CVXR/general_utils.R")
source("C:/Users/drain/Box Sync/USC_Stuff/Featureless_Flowmix/flowmix_CVXR/admm.R")
source("C:/Users/drain/Box Sync/USC_Stuff/Featureless_Flowmix/flowmix_CVXR/CVXR.R")


#' Mstep_admm_tf
#'
#' @param resp - responsbilities of each particle.
#' @param ylist - 
#' @param lambda 
#' @param l 
#' @param sigma 
#' @param sigma_eig_by_clust 
#' @param Dlm1 
#' @param Dl 
#' @param TT 
#' @param N 
#' @param dimdat 
#' @param first_iter 
#' @param mus 
#' @param Zs 
#' @param Ws 
#' @param uws 
#' @param uzs 
#' @param maxdev 
#' @param niter 
#' @param err_rel 
#' @param err_abs 
#' @param zerothresh 
#' @param local_adapt 
#' @param local_adapt_niter 
#' @param space 
#'
#' @return
#' @export
#'
#' @examples
Mstep_admm_tf <- function(resp,
                            ylist,
                            lambda = 0.5,
                            l = 3,
                            sigma,
                            sigma_eig_by_clust = NULL,
                            Dlm1, Dl,  TT, N, dimdat,
                            first_iter = TRUE,
                            
                            ## Warm startable variables
                            mus = NULL,
                            Zs = NULL,
                            Ws = NULL,
                            uws = NULL,
                            uzs = NULL,
                            ## End of warm startable variables
                            
                            maxdev = NULL,
                            niter = (if(local_adapt) 1e3 else 1e4),
                            err_rel = 1E-3,
                            err_abs = 0,
                            zerothresh = 1E-6,
                            local_adapt = FALSE,
                            local_adapt_niter = 10,
                            space = 50
){
  
  ####################
  ## Preliminaries ###
  ####################
  TT = length(ylist)
  numclust = ncol(resp[[1]])
  dimdat = ncol(ylist[[1]])
  ntlist = sapply(ylist, nrow)
  #resp.sum = t(sapply(resp, colSums))
 # resp.sum = t(sapply(resp, colSums)) ## (T x numclust)
#  resp.sum = as.matrix(resp.sum)
  resp.sum = lapply(resp, colSums) %>% do.call(rbind, .)
  #resp.sum <- lapply(resp.sum, FUN = function(r) matrix(r, ncol = numclust))
  N = sum(unlist(resp.sum)) ## NEW (make more efficient, later)
  
  # starting rho for LA-ADMM
  rho.init = 1e-1


  ## Other preliminaries
  schur_syl_A_by_clust = schur_syl_B_by_clust = term3list = list()
  ybarlist = list()
  ycentered_list = Xcentered_list = yXcentered_list = list()
  Qlist = list()
  sigmainv_list = list()
  convergences = list()
  for(iclust in 1:numclust){
    ## print(iclust, numclust, "iclust")
    ## Retrieve sigma inverse from pre-computed SVD, if necessary
    if(is.null(sigma_eig_by_clust)){
      sigmainv = solve(sigma[iclust,,])
    } else {
      sigmainv = sigma_eig_by_clust[[iclust]]$sigma_inv
    }
    
    
    resp.iclust <- lapply(resp, FUN = function(r) matrix(r[,iclust]))
    
    ## Center y and X
    obj <- weight_ylist(iclust, resp, resp.sum, ylist)
    ycentered <- obj$ycentered

    ## Form the Sylvester equation coefficients in AX + XB + C = 0
    # syl_A = rho * sigma[iclust,,]
    # Q = 1/N * t(Xcentered) %*% D %*% Xcentered
    # syl_B = Q %*% Xinv

    e_mat <- etilde_mat(TT = TT) # needed to generate B
    AB <- get_AB_mats(y = y, resp = resp.iclust, Sigma_inv = sigmainv, e_mat = e_mat, N = TT*nt, Dl = Dl,
                      Dlm1 = Dlm1, rho = (if(local_adapt) rho.init else lambda), z = NULL, w = NULL, uz = NULL, uw = NULL)
    
   # browser()
    ## Store the Schur decomposition
    schur_syl_A_by_clust[[iclust]] = myschur(AB$A)
    schur_syl_B_by_clust[[iclust]] = myschur(AB$B)
    

    ## Calculate coefficients for objective value  calculation
    # Qlist[[iclust]] = Q
    
    ## ## Also calculate some things for the objective value
    ## ylong = sweep(do.call(rbind, ylist), 2, obj$ybar)
    ## longwt = do.call(c, lapply(1:TT, function(tt){ resp[[tt]][,iclust]})) %>% sqrt()
    ## wt.long = longwt * ylong
    ## wt.ylong = longwt * ylong
    ## crossprod(wt.ylong, wt.ylong)
    
    ## Store the third term
    # term3list[[iclust]] =  1 / N * sigmainv %*% yXcentered
    # ybarlist[[iclust]] = obj$ybar
    # 
     ycentered_list[[iclust]] = ycentered
     #print(ycentered)
    # Xcentered_list[[iclust]] = Xcentered
    # yXcentered_list[[iclust]] = yXcentered
    sigmainv_list[[iclust]] = sigmainv
  }
  
  ##########################################
  ## Run ADMM separately on each cluster ##
  #########################################
  yhats = admm_niters = admm_inner_iters = vector(length = numclust, mode = "list")
  if(first_iter) mus = vector(length = numclust, mode = "list")
  if(first_iter){
    Zs <- lapply(1:numclust, function(x) matrix(0, nrow = TT, ncol = dimdat) )
    Ws <- lapply(1:numclust, function(x) matrix(0, nrow = dimdat, ncol = TT - l + 1))
    uzs <- lapply(1:numclust, function(x) matrix(0, nrow = TT, ncol = dimdat) )
    uws <- lapply(1:numclust, function(x) matrix(0, nrow = dimdat, ncol = TT - l + 1))
    
   # Zs =  Ws =  Us  = vector(length = numclust, mode = "list")
  }
  
  fits = matrix(NA, ncol = numclust, nrow = ceiling(niter / space))
  
  #browser()
  
  ## For every cluster, run LA-ADMM
  resid_mat_list = list()
  start.time = Sys.time()
  for(iclust in 1:numclust){
    resp.iclust <- lapply(resp, FUN = function(r) matrix(r[,iclust]))
    resp.sum.iclust <- lapply(resp.sum, FUN = function(r) matrix(r[iclust]))
    
    ## Possibly locally adaptive ADMM, for now just running with rho == lambda
    res = la_admm_oneclust(K = (if(local_adapt) local_adapt_niter else 1),
                           local_adapt = local_adapt,
                           iclust = iclust,
                           niter = niter,
                           TT = TT, N = N, dimdat = dimdat, maxdev = maxdev,
                           schurA = schur_syl_A_by_clust[[iclust]],
                           schurB = schur_syl_B_by_clust[[iclust]],
                           #term3 = term3list[[iclust]],
                           sigmainv = sigmainv_list[[iclust]],
                           # Xinv = Xinv,
                           # Xaug = Xaug,
                           # Xa = Xa,
                           rho =  (if(local_adapt) rho.init else lambda),
                           rhoinit = (if(local_adapt) rho.init else lambda),
                           sigma = sigma,
                           lambda = lambda,
                           resp = resp.iclust,
                           l = l,
                           Dl = Dl,
                           Dlm1 = Dlm1,
                          #resp.sum = resp.sum.iclust,
                           y = ylist, 
                           err_rel = err_rel,
                           err_abs = err_abs,
                           zerothresh = zerothresh,
                           sigma_eig_by_clust = sigma_eig_by_clust,
                           space = space,
                           objective = T, norms = F,
                           
                           ## Warm starts from previous *EM* iteration
                           first_iter = first_iter,
                          # beta = betas[[iclust]],
                           #mu = mus[[iclust]],
                           uw = uws[[iclust]],
                           uz = uzs[[iclust]],
                           z = Zs[[iclust]],
                           w = Ws[[iclust]]
    )
    
    ## Store the results
    mus[[iclust]] = res$mu
    yhats[[iclust]] = res$yhat
    ## fits[,iclust] = res$fits
    admm_niters[[iclust]] = res$kk
    admm_inner_iters[[iclust]] = res$inner.iter
    
    ## Store other things for for warmstart
    Zs[[iclust]] = res$Z
    uzs[[iclust]] = res$uz
    uws[[iclust]] = res$uw
    Ws[[iclust]] = res$W
    
    ## The upper triangular matrix remains the same.
    ## The upper triangular matrix remains the same.
    
    resid_mat_list[[iclust]] = res$resid_mat ## temporary
    convergences[[iclust]] = res$converge
  }
  
  ## Aggregate the yhats into one array
  yhats_array = array(NA, dim = c(TT, dimdat, numclust))
  for(iclust in 1:numclust){ yhats_array[,,iclust] = yhats[[iclust]] }
  
  ## Each are lists of length |numclust|.
  return(list(mus = mus,
              mns = yhats_array,
              fits = fits,
              resid_mat_list = resid_mat_list,
              convergences = convergences,
              admm_niters = admm_niters, ## Temporary: Seeing the number of
              ## outer iterations it took to
              ## converge.
              admm_inner_iters = admm_inner_iters,
              
              ## For warmstarts
              Zs = Zs,
              Ws = Ws,
              uws = uws,
              uzs = uzs,
              
              ## For using in the Sigma M step
              ycentered_list = ycentered_list,
              Xcentered_list = Xcentered_list,
              yXcentered_list = yXcentered_list,
              Qlist = Qlist
  ))
}



# Examples ----------------------------------------------------------------

library(flowmix)

# data and utilities
set.seed(0)
numclust = 4
TT = 50
nt = 100
datobj = generate_data_generic(p=5, TT=50, fac=.5, nt=nt, dimdat = 2)
ylist = datobj$ylist

# model parameters
lambda = 0.2
l = 3 
maxdev = 5
Dl = gen_diff_mat(n = TT, l = l)
Dlm1 = gen_diff_mat(n = TT, l = l-1)

# coercing the sigma's to an array
sigma.array <- array(NA, dim = c(numclust, 2, 2))
for(ii in 1:numclust){
  sigma.array[ii,,] <- datobj$sigmalist[[ii]]
}


# generating responsibilities - taking simple proportions (plus a small number to get rid of nonzero responsbiliites.)
resps <- lapply(datobj$classlist, FUN = function(classes){
  
  classes.fac <- factor(classes, levels = 1:numclust)
  
  probs.class <- (table(classes.fac) + 0.001)/(length(classes.fac)*(1 + numclust*0.001))
  
  class.df <- matrix(probs.class, nrow = numclust, ncol = length(classes)) %>% t()
  
  return(class.df)
})


m.init <- Mstep_admm_tf(resp = resps, ylist = ylist, sigma = sigma.array, local_adapt = F, maxdev = maxdev, local_adapt_niter = 5,
                        Dlm1 = Dlm1, Dl = Dl, TT = TT, N = TT*nt, lambda = lambda, l = l, niter = 1e5, first_iter = T)


library(profvis)
profvis(m.init <- Mstep_admm_tf(resp = resps, ylist = ylist, sigma = sigma.array, local_adapt = F, maxdev = maxdev, local_adapt_niter = 5,
                                Dlm1 = Dlm1, Dl = Dl, TT = TT, N = TT*nt, lambda = lambda, l = l, niter = 1e5, first_iter = T))
