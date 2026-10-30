bootstrap_nb <- function(B, obj.unpen, obj.pen, xv, isoyv, isowv, etx, adjmat, 
                         probs = c(0.025, 0.975),
                         ncores = parallel::detectCores() - 1
) {
  
  library(parallel)
  
  # 2. Mhat 
  mhat <- obj.unpen$lcd$mhat
  dhat <- etx * mhat
  
  # 3. Dispersion
  exp.ds  <- exp(outer(obj.unpen$lcd[['phix']], obj.unpen$lcd[['phir']], "+"))
  exp.ds  <- replicate(dim(etx)[1], exp.ds, simplify = "array") 
  exp.ds  <- aperm(exp.ds, c(3, 1, 2))
  
  # 3. Bootstrap (parallel)
  cl <- snow::makeCluster(ncores)
  doSNOW::registerDoSNOW(cl)
  snow::clusterEvalQ(cl, c(library(dplyr), library(MASS), library(lubridate)))
  snow::clusterExport(cl, list = c('dhat', 'exp.ds', 'etx', 'xv',
                                   'isowv', 'isoyv', 'adjmat', 'obj.unpen',
                                   'fit_wLC.DLNM.NB', 'wLC.NB', 'wLC.DLNM.NB',
                                   'llmax.beta.dlnm', 'llmax.gamma.dlnm',
                                   'llmax.beta1', 'llmax.beta2', 'llmax.beta3',
                                   'llmax.kappa2', 'llmax.kappa3', 'llmax.phir',
                                   'llmax.phix', 'dlnm_penalty', 'obj.pen'),
                      envir = environment())

  pb <- progress::progress_bar$new(format = "Bootstrapping = :letter [:bar] :elapsed | eta: :eta",
                                   total = B, width = 60)
  progress_letter <- 1:B
  progress <- function(n){
    pb$tick(tokens = list(letter = progress_letter[n]))
  }
  opts <- list(progress = progress)
  
  boot_lst <- foreach(k = 1:B, .options.snow = opts) %dopar% {
    
    # Simulate NB
    set.seed(k)
    dtx_star <- array(
      rnbinom(length(dhat), size = exp.ds, mu = dhat),
      dim = dim(dhat)
    )
    dtx_star[which(dtx_star == 0)] <- 0.001
    dimnames(dtx_star) <- dimnames(etx)
    
    # Refit model
    fit <- fit_wLC.DLNM.NB(xv = xv, isoyv = isoyv, isowv = isowv,
                           etx = etx, dtx = dtx_star, adjmat = adjmat,
                           cb = obj.unpen$lcd$cb, 
                           lambda = obj.pen$lcd$lambda, 
                           init = obj.pen, max.iter = 50,
                           update.dispersion = FALSE)
    
    # Return
    lst         <- fit$lcd
    lst$cvg     <- fit$cvg
    lst$itr     <- fit$itr
    lst$mhat    <- NULL
    lst$epsilon <- NULL
    lst$cb      <- NULL
    
    lst
  }
  stopCluster(cl)

  # 4. Output
  boot_lst
  
}