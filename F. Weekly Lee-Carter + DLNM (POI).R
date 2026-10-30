fit_wLC.DLNM.POI <- function(xv, isoyv, isowv, etx, dtx, cb, adjmat,
                             lambda, eps1 = 1e-9, eps2 = 1e-9, tol.ll = 1e-9,
                             max.iter = 50, verbose = TRUE, init = NULL) {
  
  ll.old    <- -Inf
  converged <- FALSE
  iter      <- 0
  
  while (!converged & iter < max.iter) {
    iter <- iter + 1
    
    # --- Step 1: Baseline model ---
    lcb <- wLC.POI(
      xv                = xv, 
      isoyv             = isoyv, 
      isowv             = isowv,
      etx               = etx, 
      dtx               = dtx,
      int               = xv * 0,
      constraints       = "ALL",
      eps               = eps1,
      adjmat            = adjmat,
      initialize        = if(iter == 1) init[[1]] else NULL,
      verbose           = FALSE,
      prev.dlnm         = if(iter == 1) NULL else lcd
    )
    
    # --- Step 2: DLNM model ---
    lcd <- wLC.DLNM.POI(
      xv                = xv, 
      isoyv             = isoyv, 
      isowv             = isowv,
      etx               = etx, 
      dtx               = dtx,
      int               = xv * 0,
      constraints       = "ALL",
      cb                = cb,
      adjmat            = adjmat,
      lambda            = lambda,
      eps               = eps2,
      prev.base         = lcb,
      initialize        = if(iter == 1) init[[2]] else NULL,
      verbose           = FALSE
    )
    
    # --- Step 3: Convergence check (ll.old must be updated LAST) ---
    d.ll      <- if (is.finite(ll.old)) abs(lcd$ll - ll.old) / abs(ll.old) else Inf
    if (lcd$ll < ll.old) warning("penalised log-likelihood decreased")
    if (verbose) cat(sprintf("cycle %3d | penalised loglik = %.4f | rel. change = %.2e\n",
                             iter, lcd$ll, d.ll))
    converged <- d.ll < tol.ll
    ll.old    <- lcd$ll
  }
  if (!converged) warning("no convergence after ", max.iter, " cycles")
  
  # --- Step 4: Return list ---
  list(lcb = lcb, lcd = lcd, itr = iter, cvg = converged)
}
