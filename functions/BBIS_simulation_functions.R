# simulate langevin track (High-res)
sumLTrack <- function(dt, gamma2, covlist, beta, loc0, nobs) {
  x = mvnfast::rmvn(nobs, rep(0,2), dt * gamma2 * diag(1,2,2))
  x[1, ] <- loc0 
  for (i in 2:nrow(x)) {
    grad = bilinearGradVec(matrix(x[i-1, 1:2], nrow=1), covlist)
    x[i, ]  = x[i, ] + x[i-1, ] + (dt * gamma2 / 2)*beta %*% grad[,1,]
  }
  return(x)
}

# thin high-resolution track
thinTrack <- function(x, thin) {
  n = nrow(x)
  x = x[seq(1, n, thin), ]
  return(x)
}

# simulate thinned Langevin movement model
simLMM <- function(dt, gamma2, covlist, beta, loc0, nobs){ 
  x <- sumLTrack(dt, gamma2, covlist, beta, loc0, nobs)  # Simulate high-resolution track
  #x <- thinTrack(x, delta) # Thin the track
  return(x)
}

# define cpp_path in global environment
library(Rcpp)
Rcpp::sourceCpp("functions/compute_lik_grad_full.cpp")

# bilinear interpolation gradient function
bilinearGradVec <- function(loc_mat, cov_list) {
  x_grid <- cov_list[[1]]$x
  y_grid <- cov_list[[1]]$y
  n_cov <- length(cov_list)
  n_obs <- nrow(loc_mat)
  
  ix <- findInterval(loc_mat[,1], x_grid)
  iy <- findInterval(loc_mat[,2], y_grid)
  
  valid <- ix > 0 & ix < length(x_grid) & iy > 0 & iy < length(y_grid)
  
  x1 <- x_grid[ix]
  x2 <- x_grid[ix + 1]
  y1 <- y_grid[iy]
  y2 <- y_grid[iy + 1]
  
  dx <- x2 - x1
  dy <- y2 - y1
  lx <- loc_mat[,1]
  ly <- loc_mat[,2]
  
  grad_array <- array(NA_real_, dim = c(n_cov, n_obs, 2))
  
  for (j in seq_len(n_cov)) {
    f11 <- cov_list[[j]]$z[cbind(ix,     iy)]
    f21 <- cov_list[[j]]$z[cbind(ix + 1, iy)]
    f12 <- cov_list[[j]]$z[cbind(ix,     iy + 1)]
    f22 <- cov_list[[j]]$z[cbind(ix + 1, iy + 1)]
    
    dfdx <- ((y2 - ly) * (f21 - f11) + (ly - y1) * (f22 - f12)) / (dy * dx)
    dfdy <- ((x2 - lx) * (f12 - f11) + (lx - x1) * (f22 - f21)) / (dy * dx)
    
    grad_array[j, valid, 1] <- dfdx[valid]
    grad_array[j, valid, 2] <- dfdy[valid]
  }
  
  grad_array
}

# vectorized and paralellized likelihood and gradient function using analytical gradient and no precomputed 
# Generalized likelihood function for fixed/regular sampling interval
lik_grad_regular <- function(par, cl, n_cov, chol_m,
                             X, covlist, delta, N, M, mu_x_all, mu_y_all, B, P, cpp_path,
                             ...){  
  
  l = 0
  lik_grad = rep(0, n_cov + 1)
  
  sigma <- diag(delta * par[n_cov + 1] / (N + 1), 2, 2)
  delta_par <- delta * par[n_cov + 1] / (2 * (N + 1))
  
  gamma = sqrt(par[n_cov + 1])
  chol_matrix = gamma * chol_m
  
  compute <- function(i){
    mu_x <- mu_x_all[((i - 1) * N + 1):(i * N)]
    mu_y <- mu_y_all[((i - 1) * N + 1):(i * N)]
    
    x_samples = sweep(B[1, i, 1:M, 1:N] * gamma, 2, mu_x, "+")
    y_samples = sweep(B[2, i, 1:M, 1:N] * gamma, 2, mu_y, "+")
    
    # bridge contribution to likelihood (log-scale)
    l_k <- P[i, 1:M] + log(gamma) * (2 * N)
    
    full_x <- cbind(X[i, 1], x_samples, X[i + 1, 1])
    full_y <- cbind(X[i, 2], y_samples, X[i + 1, 2])
    
    
    res_mat <- compute_log_lik_grad_full_cpp(full_x, full_y, l_k, X, i, par, delta, N, covlist)
    
    logw <- res_mat[1, ]
    a    <- max(logw)
    w    <- exp(logw - a)
    Wsum <- sum(w)
    if (!is.finite(Wsum) || Wsum <= 0) return(rep(0, n_cov + 2L))
    
    loglike_i <- a + log(Wsum) - log(M)
    
    norm_w <- w / Wsum
    
    # Compute gradients for all beta coefficients
    grad_beta <- numeric(n_cov)
    for(j in 1:n_cov) {
      grad_beta[j] <- -sum(res_mat[1 + j, ] * norm_w) / 2
    }
    
    # Gradient for variance parameter
    grad_sigma <- -sum(res_mat[n_cov + 2, ] * norm_w)
    
    return(c(loglike_i, grad_beta, grad_sigma))
  }
  
  results <- parLapply(cl, 1:(nrow(X) - 1), compute)
  
  # Extract results
  results_mat <- do.call(rbind, results)
  l = sum(results_mat[, 1])
  for(j in 1:(n_cov + 1)) {
    lik_grad[j] = sum(results_mat[, j + 1])
  }
  
  if(is.nan(l)){
    print("NaN")
    return(list(l = 1e10, g = rep(0, n_cov + 1)))
  }
  
  if(is.infinite(l)){
    print("Inf")
    return(list(l = 1e10, g = rep(0, n_cov + 1)))
  } else {
    return(list(l = -l, g = lik_grad))
  }
}

#cached version of likelihood and gradient function
cached_lik_grad_regular <- function(cl, n_cov, chol_m,
                                    X, covlist, delta, N, M, mu_x_all, mu_y_all, B, P, cpp_path,
                                    ...) {
  last_par <- NULL
  last_res <- NULL
  
  eval_lik_grad <- function(par) {
    if (is.null(last_par) || length(par) != length(last_par) || any(par != last_par)) {
      last_par <<- par
      last_res <<- lik_grad_regular(
        par, cl, n_cov, chol_m,
        X, covlist, delta, N, M, mu_x_all, mu_y_all, B, P, cpp_path, ...
      )
    }
    last_res
  }
  
  fn <- function(par) eval_lik_grad(par)$l
  gr <- function(par) eval_lik_grad(par)$g
  
  list(fn = fn, gr = gr)
}


# Generalized likelihood function for iregular sampling interval
lik_grad_irregular <- function(par, cl,
                               X, covlist, times, ID, dt_max, M, B, P, bilinearGradVec, cpp_path,
                               ...){  
  p <- length(par) - 1L
  if (p < 1L) stop("par must be [betas..., s] with at least one covariate.")
  
  gamma <- sqrt(par[p + 1L])
  
  compute <- function(i){
    if (ID[i] != ID[i + 1]) return(rep(0, p + 2L))
    
    N     <- as.numeric(ceiling((times[i + 1] - times[i]) / dt_max))
    delta <- as.numeric(times[i + 1] - times[i])
    
    if (N == 1L) {
      y <- c(X[i + 1, 1] - X[i, 1],
             X[i + 1, 2] - X[i, 2])
      
      ## gradients at the start location
      G_arr <- bilinearGradVec(matrix(X[i, ], nrow = 1), covlist)
     
      G_i <- t(drop(G_arr[, 1, ]))  # 2 x p
      
      
      if (!all(is.finite(G_i))){
        print("inf")
        return(rep(0, p + 2L))
      }
      
      beta <- par[1:p]
      s    <- gamma^2
      
      mu <- as.vector((delta * s / 2) * (G_i %*% beta))  # length 2
      e  <- y - mu
      
      ## log-likelihood for d=2, Σ = δ s I_2
      loglike_i <- -log(2 * pi) - log(delta * s) - sum(e * e) / (2 * delta * s)
      
      ## scores: wrt beta and s (= gamma^2)
      g_beta <- 0.5 * as.numeric(crossprod(G_i, e))  # p-vector
      g_s    <- -1 / s + (sum(e * e)) / (2 * delta * s^2) + as.numeric(t(e) %*% (G_i %*% beta)) / (2 * s)
      
      if (!is.finite(loglike_i) || !all(is.finite(g_beta)) || !is.finite(g_s)) {
        return(rep(0, p + 2L))
      }
      return(c(loglike_i, -g_beta, -g_s))
    }
    
    
    # brownian bridge endpoints 
    mu_x <- rep(X[i, 1], each = N) + 1:N * rep((X[i + 1, 1] - X[i, 1]), each = N) / (N + 1)
    mu_y <- rep(X[i, 2], each = N) + 1:N * rep((X[i + 1, 2] - X[i, 2]), each = N) / (N + 1)
    
    #scaled bridges
    x_samples <- sweep(B[[i]][1, 1:M, 1:N] * gamma, 2, mu_x, "+")
    y_samples <- sweep(B[[i]][2, 1:M, 1:N] * gamma, 2, mu_y, "+")
    
    #bridge contribtution to likelihood (log-scale)
    l_k <- (P[i, 1:M]) + log(gamma) * (2 * N)
    
    #bridges with endpoints
    full_x <- cbind(X[i, 1], x_samples, X[i + 1, 1])
    full_y <- cbind(X[i, 2], y_samples, X[i + 1, 2])

    #likelihood and gradient calculations using c++ function
    res_mat <- compute_log_lik_grad_full_cpp(full_x, full_y, l_k, X, i, par, delta, N, covlist)

    logw <- res_mat[1, ]
    a    <- max(logw)
    w    <- exp(logw - a)
    Wsum <- sum(w)
    if (!is.finite(Wsum) || Wsum <= 0) return(rep(0, p + 2L))

    loglike_i <- a + log(Wsum) - log(M)

    norm_w <- w / Wsum

    g_beta <- numeric(p)
    for (j in 1:p) {
      g_beta[j] <- -sum(res_mat[1 + j, ] * norm_w) / 2
    }

    g_s <- -sum(res_mat[p + 2L, ] * norm_w)

    c(loglike_i, g_beta, g_s)

    
  }
  results_list <- parLapply(cl, 1:(nrow(X) - 1L), compute)
  res_mat <- do.call(rbind, results_list)  
  
  # fold
  l_sum      <- sum(res_mat[, 1], na.rm = TRUE)
  g_betasum  <- colSums(res_mat[, 2:(p + 1L), drop = FALSE], na.rm = TRUE)
  g_ssum     <- sum(res_mat[, p + 2L], na.rm = TRUE)
  
  if (!is.finite(l_sum)) {
    return(list(l = 1e10, g = rep(0, p + 1L)))
  }
  list(l = -l_sum, g = c(g_betasum, g_ssum))
}

cached_lik_grad_irregular <- function(cl, 
                                           X, covlist, times, ID, dt_max, M, B, P, bilinearGradVec, cpp_path,
                                           ...) {
  last_par <- NULL
  last_res <- NULL
  
  eval_lik_grad <- function(par) {
    if (is.null(last_par) || length(par) != length(last_par) || any(par != last_par)) {
      last_par <<- par
      last_res <<- lik_grad_irregular(
        par, cl, X, covlist, times, ID, dt_max, M, B, P, bilinearGradVec, cpp_path, ...
      )
    }
    last_res
  }
  
  fn <- function(par) eval_lik_grad(par)$l
  gr <- function(par) eval_lik_grad(par)$g
  
  list(fn = fn, gr = gr)
}


# Wrapper function to fit the model
fit_langevin_bbis <- function(X, covlist, delta, 
                              M = 50,
                              N = NULL, 
                              dt_max = NULL, 
                              dt_units = "secs",  # can be: "secs", "mins", "hours", "days", "weeks"
                              ncores = 10, 
                              init_beta = NULL, init_sigma = 1,
                              lower_sigma = 0.0001, 
                              cpp_path = "functions/compute_lik_grad_full.cpp", verpose = TRUE,
                              fixed_sampling = NULL) {
  
  # define likelihood function to use
  if(isTRUE(fixed_sampling)){
    if(is.null(N)) stop("must define number of nodes 'N' if fixed_sampling = TRUE")
    lik_grad <- cached_lik_grad_regular
  } else if(isFALSE(fixed_sampling)){
    if(is.null(dt_max)) stop("must define number of nodes 'dt_max' if fixed_sampling = FALSE")
    lik_grad <- cached_lik_grad_irregular
  } else if(is.null(fixed_sampling)) {
    message("add code to fit_langevin_bbis to check if data is fixed or random interval")
    browser()
  }
  
  # prep data
  if ("ID" %!in% names(X)) {
    ID <- 1  # add dummy ID
  } else {
    ID <- as.integer(X$ID)
  }
  
  # prep time
  if (!fixed_sampling) {
    if("time" %!in% names(X)) stop("'time' must be column in 'X'.")
    # check if time is time formaet
    if(!inherits(X$time, "POSIXct")) stop("'time' must be class 'POSIXct'")
    times <- (X$time - min(X$time))
    units <- list(secs = 1, mins = 60, hours = 3600, days = 86400, weeks = 604800)
    times <- times/units[[dt_units]]
  }
  
  # prep location data
  X <- matrix(c(X$x, X$y), ncol = 2) 
  
  # prep covariates
  if (inherits(covlist, "SpatRaster")) {
    # convert from spatraster
    covlist <- spatRast_to_covlist(covlist)
  } 
  
  # Determine number of covariates
  n_cov <- length(covlist)
  
  # Set initial parameters if not provided
  if(is.null(init_beta)) {
    init_beta <- rep(0, n_cov)
  }
  
  if(length(init_beta) != n_cov) {
    stop("init_beta length must match number of covariates")
  }
  
  par_init <- c(init_beta, init_sigma)
  
  # Brownian bridge covariance matrix
  if (fixed_sampling) {
    sigma_matrix <- delta * outer(1 - 1:N/(N+1), 1:N/(N+1))
    sigma_matrix <- lower.tri(sigma_matrix, TRUE) * sigma_matrix +
      t(lower.tri(sigma_matrix) * sigma_matrix)
    chol_m <- chol(sigma_matrix)
    
    # Brownian bridge endpoints
    mu_x_all <- rep(X[1:(nrow(X)-1), 1], each = N) + 
      1:N * rep((X[2:nrow(X), 1] - X[1:(nrow(X)-1), 1]), each = N) / (N+1)
    mu_y_all <- rep(X[1:(nrow(X)-1), 2], each = N) + 
      1:N * rep((X[2:nrow(X), 2] - X[1:(nrow(X)-1), 2]), each = N) / (N+1)
    
    # Brownian bridge array
    B <- array(data = NA, c(2, nrow(X)-1, M, N))
    # Importance sampling weights
    P <- array(data = NA, c(nrow(X)-1, M))
    
    # Generate bridges
    if(verpose) cat("Generating Brownian bridges...\n")
    for (i in 1:(nrow(X)-1)) {
      mu_x <- mu_x_all[((i - 1) * N + 1):(i * N)]
      mu_y <- mu_y_all[((i - 1) * N + 1):(i * N)]
      
      B[1, i, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0, N), sigma = chol_m, isChol = TRUE)
      B[2, i, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0, N), sigma = chol_m, isChol = TRUE)
      
      P[i, 1:M] <- -(mvnfast::dmvn(B[1, i, 1:M, 1:N], rep(0, N), sigma = chol_m, isChol = TRUE, log = TRUE) + 
                        mvnfast::dmvn(B[2, i, 1:M, 1:N], rep(0, N), sigma = chol_m, isChol = TRUE, log =TRUE))
    }
  } else {  # irregular sampling interval
    #brownian bridge array
    B <- list()
    #importance sampling wights
    P <- array(data = NA, c(nrow(X)-1,M))
    #generating bridges
    for (i in 1:(nrow(X)-1)) {
      if(ID[i] == ID[i+1]){
        N = as.numeric(ceiling((times[i+1] - times[i])/dt_max))
        if(N != 1){
          delta = (times[i+1] - times[i])
          
          #brownian bridge covariance matrix
          sigma_matrix <- delta * outer(1 - 1:N/(N+1), 1:N/(N+1))
          sigma_matrix <- lower.tri(sigma_matrix, TRUE) * sigma_matrix +
            t(lower.tri(sigma_matrix) * sigma_matrix)
          chol_m <- chol(sigma_matrix)
          
          b <- array(data = NA, c(2, M, N))
          
          # Generate all M sample tracks at once
          b[1, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0,N), sigma = chol_m, isChol = TRUE)
          b[2, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0,N), sigma = chol_m, isChol = TRUE)
          
          B[[i]] = b
          
          P[i, 1:M] = -(mvnfast::dmvn(b[1, 1:M, 1:N], rep(0,N), sigma = chol_m, isChol = TRUE, log = TRUE) + 
                           mvnfast::dmvn(b[2, 1:M, 1:N], rep(0,N), sigma = chol_m, isChol = TRUE, log = TRUE))
        }
      } 
    }
  }
  
  # Set up parallel cluster
  if(verpose) cat("Setting up parallel cluster with", ncores, "cores...\n")
  cl <- makeCluster(ncores)
  
  # Source C++ file
  # Export variables to cluster
  if(fixed_sampling) {
    varlist <- c("X", "N", "M", "mu_x_all", "mu_y_all", 
                 "chol_m", "covlist", "delta", "B", "P", "n_cov", "cpp_path")
  } else {
    varlist <- c("X",  "M", "covlist", "bilinearGradVec", "B", "P", 
                 "times", "ID", "dt_max", "cpp_path")
  }
  clusterExport(cl, varlist = varlist, envir = environment())
  
  # Load required libraries on workers
  clusterEvalQ(cl, {
    library(mvnfast)
    library(Rcpp)
    sourceCpp(cpp_path)
  })
  
  # Run optimization
  if(verpose) cat("Running optimization...\n")
  t_start <- Sys.time()
  
  # Set bounds: no bounds on betas, lower bound on sigma
  lower_bounds <- c(rep(-Inf, n_cov), lower_sigma)
  if (fixed_sampling) {
    cached <- cached_lik_grad_regular(
      cl = cl, n_cov = n_cov, chol_m = chol_m,
      X = X, covlist = covlist, delta = delta, N = N, M = M,
      mu_x_all = mu_x_all, mu_y_all = mu_y_all, B = B, P = P,
      cpp_path = cpp_path
    )
  } else {
    cached <- cached_lik_grad_irregular(
      cl = cl,
      X = X, covlist = covlist, times = times, ID = ID, dt_max = dt_max,
      M = M, B = B, P = P, bilinearGradVec = bilinearGradVec,
      cpp_path = cpp_path
    )
  }
  
  o <- optim(
    par = par_init,
    fn = cached$fn,
    gr = cached$gr,
    method = "L-BFGS-B",
    lower = lower_bounds
  )
  
  t_elapsed <- Sys.time() - t_start
  
  # Stop cluster
  stopCluster(cl)
  
  # Prepare results
  results <- list(
    par = o$par,
    beta = o$par[1:n_cov],
    sigma = o$par[n_cov + 1],
    convergence = o$convergence,
    value = o$value,
    counts = o$counts,
    time = t_elapsed,
    delta = delta,
    N = N,
    M = M,
    n_cov = n_cov
  )
  
  # Print summary
  if(verpose) {
    cat("\n--- Optimization Results ---\n")
    cat("Time elapsed:", format(t_elapsed), "\n")
    cat("Convergence:", o$convergence, "\n")
    cat("Function evaluations:", o$counts[1], "\n")
    cat("Gradient evaluations:", o$counts[2], "\n")
    cat("\nEstimated parameters:\n")
    for(i in 1:n_cov) {
      cat(sprintf("  beta[%d] = %.4f\n", i, results$beta[i]))
    }
    cat(sprintf("  sigma   = %.4f\n", results$sigma))
  }
  
  return(results)
}




