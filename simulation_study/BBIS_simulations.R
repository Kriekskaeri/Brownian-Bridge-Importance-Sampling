# Simulation study
# prep workspace ---------------------------------------------------------- ####
library(here)
# custom functions
source(here("functions/utility_functions.R"))  # custom general perpose functions
sourceDir("functions")  # custom function to load all functions in folder
load_lib(mvnfast, parallel, terra, dplyr,
         ambient, Rcpp)  # custom function to install & load packages
cpp_path <- here("functions/compute_lik_grad_full.cpp")
Rcpp::sourceCpp(cpp_path)

output_path <- "simulation_study/outputs"
make_path(output_path)
# Define parameters up front ---------------------------------------------- ####
set.seed(123)

## track pars 
speed <- 5              # speed parameter for Langevin model
dt    <- 1/60       # temporal resolution of simulated tracks (fraction of 1 hour)
beta  <- c(4, 2, -0.1)  # covariate coefficients
loc0  <- c(0, 0)        # starting location of tracks

## default estimation pars
n_sim  <- 100    # number of simulations per simulation study
ncores <- parallel::detectCores() - 2 # number of cores used in parallel computations
thin   <- 100      # thinning
N      <- 1/(dt*15) - 1 # default nodes
M      <- 20       # default number of bridges
n_obs  <- 3*31*24 / dt # #hours / dt

## covariate pars
res  <- 1  # resolution of covariates 
ncov <- 2  # number of covariates
scal = 100 #kilometers
ext  <- c(-1, 1, -1, 1)*scal  # extent of study area
perlin_f <- 0.05  # Perlin noise frequency

# simulate covariates with Perlin noise ----------------------------------- ####
covlist <- list()
xgrid <- seq(ext[1], ext[2], by = res)
ygrid <- seq(ext[3], ext[4], by = res)
coords <- as.matrix(expand.grid(xgrid, ygrid))
for(i in 1:ncov) {
  vals <- 3*noise_perlin(c(length(xgrid), length(ygrid)), frequency = perlin_f)
  covlist[[i]] = list(x = xgrid, y = ygrid, z = matrix(vals, nrow = length(xgrid)))
}

# Include squared distance to centre of map as covariate
xgrid <- seq(ext[1], ext[2], by = res)
ygrid <- seq(ext[3], ext[4], by = res)
xygrid <- expand.grid(xgrid,ygrid)
dist2 <- ((xygrid[,1])^2+(xygrid[,2])^2)/(100)
covlist[[3]] <- list(x = xgrid, y = ygrid,
                     z = matrix(dist2, length(xgrid), length(ygrid)))

# define result data.frame column names
col_names <- c(
  "sim",
  "method",  # (eiler/bbis)
  "dt",
  "Tmax",
  "delta",
  "N",
  "M",
  "convergence",
  "iterations",
  "dt",
  paste0("beta", seq(length(beta))), 
  "gammasq"
)

# Sim 1: varying delta_t, fixed number of observations -------------------- ####
print("varying delta_t, fixed number of observations")
sim_var <- c(1.0, 2.0, 12.0, 24.0, 48.0) / dt
sim_results <- data.frame()  # refresh result target
for (ik in 1:n_sim) {
  
  #simulate track
  beta_sim <- beta
  n_obs_sim <- n_obs
  X_full <- simLMM(dt, speed, covlist, beta_sim, loc0, n_obs_sim)
  
  for (jk in seq_along(sim_var)) {
    # set up simulation parameters
    thin_sim <- sim_var[jk]
    dt_sim <- dt
    delta <- dt_sim * thin_sim
    N_sim <- delta * (60 / 15) - 1
    
    M_sim <- M
    Tmax <- n_obs_sim * thin_sim * dt_sim
    
    # thinning track
    X_thin = thinTrack(X_full, N_sim)
    
    # estimate with euler
    UD <- langevinUD(X_thin, (0:(nrow(X_thin) - 1)) * delta, 
                     grad_array = bilinearGradArray(X_thin, covlist))
    ## extract & store euler outputs  
    sim_results <- data.frame(ik, "euler",   # sim & method
                              dt_sim, Tmax,  # sim conditions
                              delta, N_sim, M_sim,   # fit conditions
                              1, NA,     # convergence, iterations
                              as.numeric( UD$time, units = "secs"),  # compute time
                              matrix(c(UD$betaHat, UD$gamma2Hat), nrow = 1)) |> 
      setNames(col_names) %>% 
      rbind(sim_results, .)
    
    # estimate with bbis
    X_thin <- data.frame(x = X_thin[, 1], y = X_thin[, 2])
    out <- fit_langevin_bbis(X_thin, covlist, delta, N = N_sim, M = M_sim,
                             ncores = ncores, fixed_sampling = TRUE)  
    
    # extract & store bbis outputs
    sim_results <- data.frame(ik, "bbis",   # sim, method
                           dt_sim, Tmax,  # sim conditions
                           delta, N_sim, M_sim,  # fit conditions
                           out$convergence,  # convergence
                           as.numeric((out$counts)[1]),  # iterations
                           as.numeric(out$time, units = "secs"),  # compute time
                           matrix(out$par, nrow = 1)) %>%   # estimates
      setNames(col_names) %>% 
      rbind(sim_results, .)
  }
  write.csv(sim_results, file = here(output_path,"varying_thin_estimates.csv"), 
            row.names = F)
}



# Sim 2: varying delta_t, fixed maximum time ------------------------------ ####
print("varying delta_t, fixed maximum time")
sim_var <- c(5, 10, 20, 50, 100) * 36
sim_results <- results_template  # refresh result target

for (ik in 1:n_sim) {
  for (jk in seq_along(sim_var)) {
    # set up simulation parameters
    beta_sim <- beta
    thin_sim <- sim_var[jk]
    dt_sim <- dt
    delta <- dt_sim*thin_sim
    N_sim <- thin_sim-1
    M_sim <- M
    Tmax <- 500
    n_obs_sim <- Tmax/(dt_sim*thin_sim)
    
    # simulating track
    X <- simLMM(delta, speed, covlist, beta_sim, loc0, n_obs_sim)
    
    # estimate with euler
    UD <- langevinUD(X, (0:(nrow(X) - 1)) * delta, 
                     grad_array = bilinearGradArray(X, covlist))
    ## extract & store euler outputs  
    sim_results <- data.frame(ik, "euler",   # sim & method
                              dt_sim, Tmax,  # sim conditions
                              delta, N_sim, M_sim,   # fit conditions
                              1, NA,     # convergence, iterations
                              as.numeric( UD$time, units = "secs"),  # compute time
                              matrix(c(UD$betaHat, UD$gamma2Hat), nrow = 1)) |> 
      setNames(col_names) %>% 
      rbind(sim_results, .)
    
    # fit model
    X = data.frame(x = X[,1], y = X[,2])
    out <- fit_langevin_bbis(X, covlist, delta, N = N_sim, M = M_sim,
                             ncores = ncores, cpp_path = cpp_path)  
    
    # extract bbis outputs
    sim_results <- data.frame(ik, "bbis",   # sim, method
                           dt_sim, Tmax,  # sim conditions
                           delta, N_sim, M_sim,  # fit conditions
                           out$convergence,  # convergence
                           as.numeric((out$counts)[1]),  # iterations
                           as.numeric(out$time, units = "secs"),  # compute time
                           matrix(out$par, nrow = 1)) |>   # estimates
      setNames(col_names) %>% 
      rbind(sim_results, .)
  }
  # save output
  write.csv(sim_results, 
            here(output_path, "varying_thin_estimates_fixed_Tmax.csv"),
            row.names = FALSE)
}

# Sim 3: varying number of bridges (M) ------------------------------------ ####
print("varying M")
sim_var <- c(5, 10, 50, 100, 200)
sim_results <- results_template  # refresh result target

for (ik in 1:n_sim) {
  beta_sim <- beta
  thin_sim <- thin
  dt_sim <- dt
  delta <- dt*thin
  N_sim <- 49
  n_obs_sim <- n_obs
  Tmax <- n_obs_sim*thin_sim*dt_sim
  
  # simulating track
  X <- simLMM(delta, speed, covlist, beta_sim, loc0, n_obs_sim)
  
  # estimate with euler
  UD <- langevinUD(X, (0:(nrow(X) - 1)) * delta, 
                   grad_array = bilinearGradArray(X, covlist))
  ## extract & store euler outputs  
  sim_results <- data.frame(ik, "euler",   # sim & method
                            dt_sim, Tmax,  # sim conditions
                            delta, N_sim, M_sim,   # fit conditions
                            1, NA,     # convergence, iterations
                            as.numeric( UD$time, units = "secs"),  # compute time
                            matrix(c(UD$betaHat, UD$gamma2Hat), nrow = 1)) |> 
    setNames(col_names) %>% 
    rbind(sim_results, .)
  
  # loop for BBIS
  for (jk in seq_along(sim_var)) {
    M_sim <- sim_var[jk]
    
    # fit model
    out <- fit_langevin_bbis(X, covlist, delta, N = N_sim, M = M_sim,
                             ncores = ncores, cpp_path = cpp_path)  
    
    # extract & store bbis outputs
    sim_results <- data.frame(ik, "bbis",   # sim, method
                           dt_sim, Tmax,  # sim conditions
                           delta, N_sim, M_sim,  # fit conditions
                           out$convergence,  # convergence
                           as.numeric((out$counts)[1]),  # iterations
                           as.numeric(out$time, units = "secs"),  # compute time
                           matrix(out$par, nrow = 1)) |>   # estimates
      setNames(col_names) %>% 
      rbind(sim_results, .)
  }
  # save output
  write.csv(sim_results, here(output_path, "varying_M_estimates.csv"),
            row.names = FALSE)
}

# Sim 4: varying number of nodes (N) -------------------------------------- ####
print("varying N")
sim_var <- c(4, 9, 49, 99)
sim_results <- results_template  # refresh result target

for (ik in 1:n_sim) {
  beta_sim <- beta
  thin_sim <- thin
  dt_sim <- dt
  delta <- dt_sim*thin_sim
  M_sim <- M
  n_obs_sim <- n_obs
  Tmax <- n_obs_sim*thin_sim*dt_sim
  # simulating track
  X <- simLMM(delta, speed, covlist, beta_sim, loc0, n_obs_sim)
  
  # estimate with euler
  UD <- langevinUD(X, (0:(nrow(X) - 1)) * delta, 
                   grad_array = bilinearGradArray(X, covlist))
  ## extract & store euler outputs  
  sim_results <- data.frame(ik, "euler",   # sim & method
                            dt_sim, Tmax,  # sim conditions
                            delta, N_sim, M_sim,   # fit conditions
                            1, NA,     # convergence, iterations
                            as.numeric( as.numeric( UD$time, units = "secs"), units = "secs"),  # compute time
                            matrix(c(UD$betaHat, UD$gamma2Hat), nrow = 1)) |> 
    setNames(col_names) %>% 
    rbind(sim_results, .)
  
  # loop for BBIS
  for (jk in seq_along(sim_var)) {
    N_sim <- sim_var[jk]
    
    # fit model
    out <- fit_langevin_bbis(X, covlist, delta, N = N_sim, M = M_sim,
                             ncores = ncores, cpp_path = cpp_path)  
    # extract & store bbis outputs
    bbis_out <- data.frame(ik, "bbis",   # sim, method
                           dt_sim, Tmax,  # sim conditions
                           delta, N_sim, M_sim_sim,  # fit conditions
                           out$convergence,  # convergence
                           as.numeric((out$counts)[1]),  # iterations
                           as.numeric(out$time, units = "secs"),  # compute time
                           matrix(out$par, nrow = 1)) |>   # estimates
      setNames(col_names) %>% 
      rbind(sim_results, .)
  }
  # save output
  write.csv(sim_results, here(output_path, "varying_N_estimates.csv"),
            row.names = FALSE)
}
