# See https://mrc-ide.github.io/mcstate/articles/nested_sir_models.html

case_compare <- function(state, observed, pars = NULL) {
  exp_noise <- 1e6
  n <- ncol(state)
  
  # incidence based on model's "n_AD_daily" from gen_sir$new(pars = list(), time = 0, n_particles = 1L)$info()
  # sir_model$info()$index$D
  incidence_modelled_1 <- state[13, , drop = TRUE] # D is extracted for each demographic group
  incidence_modelled_2 <- state[14, , drop = TRUE]
  incidence_modelled_3 <- state[15, , drop = TRUE]
  
  # incidence based on data already in x = observed$cases
  # lamb <- incidence_modelled + rexp(n, exp_noise)
  
  lambda_2 <- rexp(n, exp_noise)
  
  loglik_1 <- dpois(x = observed$cases_1,
                            lambda = incidence_modelled_1 + lambda_2,
                            log = T)
  loglik_2 <- dpois(x = observed$cases_2,
                        lambda = incidence_modelled_2 + lambda_2,
                        log = T)
  loglik_3 <- dpois(x = observed$cases_3,
                         lambda = incidence_modelled_3 + lambda_2,
                         log = T)
  
  loglik_cases <- loglik_1+loglik_2+loglik_3
  if (any(!is.finite(loglik_cases))) {
    # return -Inf to force rejection
    loglik_cases[!is.finite(loglik_cases)] <- -1e10
  }
  
  return(loglik_cases)
}

# generate index function
index_fun <- function(info){
  if (is.null(info$index)){
    info <- info[[1]]
  }
  list(run = unlist(info$index),
       state = unlist(info$index))
}

# That transform function
# https://github.com/mrc-ide/mcstate/blob/da9f79e4b5dd421fd2e26b8b3d55c78735a29c27/tests/testthat/test-if2.R#L40
# https://github.com/mrc-ide/mcstate/issues/184
parameter_transform <- function(transmission) {
  age.limits = c(0, 2, 65)
  N_age <- length(age.limits)
  contact_demographic <- socialmixr::contact_matrix(polymod,
                                                    countries = "United Kingdom",
                                                    age.limits = age.limits,
                                                    symmetric = TRUE
  )
  
  transmission <- contact_demographic$matrix /
    rep(contact_demographic$demography$population, each = ncol(contact_demographic$matrix))
  transmission
  
  transform <- function(pars){
    # re-define pars with pars that I really wanna fit only
    log_A_ini <- pars[paste0("log_A_ini_", 1:3)]
    time_shift_1 <- pars[["time_shift_1"]]
    time_shift_2 <- pars[["time_shift_2"]]
    beta_0 <- pars[["beta_0"]]
    beta_1 <- pars[["beta_1"]]
    beta_2 <- pars[["beta_2"]]
    scaled_wane <- pars[["scaled_wane"]]
    log_delta <- pars[["log_delta"]]
    psi <- pars[["psi"]]
    # sigma_2 <- pars[["sigma_2"]]
    
    pars <- list(log_A_ini = log_A_ini,
                 time_shift_1 = time_shift_1,
                 time_shift_2 = time_shift_2,
                 beta_0 = beta_0,
                 beta_1 = beta_1,
                 beta_2 = beta_2,
                 scaled_wane = scaled_wane,
                 log_delta = log_delta,
                 psi = psi
                 # sigma_2 = sigma_2
    )
    
    pars$N_ini <-  contact_demographic$demography$population
    pars$D_ini <-  c(0,0,0)
    pars$R_ini <-  c(0,0,0)
    pars$m <- transmission
    
    pars
  }
}

transform <- function(pars) {
  parameter_transform(pars)
}

# transform <- parameter_transform(transmission)

prepare_parameters <- function(initial_pars, priors, proposal, transform) {
  
  mcmc_pars <- mcstate::pmcmc_parameters$new(
    list(# mcstate::pmcmc_parameter("log_A_ini", (-5.69897), min = (-10), max = 0,
                                  # prior = priors$log_A_ini),
         mcstate::pmcmc_parameter("log_A_ini_1", -4, min = -10, max = 0,
                                  prior = priors$log_A_ini_1),
         mcstate::pmcmc_parameter("log_A_ini_2", -4, min = -10, max = 0,
                                  prior = priors$log_A_ini_2),
         mcstate::pmcmc_parameter("log_A_ini_3", -4, min = -10, max = 0,
                                  prior = priors$log_A_ini_3),
         mcstate::pmcmc_parameter("time_shift_1", 0.2, min = 0, max = 1,
                                  prior = priors$time_shift_1),
         mcstate::pmcmc_parameter("time_shift_2", 0.2, min = 0, max = 1,
                                  prior = priors$time_shift_2),
         mcstate::pmcmc_parameter("beta_0", 0.06565, min = 0, max = 0.8,
                                  prior = priors$beta_0),
         mcstate::pmcmc_parameter("beta_1", 0.07, min = 0, max = 0.8,
                                  prior = priors$beta_1),
         mcstate::pmcmc_parameter("beta_2", 0.07, min = 0, max = 0.8,
                                  prior = priors$beta_2),
         mcstate::pmcmc_parameter("scaled_wane", (0.5), min = (0), max = 1,
                                  prior = priors$scaled_wane),
         mcstate::pmcmc_parameter("log_delta", (-4.98), min = (-10), max = 0.7,
                                  prior = priors$log_delta),
         mcstate::pmcmc_parameter("psi", (1), min = (0), max = 1,
                                  prior = priors$psi)
         # mcstate::pmcmc_parameter("sigma_2", 1, min = 0, max = 10,
         #                          prior = priors$sigma_2)
    ),
    proposal = proposal,
    transform = transform)
  
}

prepare_priors <- function(pars) {
  priors <- list()
  
  priors$log_A_ini_1 <- function(s) {
    dunif(s, min = (-10), max = 0, log = TRUE)
  }
  priors$log_A_ini_2 <- function(s) {
    dunif(s, min = (-10), max = 0, log = TRUE)
  }
  priors$log_A_ini_3 <- function(s) {
    dunif(s, min = (-10), max = 0, log = TRUE)
  }
  priors$time_shift_1 <- function(s) {
    dunif(s, min = 0, max = 1, log = TRUE)
  }
  priors$time_shift_2 <- function(s) {
    dunif(s, min = 0, max = 1, log = TRUE)
  }
  priors$beta_0 <- function(s) {
    dgamma(s, shape = 1, scale = 0.1, log = TRUE)
  }
  priors$beta_1 <- function(s) {
    dgamma(s, shape = 1, scale = 0.1, log = TRUE)
  }
  priors$beta_2 <- function(s) {
    dgamma(s, shape = 1, scale = 0.1, log = TRUE)
  }
  priors$scaled_wane <- function(s) {
    dbeta(s, shape1 = 1.25, shape2 = 1.25, log = TRUE)
  }
  priors$log_delta <- function(s) {
    dunif(s, min = (-10), max = 0.7, log = TRUE)
  }
  priors$psi <- function(s) {
    dgamma(s, shape = 1, scale = 0.1, log = TRUE)
  }
  # priors$sigma_2 <- function(s) {
  #   dgamma(s, shape = 1, scale = 1, log = TRUE)
  # }
  priors
}

pmcmc_further_process <- function(n_steps, pmcmc_result) {
  processed_chains <- mcstate::pmcmc_thin(pmcmc_result, burnin = n_steps*0.8, thin = NULL)
  parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
  parameter_mean_hpd
  
  mcmc1 <- coda::as.mcmc(cbind(pmcmc_result$probabilities, pmcmc_result$pars))
  mcmc1
}

ess_calculation <- function(mcmc1){
  # compile par names & generate switch
  par_names <- colnames(mcmc1)
  
  ess_values <- sapply(par_names, function(p){
    trace <- mcmc1[, p]
    if (var(trace) == 0) {
      warning(sprintf("Parameter '%s' has zero variance. ESS set to NA.", p))
      return(NA)
    } else {
      return(coda::effectiveSize(trace))
    }
  })
  acceptance_rate = 1 - coda::rejectionRate(mcmc1)
  
  list(
    ess = ess_values,
    acceptance_rate = acceptance_rate
  )
}

pmcmc_trace <- function(mcmc1) {
  plot(mcmc1) # to save the figures into pdf
  # png("pictures/mcmc1.png", res = 1200)
  # plot(mcmc1)
  # dev.off()
}

################################################################################
# Tuning functions
tuning_pmcmc_further_process <- function(n_steps, tune_pmcmc_result) {
  processed_chains <- mcstate::pmcmc_thin(tune_pmcmc_result, burnin = n_steps*0.8, thin = 2)
  parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
  parameter_mean_hpd
  
  tune_pmcmc_result <- coda::as.mcmc(cbind(processed_chains$probabilities, processed_chains$pars))
  tune_pmcmc_result
}

################################################################################
# MCMC Diagnostics
# 1. Gelman-Rubin Diagnostic
# https://cran.r-project.org/web/packages/coda/coda.pdf

diag_init_gelman_rubin <- function(tune_pmcmc_result){
  n_chains <- 4 # tune_control$n_chains
  n_samples <- nrow(tune_pmcmc_result$pars)/n_chains
  
  # Split the parameter samples and probabilities by chains
  chains <- lapply(1:n_chains, function(i) {
    start <- (i - 1) * n_samples + 1
    end <- i * n_samples
    list(
      pars = tune_pmcmc_result$pars[start:end, ],
      probabilities = tune_pmcmc_result$probabilities[start:end, ]
    )
  })
  
  # Convert chains to mcmc objects
  mcmc_chains <- lapply(chains, function(chain) {
    as.mcmc(cbind(chain$probabilities, chain$pars))
  })
  
  # Combine chains into a list
  mcmc_chains_list <- do.call(mcmc.list, mcmc_chains)
  mcmc_chains_list
}

diag_cov_mtx <- function(mcmc_chains_list) {
  # print("Covariance matrix of mcmc2")
  cov(as.matrix(mcmc_chains_list))
}

diag_gelman_rubin <- function(mcmc_chains_list) {
  # print("Gelman-Rubin diagnostic")
  gelman_plot <- coda::gelman.plot(mcmc_chains_list,
                                   bin.width = 10,
                                   max.bins = 50,
                                   confidence = 0.95,
                                   transform = FALSE,
                                   autoburnin=TRUE,
                                   auto.layout = TRUE)
  # ask, col, lty, xlab, ylab, type, ...)
  
  coda::gelman.diag(mcmc_chains_list,
                    confidence = 0.95,
                    transform=FALSE,
                    autoburnin=TRUE,
                    multivariate=F) # Change multivariate = F instead of T
  
}

# 2. Autocorrelation plots
diag_aucorr <- function(mcmc2){
  for (name in colnames(mcmc2)){
    print(coda::acfplot(mcmc2[, name], main = name))
  }
}