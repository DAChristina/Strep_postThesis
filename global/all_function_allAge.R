# See https://mrc-ide.github.io/mcstate/articles/nested_sir_models.html
ll_nbinom <- function(data, model, kappa, exp_noise) {
  if (is.na(data)) {
    return(numeric(length(model)))
  }
  mu <- model + rexp(length(model), rate = exp_noise)
  dnbinom(data, kappa, mu = mu, log = TRUE)
}

case_compare <- function(state, observed, pars = NULL) {
  exp_noise <- 1e6
  n <- ncol(state)
  
  # sir_model$info()$index$n_AD_weekly
  # model_55_all <- state[6, , drop = TRUE]
  model_55_1 <- state[9, , drop = TRUE]
  model_55_2 <- state[10, , drop = TRUE]
  model_55_3 <- state[11, , drop = TRUE]
  
  # if (is.na(observed$count_55_all)) {
  #   ll_55_all <- numeric(n)
  # } else {
  #   ll_55_all <- ll_nbinom(data = observed$count_55_all,
  #                          model = model_55_all,
  #                          kappa = pars$kappa_55,
  #                          exp_noise = exp_noise)
  # }
  
  if (is.na(observed$count_55_1)) {
    ll_55_1 <- numeric(n)
  } else {
    ll_55_1 <- ll_nbinom(data = observed$count_55_1,
                         model = model_55_1,
                         kappa = pars$kappa_55,
                         exp_noise = exp_noise)
  }
  
  if (is.na(observed$count_55_2)) {
    ll_55_2 <- numeric(n)
  } else {
    ll_55_2 <- ll_nbinom(data = observed$count_55_2,
                         model = model_55_2,
                         kappa = pars$kappa_55,
                         exp_noise = exp_noise)
  }
  
  if (is.na(observed$count_55_3)) {
    ll_55_3 <- numeric(n)
  } else {
    ll_55_3 <- ll_nbinom(data = observed$count_55_3,
                         model = model_55_3,
                         kappa = pars$kappa_55,
                         exp_noise = exp_noise)
  }
  
  ll <- ll_55_1 + ll_55_2 + ll_55_3
  
  if (any(!is.finite(ll))) {
    # return -Inf to force rejection
    ll[!is.finite(ll)] <- -1e10
  }
  return(ll)
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
parameter_transform <- function(pars) {
  log_A_ini <- pars[["log_A_ini"]]
  time_shift_1 <- pars[["time_shift_1"]]
  beta_0 <- pars[["beta_0"]]
  beta_1 <- pars[["beta_1"]]
  # log_delta <- pars[["log_delta"]]
  log_delta1 <- pars[["log_delta1"]]
  log_delta2 <- pars[["log_delta2"]]
  log_delta3 <- pars[["log_delta3"]]
  kappa_55 <- pars[["kappa_55"]]
  
  list(log_A_ini = log_A_ini,
       time_shift_1 = time_shift_1,
       beta_0 = beta_0,
       beta_1 = beta_1,
       # log_delta = log_delta,
       log_delta1 = log_delta1,
       log_delta2 = log_delta2,
       log_delta3 = log_delta3,
       kappa_55 = kappa_55
  )
  
}

transform <- function(pars) {
  parameter_transform(pars)
}

prepare_parameters <- function(initial_pars, priors, proposal, transform) {
  
  mcmc_pars <- mcstate::pmcmc_parameters$new(
    list(mcstate::pmcmc_parameter("log_A_ini", (0.55), min = (0), max = 1,
                                  prior = priors$log_A_ini),
         mcstate::pmcmc_parameter("time_shift_1", 0.1, min = 0, max = 1,
                                  prior = priors$time_shifts),
         mcstate::pmcmc_parameter("beta_0", 0.04, min = 0, max = 0.8,
                                  prior = priors$betas),
         mcstate::pmcmc_parameter("beta_1", 0.2, min = 0, max = 0.7,
                                  prior = priors$betas),
         # mcstate::pmcmc_parameter("log_delta", (-4.55), min = (-10), max = 0.7,
         #                          prior = priors$log_delta),
         mcstate::pmcmc_parameter("log_delta1", (-5), min = (-10), max = 0.7,
                                  prior = priors$log_delta),
         mcstate::pmcmc_parameter("log_delta2", (-4.8), min = (-10), max = 0.7,
                                  prior = priors$log_delta),
         mcstate::pmcmc_parameter("log_delta3", (-4.8), min = (-10), max = 0.7,
                                  prior = priors$log_delta),
         mcstate::pmcmc_parameter("kappa_55", 6, min = 0,
                                  prior = priors$kappas) #function(p) log(1e-10))
    ),
    proposal = proposal,
    transform = transform
    )
}

prepare_priors <- function(pars) {
  priors <- list()
  
  priors$log_A_ini <- function(s) {
    dbeta(s, shape1 = 3.5, shape2 = 3.5, log = TRUE)
  }
  priors$time_shifts <- function(s) {
    dunif(s, min = 0, max = 1, log = TRUE)
  }
  priors$betas <- function(s) {
    dgamma(s, shape = 1, scale = 0.1, log = TRUE)
  }
  priors$log_delta <- function(s) {
    stabledist::dstable(s, alpha = 2, beta = 0, gamma = 0.8, delta = -4.8, log = TRUE)
    # dunif(s, min = (-10), max = 0.7, log = TRUE)
  }
  priors$deltas <- function(s) {
    dunif(s, min = (0), max = 10, log = TRUE)
  }
  priors$kappas <- function(s) {
    # dunif(s, min = 0, log = TRUE)
    stabledist::dstable(s, alpha = 2, beta = 0, gamma = 1, delta = 5, log = TRUE)
  }
  
  priors
}


pmcmc_further_process <- function(n_steps, pmcmc_result) {
  processed_chains <- mcstate::pmcmc_thin(pmcmc_result, burnin = round(n_steps*0.5), thin = NULL)
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
}

################################################################################
# Tuning functions
tuning_pmcmc_further_process <- function(n_steps, tune_pmcmc_result) {
  processed_chains <- mcstate::pmcmc_thin(tune_pmcmc_result,
                                          burnin = round(n_steps*0.5),
                                          thin = 2)
  parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
  parameter_mean_hpd
  
  tune_pmcmc_result <- coda::as.mcmc(cbind(processed_chains$probabilities,
                                           processed_chains$pars))
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
    coda::as.mcmc(cbind(chain$probabilities, chain$pars))
  })
  
  # Combine chains into a list
  mcmc_chains_list <- do.call(coda::mcmc.list, mcmc_chains)
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
  
  coda::gelman.diag(mcmc_chains_list,
                    confidence = 0.95,
                    transform=FALSE,
                    autoburnin=TRUE,
                    multivariate=F)
  
}

# 2. Autocorrelation plots
diag_aucorr <- function(mcmc2){
  for (name in colnames(mcmc2)){
    print(coda::acfplot(mcmc2[, name], main = name))
  }
}

################################################################################
# Particle samples (adapted from Lilith's)
observe_pois <- function(lambda) {
  n_par <- nrow(lambda)
  n_obs <- ncol(lambda)
  ret <- vapply(seq_len(n_par), function(i) {
    rpois(n_obs, lambda[i, ])}, numeric(n_obs))
  t(ret)
}

observe <- function(pmcmc_samples) {
  
  state <- pmcmc_samples$trajectories$state
  pars <- apply(pmcmc_samples$pars, MARGIN = 1, pmcmc_samples$predict$transform)
  time <- pmcmc_samples$trajectories$time
  
  ## extract model outputs
  # model state refers to n_AD_weekly (not the D compartment)
  model_55 <- (state[9, , , drop = TRUE]+state[10, , , drop = TRUE]+state[11, , , drop = TRUE])
  
  observed <- list()
  observed$cases_child_GPSC55 <- observe_pois(model_55)
  
  abind::abind(c(list(state), observed), along = 1)
}

plot_states <- function(state, data) {
  col <- grey(0.3, 0.1)
  # model state refers to n_AD_weekly (not the D compartment)
  matplot(data$yearWeek, t((state["n_AD1_weekly", , -1]+state["n_AD2_weekly", , -1]+state["n_AD3_weekly", , -1])),
          type = "l", lty = 1, col = col,
          xlab = "", ylab = "GPSC55 cases")
  points(data$yearWeek, data$count_55_all, col = 3, pch = 20)
  points(data$yearWeek, data$count_12F_all, col = 4, type = "l")
  
  matplot(data$yearWeek, xlab = "", t(state["S", , -1]),
          type = "l", lty = 1, col = 2, ylab = "%", ylim = c(0, 6.7e7), yaxt = "n")
  axis(side = 2, at = seq(0, 6e7, length.out = 5),
       labels = seq(0, 100, length.out = 5))
  
  matlines(data$yearWeek, t(state["A", , -1]), lty = 1, col = 1)
  matlines(data$yearWeek, t(state["D", , -1]), lty = 1, col = 4)
  matlines(data$yearWeek, t(state["R", , -1]), lty = 1, col = 5)
  legend("right", bty = "n", fill = 2:4, legend = c("S", "A", "D", "R"))
  # 
  # matplot(data$yearWeek, xlab = "", t(state["I_tot", , -1]),
  #         type = "l", lty = 1, col = 3, ylab = "carriers")
}




