# 2. Data Fitting ##############################################################
library(mcstate)
library(coda)
library(odin.dust)
library(dust)
library(GGally)
library(socialmixr)

source("global/all_function.R")

# The anatomy of an mcstate particle filter, as noted above, consists of three main components: \n 
# 1. A set of observations to fit the model to, generated using mcstate::particle_filter_data(). \n 
# 2. A model to fit, which must be a dust generator, either dust::dust() or odin.dust::odin_dust(). \n 
# 3. A comparison function, which is an R function which calculates the likelihood of the state given the data at one time point.

# There is a calibration function in mcstate to fit our model to data.
# https://mrc-ide.github.io/mcstate/articles/sir_models.html

# To make my life easier I compile the Serotype 1 cases into a new object called sir_data
# data is fed as an input to mcstate::particle_filter_data
incidence <- read.csv("inputs/incidence_week_12F_allAge.csv") # %>% 
  # dplyr::mutate(across(everything(), ~ tidyr::replace_na(.x, 0)))

dt <- (1/7) # rate must be an integer; 0.25 to make it 4 days, I make it 1/7
sir_data <- mcstate::particle_filter_data(data = incidence,
                                          time = "week",
                                          rate = 1 / dt,
                                          initial_time = 0) # Initial time makes t0 start from 0 (not 1)

# Annotate the data so that it is suitable for the particle filter to use
rmarkdown::paged_table(sir_data)

# Contact matrix:
# Create contact_matrix 5 demographic groups:
# > 2
# 2-64
# 65+
# age.limits = c(0, 2, 65)
# N_age <- length(age.limits)
# 
# contact_demographic <- socialmixr::contact_matrix(polymod,
#                                                   countries = "United Kingdom",
#                                                   age.limits = age.limits,
#                                                   symmetric = TRUE
# )
# 
# transmission <- contact_demographic$matrix /
#   rep(contact_demographic$demography$population, each = ncol(contact_demographic$matrix))
# transmission


## 2a. Model Load ##############################################################
# The model below is stochastic, closed system SADR model that I have created before
# I updated the code, filled the parameters with numbers;
# e.g.dt <- user(0) because if dt <- user() generates error during MCMC run
gen_sir <- odin.dust::odin_dust("model/sir_stochastic.R")

# This is part of sir odin model:
pars <- list(log_A_ini = (-5.69897), # S_ini*10^(log10(-5.69897)) = 120 people; change A_ini into log10(A_ini)
             time_shift_1 = 0.2,
             time_shift_2 = 0.2,
             beta_0 = 0.06565,
             beta_1 = 0.07, # in toy data the real value of beta_1 = 0.07
             beta_2 = 0.2,
             max_wane = (-0.5),
             min_wane = (-4),
             scaled_wane = (0.5),
             log_delta = (-4.98),
             sigma_2 = 1
)

# https://mrc-ide.github.io/odin-dust-tutorial/mcstate.html#/the-model-over-time
# n_particles <- 50 # Trial n_particles = 50
# filter <- mcstate::particle_filter$new(data = sir_data,
#                                        model = gen_sir, # Use odin.dust input
#                                        n_particles = n_particles,
#                                        compare = case_compare,
#                                        seed = 1L)
# 
# filter$run(pars)

# Variance and particles estimation (as suggested by Rich)
# parallel::detectCores() # 4 cores
# x <- replicate(30, filter$run(pars))
# var(x)
# [1] 3520.937
# Trial 320000 particles to get var(x) = 1 on 4 chains/4 nodes/4 cores per-node
# (320000/4/4/4)/3520.937
# Trial 320000 particles to get var(x) = 1 on 4 chains/1 nodes/20 cores per-node
# 3520.937*(4*1*20) # ~ 281675
# (281675/4/1/20)/3520.937

# Update n_particles based on calculation in 4 cores with var(x) ~ 3520.937: 281675
# 
priors <- prepare_priors(pars)
proposal_matrix <- diag(300, 8)
proposal_matrix <- (proposal_matrix + t(proposal_matrix)) / 2
rownames(proposal_matrix) <- c("log_A_ini", "time_shift_1", "time_shift_2", "beta_0", "beta_1", "beta_2", "scaled_wane", "log_delta")
colnames(proposal_matrix) <- c("log_A_ini", "time_shift_1", "time_shift_2", "beta_0", "beta_1", "beta_2", "scaled_wane", "log_delta")

transform <- function(pars) {
  parameter_transform(pars)
}
mcmc_pars <- prepare_parameters(initial_pars = pars,
                                priors = priors,
                                proposal = proposal_matrix,
                                transform = transform)

# Check the transform function:
# transform(mcmc_pars$initial())

# n_steps <- 100 #1e6

# I change pmcmc_run into a function that involve control inside:
# pmcmc_run <- mcstate::pmcmc(mcmc_pars, filter_deterministic, control = control)

# Directory for saving the outputs
dir.create("outputs/non_heterogeneity/trial_deterministic_1e3/figs", FALSE, TRUE)

# Trial combine pMCMC + tuning #################################################
pmcmc_run_plus_tuning <- function(n_pars, n_sts){
  filter <- mcstate::particle_filter$new(data = sir_data,
                                         model = gen_sir, # Use odin.dust input
                                         n_particles = n_pars,
                                         compare = case_compare,
                                         seed = 1L)
  
  # Use deterministic model by add filter_deterministic
  # https://mrc-ide.github.io/mcstate/articles/deterministic.html
  # Index function is optional when only a small number of states are used in comparison function.
  filter_deterministic <- mcstate::particle_deterministic$new(data = sir_data,
                                                              model = gen_sir,
                                                              compare = case_compare,
                                                              index = index_fun
  )
  
  
  control <- mcstate::pmcmc_control(n_steps = n_sts,
                                    rerun_every = 50,
                                    rerun_random = TRUE,
                                    progress = TRUE)
  
  # The pmcmc
  pmcmc_result <- mcstate::pmcmc(mcmc_pars, filter_deterministic, control = control)
  pmcmc_result
  saveRDS(pmcmc_result, "outputs/heterogeneity/pmcmc_result.rds")
  
  new_proposal_mtx <- cov(pmcmc_result$pars)
  write.csv(new_proposal_mtx, "outputs/heterogeneity/new_proposal_mtx.csv", row.names = TRUE)
  
  lpost_max <- which.max(pmcmc_result$probabilities[, "log_posterior"])
  write.csv(as.list(pmcmc_result$pars[lpost_max, ]),
            "outputs/heterogeneity/initial.csv", row.names = FALSE)
  
  # Further processing for thinning chains
  mcmc1 <- pmcmc_further_process(n_sts, pmcmc_result)
  write.csv(mcmc1, "outputs/heterogeneity/mcmc1.csv", row.names = TRUE)
  
  # Calculating ESS & Acceptance Rate
  calc_ess <- ess_calculation(mcmc1)
  write.csv(calc_ess, "outputs/heterogeneity/calc_ess.csv", row.names = TRUE)
  
  # Figures! (still failed, margin error)
  fig <- pmcmc_trace(mcmc1)
  # trial recursively save figs
  png("outputs/heterogeneity/trial_deterministic_5e3/figs/mcmc1_%02d.png", width = 17, height = 17, unit = "cm", res = 600)
  pmcmc_trace(mcmc1)
  dev.off()
  
  Sys.sleep(10) # wait 10 secs before conducting tuning
  
  # New proposal matrix
  new_proposal_matrix <- as.matrix(read.csv("outputs/heterogeneity/new_proposal_mtx.csv"))
  new_proposal_matrix <- new_proposal_matrix[, -1]
  new_proposal_matrix <- apply(new_proposal_matrix, 2, as.numeric)
  new_proposal_matrix <- new_proposal_matrix/1e3 # 100 resulted in bad chains while lower denominators resulted in jumpy steps among chains
  new_proposal_matrix <- (new_proposal_matrix + t(new_proposal_matrix)) / 2
  rownames(new_proposal_matrix) <- c("log_A_ini", "time_shift_1", "time_shift_2", "beta_0", "beta_1", "beta_2", "scaled_wane", "log_delta")
  colnames(new_proposal_matrix) <- c("log_A_ini", "time_shift_1", "time_shift_2", "beta_0", "beta_1", "beta_2", "scaled_wane", "log_delta")
  # isSymmetric(new_proposal_matrix)
  
  tune_mcmc_pars <- prepare_parameters(initial_pars = pars, priors = priors, proposal = new_proposal_matrix, transform = transform)
  
  # Including adaptive proposal control
  # https://mrc-ide.github.io/mcstate/reference/adaptive_proposal_control.html
  tune_control <- mcstate::pmcmc_control(n_steps = n_sts,
                                         n_chains = 4,
                                         rerun_every = 50,
                                         rerun_random = TRUE,
                                         progress = TRUE,
                                         adaptive_proposal = adaptive_proposal_control(initial_vcv_weight = 10,
                                                                                       initial_scaling = 0.7,
                                                                                       scaling_increment = NULL,
                                                                                       log_scaling_update = T,
                                                                                       acceptance_target = 0.234,
                                                                                       forget_rate = 0.1,
                                                                                       forget_end = n_sts*0.75,
                                                                                       adapt_end = n_sts*0.95,
                                                                                       pre_diminish = n_sts*0.1)
                                         )
  
  filter <- mcstate::particle_filter$new(data = sir_data,
                                         model = gen_sir, # Use odin.dust input
                                         n_particles = n_pars,
                                         compare = case_compare,
                                         seed = 1L
  )
  
  # The pmcmc
  tune_pmcmc_result <- mcstate::pmcmc(tune_mcmc_pars, filter_deterministic, control = tune_control)
  tune_pmcmc_result
  saveRDS(tune_pmcmc_result, "outputs/heterogeneity/tune_pmcmc_result.rds")
  
  new_proposal_mtx <- cov(pmcmc_result$pars)
  write.csv(new_proposal_mtx, "outputs/heterogeneity/new_proposal_mtx.csv", row.names = TRUE)
  
  tune_lpost_max <- which.max(tune_pmcmc_result$probabilities[, "log_posterior"])
  write.csv(as.list(tune_pmcmc_result$pars[tune_lpost_max, ]),
            "outputs/heterogeneity/tune_initial.csv", row.names = FALSE)
  
  # Further processing for thinning chains
  mcmc2 <- tuning_pmcmc_further_process(n_sts, tune_pmcmc_result)
  mcmc2 <- coda::as.mcmc(cbind(
    tune_pmcmc_result$probabilities, tune_pmcmc_result$pars))
  write.csv(mcmc2, "outputs/heterogeneity/mcmc2.csv", row.names = TRUE)
  
  # Calculating ESS & Acceptance Rate
  tune_calc_ess <- ess_calculation(mcmc2)
  write.csv(tune_calc_ess, "outputs/heterogeneity/tune_calc_ess.csv", row.names = TRUE)
  
  # Figures! (still failed, margin error)
  fig <- pmcmc_trace(mcmc2)
  
  png("outputs/heterogeneity/trial_deterministic_5e3/figs/mcmc2_%02d.png", width = 17, height = 17, unit = "cm", res = 600)
  pmcmc_trace(mcmc2)
  dev.off()
  
  ##############################################################################
  # MCMC Diagnostics
  
  # 1. Gelman-Rubin Diagnostic
  # https://cran.r-project.org/web/packages/coda/coda.pdf
  # png("pictures/diag_gelman_rubin.png", width = 17, height = 12, unit = "cm", res = 1200)
  figs_gelman_init <- diag_init_gelman_rubin(tune_pmcmc_result)
  fig <- diag_cov_mtx(figs_gelman_init)
  fig <- diag_gelman_rubin(figs_gelman_init)
  # dev.off()
  
  png("outputs/heterogeneity/trial_deterministic_5e3/figs/mcmc2_diag_gelmanRubin_%02d.png", width = 17, height = 17, unit = "cm", res = 600)
  diag_gelman_rubin(figs_gelman_init)
  dev.off()
  
  # 2. Autocorrelation
  # png("pictures/diag_aucorr.png", width = 17, height = 12, unit = "cm", res = 1200)
  fig <- diag_aucorr(mcmc2)
  # dev.off()
  
  png("outputs/heterogeneity/trial_deterministic_5e3/figs/mcmc2_diag_auCorr_%02d.png", width = 17, height = 17, unit = "cm", res = 600)
  diag_aucorr(mcmc2)
  dev.off()
  
  # png("outputs/heterogeneity/temporary_deterministic_1e3/figs/mcmc2_ggpairs_%03d.png", width = 20, height = 20, unit = "cm", res = 600)
  fig <- GGally::ggpairs(as.data.frame(tune_pmcmc_result$pars))
  # dev.off()
  
  png("outputs/heterogeneity/trial_deterministic_5e3/figs/mcmc2_diag_ggPairs_%02d.png", width = 17, height = 17, unit = "cm", res = 600)
  GGally::ggpairs(as.data.frame(tune_pmcmc_result$pars))
  dev.off()
  
}

# pmcmc_run_plus_tuning(320000, 1000)
