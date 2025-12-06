source("global/all_function_allAge.R")

# generate post-mcmc picts
post_pmcmc_pics <- function(n_sts){
  dir_name <- paste0("outputs/genomics/trial_", n_sts, "/")
  dir.create(paste0(dir_name, "/figs"), FALSE, TRUE)
  mcmc1 <- read.csv(paste0(dir_name, "mcmc1.csv"))
  mcmc2 <- read.csv(paste0(dir_name, "mcmc2.csv"))
  mcmc2_burnedin <- read.csv(paste0(dir_name, "mcmc2_burnedin.csv"))
  # tune_pmcmc_result <- readRDS(paste0(dir_name, "tune_pmcmc_result.rds"))
  
  # mcmc1
  # fig <- pmcmc_trace(coda::as.mcmc(mcmc1))
  png(paste0(dir_name, "figs/mcmc1_%02d.png"),
      width = 17, height = 17, unit = "cm", res = 600)
  pmcmc_trace(coda::as.mcmc(mcmc1))
  dev.off()
  
  # mcmc2
  # fig <- pmcmc_trace(coda::as.mcmc(mcmc2))
  png(paste0(dir_name, "figs/mcmc2_%02d.png"),
      width = 17, height = 17, unit = "cm", res = 600)
  pmcmc_trace(coda::as.mcmc(mcmc2))
  dev.off()
  
  # mcmc2 burned in
  # fig <- pmcmc_trace(coda::as.mcmc(mcmc2))
  png(paste0(dir_name, "figs/mcmc2_burnedin_%02d.png"),
      width = 17, height = 17, unit = "cm", res = 600)
  pmcmc_trace(coda::as.mcmc(mcmc2_burnedin))
  dev.off()
  
  # final parameters with CI
  # tune_lpost_max <- which.max(tune_pmcmc_result$probabilities[, "log_posterior"])
  # mcmc_lo_CI <- apply(tune_pmcmc_result$pars, 2, function(x) quantile(x, probs = 0.025))
  # mcmc_hi_CI <- apply(tune_pmcmc_result$pars, 2, function(x) quantile(x, probs = 0.975))
  # 
  # binds_tune_initial <- rbind(as.list(tune_pmcmc_result$pars[tune_lpost_max, ]),
  #                             mcmc_lo_CI, mcmc_hi_CI)
  # binds_tune_initial2 <- cbind(binds_tune_initial,
  #                              log_prior = tune_pmcmc_result$probabilities[tune_lpost_max,
  #                                                                          "log_prior"],
  #                              log_likelihood = tune_pmcmc_result$probabilities[tune_lpost_max,
  #                                                                               "log_likelihood"],
  #                              log_posterior = tune_pmcmc_result$probabilities[tune_lpost_max,
  #                                                                              "log_posterior"])
  # t_tune_initial <- t(binds_tune_initial2)
  # colnames(t_tune_initial) <- c("values", "low_CI", "high_CI")
  # 
  # write.csv(t_tune_initial,
  #           paste0(dir_name, "tune_initial_with_CI.csv"), row.names = T)
  
  # MCMC diagnostics
  # 1. Gelman-Rubin
  # figs_gelman_init <- diag_init_gelman_rubin(tune_pmcmc_result)
  # fig <- diag_cov_mtx(figs_gelman_init)
  # fig <- diag_gelman_rubin(figs_gelman_init)
  
  # png(paste0(dir_name, "figs/mcmc2_diag_gelmanRubin_%02d.png"),
  #     width = 17, height = 17, unit = "cm", res = 600)
  # diag_gelman_rubin(figs_gelman_init)
  # dev.off()
  
  # 2. Autocorrelation
  png(paste0(dir_name, "figs/mcmc2_diag_auCorr_%02d.png"),
      width = 17, height = 17, unit = "cm", res = 600)
  diag_aucorr(coda::as.mcmc(mcmc2))
  dev.off()
  
  # 3. ggpairs
  # png(paste0(dir_name, "figs/mcmc2_diag_ggPairs_%02d.png"),
  #     width = 17, height = 17, unit = "cm", res = 600)
  # p <- GGally::ggpairs(as.data.frame(tune_pmcmc_result$pars))
  # print(p)
  # dev.off()
}

# a slight modification for college's HPC
args <- commandArgs(trailingOnly = T)
n_sts <- as.numeric(args[which(args == "--n_steps") + 1])

post_pmcmc_pics(n_sts)

