source("global/all_function_allAge.R")

# load chains
n_sts <- 1100
dir_name <- paste0("outputs/genomics/trial_stochastic_", n_sts, "/")
dir.create(paste0(dir_name, "/figs"), FALSE, TRUE)
mcmc1 <- read.csv(paste0(dir_name, "mcmc1.csv"))
mcmc2 <- read.csv(paste0(dir_name, "mcmc2.csv"))
mcmc2_burnedin <- read.csv(paste0(dir_name, "mcmc2_burnedin.csv"))
tune_pmcmc_result <- readRDS(paste0(dir_name, "tune_pmcmc_result.rds"))

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


# MCMC diagnostics
# 1. Gelman-Rubin
figs_gelman_init <- diag_init_gelman_rubin(tune_pmcmc_result)
# fig <- diag_cov_mtx(figs_gelman_init)
# fig <- diag_gelman_rubin(figs_gelman_init)

png(paste0(dir_name, "figs/mcmc2_diag_gelmanRubin_%02d.png"),
    width = 17, height = 17, unit = "cm", res = 600)
diag_gelman_rubin(figs_gelman_init)
dev.off()

# 2. Autocorrelation
png(paste0(dir_name, "figs/mcmc2_diag_auCorr_%02d.png"),
    width = 17, height = 17, unit = "cm", res = 600)
diag_aucorr(coda::as.mcmc(mcmc2))
dev.off()

# 3. ggpairs
png(paste0(dir_name, "figs/mcmc2_diag_ggPairs_%02d.png"),
    width = 17, height = 17, unit = "cm", res = 600)
GGally::ggpairs(as.data.frame(tune_pmcmc_result$pars))
dev.off()