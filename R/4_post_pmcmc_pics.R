source("global/all_function_allAge.R")

# load chains
n_sts <- 2500
dir_name <- paste0("outputs/genomics/trial_", n_sts, "/")
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

# final parameters with CI
tune_lpost_max <- which.max(tune_pmcmc_result$probabilities[, "log_posterior"])
mcmc_lo_CI <- apply(tune_pmcmc_result$pars, 2, function(x) quantile(x, probs = 0.025))
mcmc_hi_CI <- apply(tune_pmcmc_result$pars, 2, function(x) quantile(x, probs = 0.975))

binds_tune_initial <- rbind(as.list(tune_pmcmc_result$pars[tune_lpost_max, ]), mcmc_lo_CI, mcmc_hi_CI)
t_tune_initial <- t(binds_tune_initial)
colnames(t_tune_initial) <- c("values", "low_CI", "high_CI")

write.csv(t_tune_initial,
          paste0(dir_name, "tune_initial_with_CI.csv"), row.names = T)

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


################################################################################
# additional analysis for pmcmc samples
pmcmc_samples <- readRDS(paste0(dir_name, "pmcmc_samples.rds"))
data <- readRDS("inputs/pmcmc_data_week_allAge.rds")

# gen_sir$new(pars = list(), time = 0, n_particles = 1L)$info()
plot_states <- function(state, data) {
  col <- grey(0.3, 0.1)
  matplot(data$yearWeek, t(state["cases_55", , -1]),
          type = "l", lty = 1, col = col,
          xlab = "", ylab = "GPSC55 cases",
          xlim = as.Date(c(min(data$yearWeek), max(data$yearWeek))),
          ylim = c(0, max(data$count_WGS_GPSC55, na.rm = TRUE))
          )
  lines(data$yearWeek, data$count_WGS_GPSC55, col = 3, lty = 1, lwd = 2)

  matplot(data$yearWeek, t(state["cases_12F", , -1]),
          type = "l", lty = 1, col = col,
          xlab = "", ylab = "12F cases",
          xlim = as.Date(c(min(data$yearWeek), max(data$yearWeek))),
          ylim = c(0, max(data$count_serotype, na.rm = TRUE))
          )
  lines(data$yearWeek, data$count_serotype, col = 2, lty = 1, lwd = 2)
  lines(data$yearWeek, data$count_WGS_GPSC55 + data$count_WGS_non55,
         col = 3, lty = 1, lwd = 2)

  matplot(data$yearWeek, t(state["cases_non55", , -1]),
          type = "l", lty = 1, col = col,
          xlab = "", ylab = "non55 cases",
          ylim = c(0, max(data$count_WGS_non55, na.rm = TRUE))
          )
  lines(data$yearWeek, data$count_WGS_non55, col = 2,lty = 1, lwd = 2)

  matplot(data$yearWeek, t(state["n_AD_weekly", , -1]),
          type = "l", lty = 1, col = col,
          xlab = "", ylab = "Diseased cases (Infections)")

  matplot(data$yearWeek, t(state["Ne", , -1]),
          type = "l", lty = 1, col = col,
          xlab = "", ylab = "Ne")
  lines(data$yearWeek, data$Ne, col = 2, lty = 1, lwd = 2)

  # matplot(data$yearWeek, xlab = "", t(state["S_tot", , -1]),
  #         type = "l", lty = 1, col = 2, ylab = "%", ylim = c(0, 6e7), yaxt = "n")
  # axis(side = 2, at = seq(0, 6e7, length.out = 5),
  #      labels = seq(0, 100, length.out = 5))

  # matlines(data$yearWeek, t(state["I_tot", , -1]), lty = 1, col = 3)
  # matlines(data$yearWeek, t(state["R_tot", , -1]), lty = 1, col = 4)
  # legend("left", bty = "n", fill = 2:4, legend = c("S", "I", "R"))
  # 
  # matplot(data$yearWeek, xlab = "", t(state["D_tot", , -1]),
  #         type = "l", lty = 1, col = 3, ylab = "carriers")
  
  
}














png("figs/traces.png", width = 17, height = 12, unit = "cm", res = 1200)
par(mfrow = c(3, 4), mar = c(3, 3, 1, 1), mgp = c(1.7, 0.7, 0), bty = "n")
plot_traces(pmcmc_results$pars, metrics = pmcmc_results$metrics)
plot_traces(pmcmc_results$probabilities)
dev.off()


png("figs/posteriors.png", width = 17, height = 8, unit = "cm", res = 1200)
par(mfrow = c(3, 3), mar = c(3, 3, 1, 1), mgp = c(1.7, 0.7, 0), bty = "n")
for (nm in names(initial_pars)) {
  hist(pmcmc_samples$pars[, nm], xlab = nm, main = "", freq = FALSE)
  if (nm %in% names(priors)) {
    curve(exp(priors[[nm]](x)), add = TRUE, col = 2)
  }
}
dev.off()


# fixed pics
png(paste0(dir_name, "figs/pmcmc_samples_fits.png"),
    width = 24, height = 17, unit = "cm", res = 600)
par(mfrow = c(3, 2), bty = "n", mar = c(3, 3, 1, 1), mgp = c(1.7, 0.7, 0))
plot_states(pmcmc_samples$trajectories$state, data)
dev.off()

png(paste0(dir_name, "figs/pmcmc_samples_diag_ggPairs_%02d.png"),
    width = 17, height = 17, unit = "cm", res = 600)
par(mar = c(3, 3, 1, 1), mgp = c(1.7, 0.7, 0))
GGally::ggpairs(as.data.frame(pmcmc_samples$pars))
dev.off()




































