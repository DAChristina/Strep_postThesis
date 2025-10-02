# A_ini converter
N <- 6.7e7 # FIXED England's pop size is roughly 67,000,000
max_A_ini <- 0
min_A_ini <- -10
# scaled_A_ini <- user(0)
# log_A_ini <- scaled_A_ini*(max_A_ini-min_A_ini)+min_A_ini # scaled_A_ini*(max_A_ini−min_A_ini)+min_A_ini; rescaled using (A_ini-A_ini_min)/(A_ini_max-A_ini_min)

# directly test log_A_ini as scaled
dir_name <- paste0("outputs/genomics/trial_", "5e+05_Aini_dgamma_shape3_scale0.03", "/")
results <- read.csv(paste0(dir_name, "tune_initial_with_CI.csv")) %>% 
  glimpse()

log_A_ini <- c(results[1,2], results[1,3], results[1,4])
rescaled_log_A_ini <- log_A_ini*(max_A_ini-min_A_ini)+min_A_ini
rescaled_log_A_ini

A_ini <- 10^(rescaled_log_A_ini)*N
cat("Initial asymptomatic count:","\n",
    round(A_ini[1], 2), " (",
    round(A_ini[2], 2), "-",
    round(A_ini[3], 2), ")"
    )


# compare 2 pmcmc simulations
library(tidyverse)
scale_3e <- readRDS("outputs/genomics/trial_5e+05_Aini_dgamma_shape3_scale0.03/pmcmc_samples.rds") %>% 
  glimpse()
scale_5e <- readRDS("outputs/genomics/trial_5e+05_Aini_dgamma_shape3_scale0.05/pmcmc_samples.rds") %>% 
  glimpse()

# chain comparison
data.frame(
  chain = c(scale_3e$chain, scale_5e$chain),
  gamma_scale = rep(c("Scale 0.03", "Scale 0.05"),
              c(length(scale_3e$chain), length(scale_5e$chain)))
) %>% 
  ggplot(., aes(x = chain, fill = gamma_scale)) +
  geom_density(alpha = 0.4) +
  labs(title = "Marginal distributions (500,000 steps with different gamma's Scale)") +
  theme_bw()
  

# A_ini comparison
dplyr::bind_rows(
  as.data.frame(scale_3e$pars) %>% 
    dplyr::mutate(
      chain = scale_3e$chain,
      iteration = scale_3e$iteration,
      gamma_scale = "Scale 0.03"
    )
  ,
  as.data.frame(scale_5e$pars) %>% 
    dplyr::mutate(
      chain = scale_5e$chain,
      iteration = scale_5e$iteration,
      gamma_scale = "Scale 0.05"
    )
) %>% 
  ggplot(., aes(x = log_A_ini, fill = gamma_scale)) +
  geom_density(alpha = 0.4) +
  labs(title = "Trace plot of scaled log(A_ini)") +
  theme_bw()


# all parameters combined in one plot
dplyr::bind_rows(
  as.data.frame(scale_3e$pars) %>% 
    dplyr::mutate(
      gamma_scale = "Scale 0.03"
    )
  ,
  as.data.frame(scale_5e$pars) %>% 
    dplyr::mutate(
      gamma_scale = "Scale 0.05"
    )
) %>% 
  tidyr::pivot_longer(
    cols = c("log_A_ini", "time_shift_1", "beta_0", "beta_1",
             "log_delta1", "log_delta2", "kappa_55"),
    names_to = "parameters"
  ) %>% 
  ggplot(., aes(x = value, fill = gamma_scale)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ parameters,
             scales = "free") +
  theme_bw()

# probabilities comparison
prob_1 <- data.frame(
  log_prior = c(scale_3e$probabilities[, "log_prior"],
                scale_5e$probabilities[, "log_prior"]
  ),
  gamma_scale = rep(c("Scale 0.03", "Scale 0.05"),
                    c(length(scale_3e$chain), length(scale_5e$chain)))
) %>% 
  ggplot(., aes(x = log_prior, fill = gamma_scale)) +
  geom_density(alpha = 0.4) +
  labs(title = "Marginal distribution of Log Prior\n(500,000 steps with different gamma's Scale)") +
  theme_bw()+
  theme(legend.position = "none")

prob_2 <- data.frame(
  log_lik   = c(scale_3e$probabilities[, "log_likelihood"],
                scale_5e$probabilities[, "log_likelihood"]
  ),
  gamma_scale = rep(c("Scale 0.03", "Scale 0.05"),
                    c(length(scale_3e$chain), length(scale_5e$chain)))
) %>% 
  ggplot(., aes(x = log_lik, fill = gamma_scale)) +
  geom_density(alpha = 0.4) +
  labs(title = "Marginal distribution of Log Likelihood\n(500,000 steps with different gamma's Scale)") +
  theme_bw()+
  theme(legend.position = "none")

prob_3 <- data.frame(
  log_post  = c(scale_3e$probabilities[, "log_posterior"],
                scale_5e$probabilities[, "log_posterior"]
  ),
  gamma_scale = rep(c("Scale 0.03", "Scale 0.05"),
                    c(length(scale_3e$chain), length(scale_5e$chain)))
) %>% 
  ggplot(., aes(x = log_post, fill = gamma_scale)) +
  geom_density(alpha = 0.4) +
  labs(title = "Marginal distribution of Log Posterior\n(500,000 steps with different gamma's Scale)") +
  theme_bw()+
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(size = 9),
        legend.text  = element_text(size = 9))

cowplot::plot_grid(prob_1, prob_2, prob_3,
                   nrow = 3,
                   labels = c("A", "B", "C"))

