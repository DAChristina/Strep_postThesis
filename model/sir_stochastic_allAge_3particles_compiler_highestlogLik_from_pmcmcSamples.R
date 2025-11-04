# highest likelihood from pmcmc_samples

library(tidyverse)
library(odin.dust)

# I update odin.dust by force
# remotes::install_github("mrc-ide/odin.dust")
source("global/all_function_allAge.R")

dir_name <- paste0("outputs/genomics/trial_", "5e+05_4with_mcmc_limitations_adjusted_beta0_25_0.01_FINAL", "/")
dir.create(paste0(dir_name, "/figs"), FALSE, TRUE)

pmcmc_samples <- readRDS(paste0(dir_name,
                                "pmcmc_samples.rds"))

top10_idx <- order(pmcmc_samples$probabilities[, "log_likelihood"],
                   decreasing = TRUE)[1:10]
top10_pars <- pmcmc_samples$pars[top10_idx, ]

# Combine with their log posterior values for clarity
top10_final <- cbind(particle = paste0("particle_", as.numeric(1:10)),
                     top10_pars, 
                     log_prior = pmcmc_samples$probabilities[top10_idx,
                                                             "log_prior"],
                     log_likelihood = pmcmc_samples$probabilities[top10_idx,
                                                                  "log_likelihood"],
                     log_posterior = pmcmc_samples$probabilities[top10_idx,
                                                                 "log_posterior"]) %>% 
  glimpse()

# write.csv(top10_final,
#           paste0(dir_name, "tune_top10_particles.csv"), row.names = F)

t(top10_final)

results <- read.csv(paste0(dir_name, "tune_initial_with_CI.csv")) %>% 
  dplyr::right_join(
    t(top10_final) %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column("X")
    ,
    by = "X"
  ) %>% 
  dplyr::transmute(parameters = X,
                   low_CI = as.numeric(low_CI),
                   high_CI = as.numeric(high_CI),
                   parameter_values = as.numeric(values),
                   particle_1  = as.numeric(V1),
                   particle_2  = as.numeric(V2),
                   particle_3  = as.numeric(V3),
                   particle_4  = as.numeric(V4),
                   particle_5  = as.numeric(V5),
                   particle_6  = as.numeric(V6),
                   particle_7  = as.numeric(V7),
                   particle_8  = as.numeric(V8),
                   particle_9  = as.numeric(V9),
                   particle_10 = as.numeric(V10)
  ) %>% 
  dplyr::filter(parameters != "particle") %>% 
  glimpse()

# write.csv(results,
#           paste0(dir_name, "tune_initial_with_CI_combined_top10_particles.csv"),
#           row.names = F)


gen_sir <- odin.dust::odin_dust("model/sir_stochastic_allAge.R")

# Running the SIR model with dust
pars <- list(
  N_ini = 6.7e7,
  log_A_ini = as.numeric(results[1, 5]),
  D_ini = 0,
  R_ini = 0,
  time_shift_1 = as.numeric(results[2, 5]),
  beta_0 = as.numeric(results[3, 5]),
  beta_1 = as.numeric(results[4, 5]),
  hypo_sigma_2 = 1,
  log_delta1 = as.numeric(results[5, 5]),
  log_delta2 = as.numeric(results[6, 5])
)

time_points <- round(seq(0, by = (365/52), length.out = 52*3)) # per-week, 22 years
# n_times <- length(time_points)
n_times <- 13800 # roughly from 1987-2025 in days
sir_model <- gen_sir$new(pars = pars,
                         time = 1,
                         n_particles = 1L,
                         n_threads = 4L,
                         seed = 1L)


n_particles <- 1L
model <- array(NA, dim = c(sir_model$info()$len, n_particles, n_times))

for (t in seq_len(n_times)) {
  model[ , , t] <- sir_model$run(t)
}

data <- readRDS("inputs/pmcmc_data_week_ageGroup12F.rds") %>% 
  glimpse()

sir_data <- dplyr::bind_rows(
  data %>% 
    dplyr::transmute(
      replicate = 1,
      # steps = time_start+1,
      weekly = seq_along(replicate),
      value = count_55_all,
      compartment = "data_count_55_all"
    )
  ,
  data %>% 
    dplyr::transmute(
      replicate = 1,
      # steps = time_start+1,
      weekly = seq_along(replicate),
      value = count_55_1,
      compartment = "data_count_55_1"
    )
  ,
  data %>% 
    dplyr::transmute(
      replicate = 1,
      # steps = time_start+1,
      weekly = seq_along(replicate),
      value = count_55_2,
      compartment = "data_count_55_2"
    )
) %>%
  glimpse()

all_dates <- data %>%
  dplyr::select(yearWeek) %>% 
  dplyr::mutate(
    weekly = seq_along(yearWeek)
  ) %>%
  glimpse()

# focused on n_AD_weekly (already in weeks)
incidence_modelled <- 
  reshape2::melt(model) %>% 
  dplyr::rename(index = Var1,     # Var1 = dimension that stored SADR values
                replicate = Var2, # Var2 = particles
                steps = Var3       # Var3 = steps are in days, but n_AD_weekly is aggregated in weeks
  ) %>% 
  # dplyr::filter(index < 5) %>%
  dplyr::mutate(compartment = 
                  dplyr::case_when(index == 1 ~ "Time",
                                   index == 2 ~ "A",
                                   index == 3 ~ "model_D1",
                                   index == 4 ~ "model_D2",
                                   index == 5 ~ "model_D",
                                   index == 6 ~ "S",
                                   index == 7 ~ "R",
                                   index == 8 ~ "model_n_AD1_weekly",
                                   index == 9 ~ "model_n_AD2_weekly"
                  )) %>% 
  dplyr::select(-index) %>%
  dplyr::mutate(weekly = ceiling(steps/7)) %>% 
  dplyr::group_by(replicate, weekly, compartment) %>% 
  dplyr::summarise(value = sum(value, na.rm = T),
                   # date = max(date),
                   .groups = "drop") %>% 
  dplyr::ungroup() %>% 
  dplyr::bind_rows(sir_data) %>%
  # add 12F data for comparison
  dplyr::bind_rows(
    data %>% 
      dplyr::transmute(
        replicate = 1,
        weekly = seq_along(replicate),
        value = count_12F_1 + count_12F_2,
        compartment = "data_count_12F_all"
      )
  ) %>% 
  dplyr::bind_rows(
    data %>% 
      dplyr::transmute(
        replicate = 1,
        weekly = seq_along(replicate),
        value = count_12F_1,
        compartment = "data_count_12F_1"
      )
  ) %>% 
  dplyr::bind_rows(
    data %>% 
      dplyr::transmute(
        replicate = 1,
        weekly = seq_along(replicate),
        value = count_12F_2,
        compartment = "data_count_12F_2"
      )
  ) %>% 
  dplyr::full_join(
    all_dates
    ,
    by = "weekly"
  ) %>%
  glimpse()

pdat <- incidence_modelled %>%
  dplyr::filter(compartment %in% c("model_D",
                                   "model_D1",
                                   "model_D2",
                                   "data_count_55_1",
                                   "data_count_55_2",
                                   "data_count_55_all",
                                   "data_count_12F_1",
                                   "data_count_12F_2",
                                   "data_count_12F_all")) %>%
  tidyr::pivot_wider(
    id_cols = yearWeek,
    names_from = compartment,
    values_from = value
  ) %>% 
  tidyr::unnest(cols = c(model_D, model_D1, model_D2,
                         data_count_55_1,
                         data_count_55_2,
                         data_count_55_all,
                         data_count_12F_1,
                         data_count_12F_2,
                         data_count_12F_all)) %>% 
  dplyr::mutate(kappa = results[7,5],
                model_D_low = qnbinom(0.025,
                                      size = kappa,
                                      mu = model_D),
                model_D_high = qnbinom(0.975,
                                       size = kappa,
                                       mu = model_D),
                model_D1_low = qnbinom(0.025,
                                      size = kappa,
                                      mu = model_D1),
                model_D1_high = qnbinom(0.975,
                                       size = kappa,
                                       mu = model_D1),
                model_D2_low = qnbinom(0.025,
                                      size = kappa,
                                      mu = model_D2),
                model_D2_high = qnbinom(0.975,
                                       size = kappa,
                                       mu = model_D2),
                
  ) %>% 
  tidyr::complete(yearWeek,
                  fill = list(data_count_55_all = 0,
                              data_count_55_1 = 0,
                              data_count_55_2 = 0,
                              data_count_12F_all = 0,
                              data_count_12F_1 = 0,
                              data_count_12F_2 = 0)) %>% 
  glimpse()



ggplot(pdat, aes(x = yearWeek)) +
  geom_ribbon(aes(ymin = model_D_low,
                  ymax = model_D_high), 
              fill = "skyblue", alpha = 0.6) +
  geom_line(aes(y = model_D, colour = "model_D"), 
            size = 0.4) +
  geom_line(aes(y = data_count_12F_all, colour = "data_count_12F_all"), 
            size = 0.4) +
  geom_line(aes(y = data_count_55_all, colour = "data_count_55_all"), 
            size = 0.4) +
  
  # 10 particle lines
  # geom_line(aes(y = particle_value, group = particle, colour = "particles"),
  #           alpha = 0.4, size = 0.2) +
  scale_x_date(
    date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title = "Model vs. Data with 95% CrI (Negative Binomial Observation Distribution)",
    x = "Time",
    y = "Number of People",
    colour = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.1, 0.85),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key.size = unit(0.8, "lines"),
    legend.text = element_text(size = 10)
  )

p_D1 <- ggplot(pdat, aes(x = yearWeek)) +
  geom_ribbon(aes(ymin = model_D1_low,
                  ymax = model_D1_high), 
              fill = "skyblue", alpha = 0.6) +
  geom_line(aes(y = model_D1, colour = "model_D1"), 
            size = 0.4) +
  geom_line(aes(y = data_count_12F_1, colour = "data_count_12F_1"), 
            size = 0.4) +
  geom_line(aes(y = data_count_55_1, colour = "data_count_55_1"), 
            size = 0.4) +
  scale_x_date(
    date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title = "Model vs. Data (age group 0-44)\nwith 95% CrI (Negative Binomial Observation Distribution)",
    x = "Time",
    y = "Number of People",
    colour = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.1, 0.85),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key.size = unit(0.8, "lines"),
    legend.text = element_text(size = 10)
  )
p_D1

p_D2 <- ggplot(pdat, aes(x = yearWeek)) +
  geom_ribbon(aes(ymin = model_D2_low,
                  ymax = model_D2_high), 
              fill = "skyblue", alpha = 0.6) +
  geom_line(aes(y = model_D2, colour = "model_D2"), 
            size = 0.4) +
  geom_line(aes(y = data_count_12F_2, colour = "data_count_12F_2"), 
            size = 0.4) +
  geom_line(aes(y = data_count_55_2, colour = "data_count_55_2"), 
            size = 0.4) +
  scale_x_date(
    date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title = "Model vs. Data (age group 45+)\nwith 95% CrI (Negative Binomial Observation Distribution)",
    x = "Time",
    y = "Number of People",
    colour = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.1, 0.85),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key.size = unit(0.8, "lines"),
    legend.text = element_text(size = 10)
  )
p_D2


# png(paste0(dir_name, "figs/model_vs_data_withCrI_and10Particles_separated_2ageGroups.png"),
#     width = 24, height = 24, unit = "cm", res = 600)
cowplot::plot_grid(p_D1, p_D2,
                   nrow = 2,
                   labels = c("A", "B")
)
# dev.off()
