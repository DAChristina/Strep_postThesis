library(tidyverse)
library(odin.dust)

# I update odin.dust by force
# remotes::install_github("mrc-ide/odin.dust")
source("global/all_function_allAge.R")

dir_name <- paste0("outputs/genomics/trial_", "5e+05_serotype1_FINAL", "/")
dir.create(paste0(dir_name, "/figs"), FALSE, TRUE)

pmcmc_samples <- readRDS(paste0(dir_name,
                                "pmcmc_samples.rds"))

top10_idx <- order(pmcmc_samples$probabilities[, "log_posterior"],
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

write.csv(top10_final,
          paste0(dir_name, "tune_top10_particles.csv"), row.names = F)

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

write.csv(results,
          paste0(dir_name, "tune_initial_with_CI_combined_top10_particles.csv"),
          row.names = F)


# gen_sir <- odin.dust::odin_dust("model/sir_basic_trial.R")
gen_sir <- odin.dust::odin_dust("model/sir_stochastic_allAge.R")

data_storage <- data.frame()
particles <- colnames(results)[4:ncol(results)]

sir_model <- list()
model <- list()
pars <- list()
n_particles <- 1L
n_times <- 7669 # roughly from 2003 to 2023 in days

data <- readRDS("inputs/pmcmc_data_week_allAge_ser1.rds") %>% 
  glimpse()

sir_data <- data %>% 
  dplyr::transmute(
    replicate = 1,
    weekly = seq_along(replicate),
    value = count_serotype,
    compartment = "data_count_serotype"
  ) %>%
  glimpse()

all_dates <- data %>%
  dplyr::select(yearWeek) %>% 
  dplyr::mutate(
    weekly = seq_along(yearWeek)
  ) %>%
  glimpse()

for (i in particles){
  # Running the SIR model with dust
  pars[[i]] <- list(
    N_ini = 6.7e7,
    log_A_ini = as.numeric(results[1, i]),
    D_ini = 0,
    R_ini = 0,
    time_shift_1 = as.numeric(results[2, i]),
    beta_0 = as.numeric(results[3, i]),
    beta_1 = as.numeric(results[4, i]),
    hypo_sigma_2 = 1,
    log_delta = as.numeric(results[5, i])
  )
  
  sir_model[[i]] <- gen_sir$new(pars = pars[[i]],
                                time = 1,
                                n_particles = n_particles,
                                n_threads = 4L,
                                seed = 1L)
  
  # compartment position check
  # sir_model$info()
  
  model[[i]] <- array(NA,
                      dim = c(sir_model[[i]]$info()$len,
                              n_particles, n_times))
  particle_storage <- i
  
  for (t in seq_len(n_times)) {
    model[[i]][ , , t] <- sir_model[[i]]$run(t)
  }
  
  # focused on n_AD_weekly (already in weeks)
  incidence_modelled <- 
    reshape2::melt(model[[i]]) %>% 
    dplyr::rename(index = Var1,     # Var1 = dimension that stored SADR values
                  replicate = Var2, # Var2 = particles
                  steps = Var3       # Var3 = steps are in days, but n_AD_weekly is aggregated in weeks
    ) %>% 
    # dplyr::filter(index < 5) %>%
    dplyr::mutate(compartment = 
                    dplyr::case_when(index == 1 ~ "Time",
                                     index == 2 ~ "A",
                                     index == 3 ~ "D",
                                     index == 4 ~ "S",
                                     index == 5 ~ "R",
                                     index == 6 ~ "model_n_AD_weekly",
                                     index == 7 ~ "Ne",
                                     index == 8 ~ "cases_55",
                                     index == 9 ~ "cases_non55",
                                     index == 10 ~ "cases_12F"
                    ),
                  particles = particle_storage,
                  weekly = ceiling(steps / 7)
                  ) %>% 
    dplyr::select(-index) %>%
    dplyr::group_by(particles, replicate, weekly, compartment) %>% 
    dplyr::summarise(value = sum(value, na.rm = T),
                     # date = max(date),
                     .groups = "drop") %>% 
    dplyr::ungroup() %>% 
    dplyr::bind_rows(sir_data) %>%
    dplyr::full_join(
      all_dates
      ,
      by = "weekly"
    ) %>%
    glimpse()
  
  # write.csv(incidence_modelled,
  #           paste0(dir_name, "incidence_modelled_serotype1.csv"),
  #           row.names = FALSE)
  
  data_storage <- dplyr::bind_rows(data_storage,
                                   incidence_modelled) %>% 
    dplyr::mutate(particles = ifelse(is.na(particles), compartment, particles))
}

glimpse(data_storage)

# split 2 data types for particles and CrI
pdat <- data_storage %>%
  dplyr::filter(compartment %in% c("D", "data_count_serotype")) %>%
  tidyr::pivot_wider(
    id_cols = yearWeek,
    names_from = particles,
    values_from = value
  ) %>% 
  tidyr::unnest(cols = c(parameter_values,
                         data_count_serotype,
                         particle_1, particle_2, particle_3,
                         particle_4, particle_5, particle_6,
                         particle_7, particle_8, particle_9,
                         particle_10)) %>%
  dplyr::rename(model_D = parameter_values) %>%
  dplyr::mutate(model_D_low = qnbinom(0.025,
                                      size = results[6,4],
                                      mu = model_D),
                model_D_high = qnbinom(0.975,
                                       size = results[6,4],
                                       mu = model_D)
                ) %>% 
  tidyr::pivot_longer(
    cols = starts_with("particle_"),
    names_to = "particle",
    values_to = "particle_value"
  ) %>%
  tidyr::complete(yearWeek,
                  fill = list(data_count_serotype = 0)) %>% 
  glimpse()

write.csv(pdat,
          paste0(dir_name,
                 "incidence_modelled_serotype1_with_top10_particles.csv"),
          row.names = F)

# separate into 2 plots
# separate into 2 plots
p1 <- ggplot(pdat, aes(x = yearWeek)) +
  geom_line(aes(y = model_D, colour = "model_D"), 
            size = 0.3) +
  geom_line(aes(y = data_count_serotype, colour = "data_count_serotype"), 
            size = 0.3) +
  
  # 10 particle lines
  # geom_line(aes(y = particle_value, group = particle, colour = "particles"),
  #           alpha = 0.4, size = 0.2) +
  geom_ribbon(aes(ymin = model_D_low,
                  ymax = model_D_high), 
              fill = "skyblue", alpha = 0.5) +
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
p1

p2 <- ggplot(pdat, aes(x = yearWeek)) +
  geom_line(aes(y = model_D, colour = "model_D"), 
            size = 0.3) +
  geom_line(aes(y = data_count_serotype, colour = "data_count_serotype"), 
            size = 0.3) +
  
  # 10 particle lines
  geom_line(aes(y = particle_value, group = particle, colour = "particles"),
            alpha = 0.4, size = 0.2) +
  # geom_ribbon(aes(ymin = model_D_low,
  #                 ymax = model_D_high), 
  #             fill = "skyblue", alpha = 0.5) +
  scale_x_date(
    date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title = "Model vs. Data with Top 10 Particles",
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
p2

png(paste0(dir_name, "figs/model_vs_data_withCrI_and10Particles.png"),
    width = 24, height = 24, unit = "cm", res = 600)
cowplot::plot_grid(p1, p2,
                   nrow = 2,
                   labels = c("A", "B")
)

dev.off()

