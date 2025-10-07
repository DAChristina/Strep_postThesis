source("global/all_function_allAge.R")
source("model/sir_stochastic_allAge_2post_pmcmc_picts.R")

source("R/3_pmcmc.R")
source("R/4_post_pmcmc_pics.R")
source("R/5_post_pmcmc_samples_pics.R")
source("R/6_post_pmcmc_age_validation.R")

pmcmc_run_plus_tuning(n_pars = 10, n_sts = 600,
                      run1_stochastic = F, run2_stochastic = F, ncpus = 4)
post_pmcmc_pics(600)
model_vs_data(600)
post_particle_pics(600)
age_validation(600)

pmcmc_run_plus_tuning(n_pars = 10, n_sts = 2000,
                      run1_stochastic = F, run2_stochastic = F, ncpus = 4)
post_pmcmc_pics(2000)
model_vs_data(2000)
post_particle_pics(2000)
age_validation(2000)






# pmcmc_run2_only(n_pars = 10, n_sts = 5250,
#                 run2_stochastic = F, ncpus = 4)



post_pmcmc_pics(5020)
post_pmcmc_pics(5030)

post_particle_pics(5020)
post_particle_pics(5030)

model_vs_data(5020)
model_vs_data(5030)

model_vs_data(5250)
model_vs_data(10250)
model_vs_data(100250)

post_particle_pics(10250)



dir_name <- paste0("outputs/genomics/trial_", 600, "/")
dir.create(paste0(dir_name, "/figs"), FALSE, TRUE)
# run 4_post_pmcmc_pics.R first
results <- read.csv(paste0(dir_name, "tune_initial_with_CI.csv")) %>% 
  glimpse()


# gen_sir <- odin.dust::odin_dust("model/sir_basic_trial.R")
gen_sir <- odin.dust::odin_dust("model/sir_stochastic_allAge.R")

# Running the SIR model with dust
pars <- list(N_ini = 6.7e7,
             log_A_ini = results[1,2],
             D_ini = 0,
             R_ini = 0,
             time_shift_1 = results[2,2],
             # time_shift_2 = results[3,2],
             beta_0 = results[3,2],
             beta_1 = results[4,2],
             # beta_2 = results[6,2],
             # scaled_wane = results[7,2],
             # psi = (0.5),
             hypo_sigma_2 = (1),
             log_delta = results[5,2]
             # alpha = results[6,2],
             # gamma_annual = results[7,2],
             # nu_annual = results[8,2]
)

time_points <- round(seq(0, by = (365/52), length.out = 52*3)) # per-week, 22 years
# n_times <- length(time_points)
n_times <- 13800 # roughly from 1987-2025 in days
sir_model <- gen_sir$new(pars = pars,
                         time = 1,
                         n_particles = 1L,
                         n_threads = 4L,
                         seed = 1L)

# compartment position check
# sir_model$info()
# sir_model$info()$index$n_AD_weekly
# update_state is required "every single time" to run & produce matrix output (don't know why)
# sir_model$update_state(pars = pars,
#                        time = 0) # make sure time is 0

# all_date <- incidence$day
# all_date <- data.frame(col = integer(4745))
# incidence <- read.csv("inputs/incidence_week_12F_allAge.csv") %>% 
#   dplyr::mutate(day = week*7) 
n_particles <- 1L
model <- array(NA, dim = c(sir_model$info()$len, n_particles, n_times))

# R0 estimation
R0 <- pars$beta_0/pars$sigma_2
R0

for (t in seq_len(n_times)) {
  model[ , , t] <- sir_model$run(t)
}
# time <- x[1, 1, ] # because in the position of [1, 1, ] is time
# x <- x[-1, , ] # compile all matrix into 1 huge df, delete time (position [-1, , ])

data <- readRDS("inputs/pmcmc_data_week_allAge_nonGAM.rds") %>% 
  glimpse()

sir_data <- data %>% 
  dplyr::transmute(
    replicate = 1,
    # steps = time_start+1,
    weekly = seq_along(replicate),
    value = count_WGS_GPSC55,
    compartment = "data_count_WGS_GPSC55"
  ) %>%
  glimpse()

# all_dates <- data.frame(date = seq(min(data$yearWeek), max(data$yearWeek), by = "day")) %>%
#   dplyr::mutate(
#     steps = seq_along(date)
#   ) %>%
#   glimpse()
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
                                   index == 3 ~ "D",
                                   index == 4 ~ "S",
                                   index == 5 ~ "R",
                                   index == 6 ~ "model_n_AD_weekly",
                                   index == 7 ~ "Ne",
                                   index == 8 ~ "cases_55",
                                   index == 9 ~ "cases_non55",
                                   index == 10 ~ "cases_12F"
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
        value = count_serotype,
        compartment = "data_count_12F"
      )
  ) %>% 
  dplyr::bind_rows(
    data %>% 
      dplyr::transmute(
        replicate = 1,
        weekly = seq_along(replicate),
        value = count_serotype,
        compartment = "data_count_WGS_non55"
      )
  ) %>% 
  dplyr::bind_rows(
    data %>% 
      dplyr::transmute(
        replicate = 1,
        weekly = seq_along(replicate),
        value = count_serotype,
        compartment = "data_Ne"
      )
  ) %>% 
  dplyr::full_join(
    all_dates
    ,
    by = "weekly"
  ) %>%
  glimpse()


ser12 <- ggplot(incidence_modelled %>% 
                  dplyr::filter(
                    compartment %in% c("data_count_WGS_GPSC55", "data_count_12F",
                                       "cases_12F"),
                    compartment != "Time"
                  )
                ,
                aes(x = yearWeek, y = value,
                    group = interaction(compartment,replicate),
                    colour = compartment)) +
  geom_line() +
  # scale_y_continuous(trans = "log1p") +
  # scale_y_continuous(limits = c(0, 50)) +
  # scale_x_continuous(limits = c(0, 700)) +
  scale_x_date(limits = c(as.Date(min(all_dates$yearWeek)), as.Date(max(all_dates$yearWeek))),
               date_breaks = "year",
               date_labels = "%Y") +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))
  
gps55 <- ggplot(incidence_modelled %>% 
                  dplyr::filter(
                    compartment %in% c("data_count_WGS_GPSC55", "data_count_12F",
                                       "cases_55"),
                    compartment != "Time"
                  )
                ,
                aes(x = yearWeek, y = value,
                    group = interaction(compartment,replicate),
                    colour = compartment)) +
  geom_line() +
  # scale_y_continuous(trans = "log1p") +
  # scale_y_continuous(limits = c(0, 50)) +
  # scale_x_continuous(limits = c(0, 700)) +
  scale_x_date(limits = c(as.Date(min(all_dates$yearWeek)), as.Date(max(all_dates$yearWeek))),
               date_breaks = "year",
               date_labels = "%Y") +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))

non55 <- ggplot(incidence_modelled %>% 
                  dplyr::filter(
                    compartment %in% c("data_count_WGS_GPSC55", "data_count_12F",
                                       "cases_non55"),
                    compartment != "Time"
                  )
                ,
                aes(x = yearWeek, y = value,
                    group = interaction(compartment,replicate),
                    colour = compartment)) +
  geom_line() +
  # scale_y_continuous(trans = "log1p") +
  # scale_y_continuous(limits = c(0, 50)) +
  # scale_x_continuous(limits = c(0, 700)) +
  scale_x_date(limits = c(as.Date(min(all_dates$yearWeek)), as.Date(max(all_dates$yearWeek))),
               date_breaks = "year",
               date_labels = "%Y") +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))

Ne <- ggplot(incidence_modelled %>% 
               dplyr::filter(
                 compartment %in% c("data_count_WGS_GPSC55", "data_count_12F",
                                    "Ne", "data_Ne"),
                 compartment != "Time"
               )
             ,
             aes(x = yearWeek, y = value,
                 group = interaction(compartment,replicate),
                 colour = compartment)) +
  geom_line() +
  # scale_y_continuous(trans = "log1p") +
  # scale_y_continuous(limits = c(0, 50)) +
  # scale_x_continuous(limits = c(0, 700)) +
  scale_x_date(limits = c(as.Date(min(all_dates$yearWeek)), as.Date(max(all_dates$yearWeek))),
               date_breaks = "year",
               date_labels = "%Y") +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))

png(paste0(dir_name, "figs/model_vs_data3.png"),
    width = 24, height = 24, unit = "cm", res = 600)
p <- cowplot::plot_grid(ser12, gps55, non55, Ne,
                        nrow = 4,
                        labels = c("A", "B", "C", "D"))
print(p)
dev.off()


ggplot(incidence_modelled %>% 
              dplyr::filter(
                grepl("cases|D|data", compartment),
                # compartment %in% c("D", "model_n_AD_weekly", "data_count_WGS_GPSC55"), # redesign the model, would rather fit to D
                # compartment %in% c("model_n_AD_weekly", "data_count_WGS_GPSC55"),
                # compartment %in% c("D", "n_AD_weekly"),
                # compartment %in% c("D", "data_count_WGS_GPSC55", "data_count_12F"),
                # grepl("cases|data", compartment),
                # compartment %in% c("D"),
                compartment != "Time",
                # compartment %in% c("S")
              )
            ,
            aes(x = yearWeek, y = value,
                group = interaction(compartment,replicate),
                colour = compartment)) +
  geom_line() +
  # scale_y_continuous(trans = "log1p") +
  # scale_y_continuous(limits = c(0, 50)) +
  # scale_x_continuous(limits = c(0, 700)) +
  scale_x_date(limits = c(as.Date(min(all_dates$yearWeek)), as.Date(max(all_dates$yearWeek))),
               date_breaks = "year",
               date_labels = "%Y") +
  ggtitle("Cases (Aggregated by Week)") +
  xlab("Time") +
  ylab("Number of People") +
  theme_bw() +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))
