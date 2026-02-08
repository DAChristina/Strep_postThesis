library(tidyverse)
library(odin.dust)
library(socialmixr)

# I update odin.dust by force
# remotes::install_github("mrc-ide/odin.dust")

gen_sir <- odin.dust::odin_dust("model/sir_basic_trial.R")

# Create contact_matrix 5 demographic groups:
# > 5
# 5-18
# 19-30
# 31-64
# 65+
# age.limits = c(0, 5, 19, 31, 65)

# Create contact_matrix 2 demographic groups:
# < 15
# 15+
age.limits = c(0, 15)
N_age <- length(age.limits)

contact_2_demographic <- suppressMessages(
  socialmixr::contact_matrix(polymod,
                             countries = "United Kingdom",
                             age.limits = age.limits,
                             symmetric = TRUE
  ))

transmission <- contact_2_demographic$matrix /
  rep(contact_2_demographic$demography$population,
      each = ncol(contact_2_demographic$matrix))
t_norm <- transmission/max(transmission)

pars <- list(m = t_norm,
             N_ini = contact_2_demographic$demography$population,
             log_E_ini = c(0.5, 0.5), # test c(0.65, 0.35),
             I_ini = c(3,3),
             log_delta1 = -5,
             log_delta2 = -5,
             hypo_sigma_1_day = 16,
             time_shift_1 = 0.1254,
             beta_0 = 0.5,
             beta_1 = 0.3
)

n_times <- 7500*24 # 500 for trial
n_pars <- 1L
sir_model <- gen_sir$new(pars = pars,
                         time = 1,
                         n_particles = n_pars,
                         n_threads = 4L,
                         seed = 1L)

# compartment position check
# sir_model$info()
# sir_model$info()$index$n_AD_weekly
# update_state is required "every single time" to run & produce matrix output (don't know why)
model <- array(NA, dim = c(sir_model$info()$len, n_pars, n_times))

for (t in seq_len(n_times)) {
  model[ , , t] <- sir_model$run(t)
}
# time <- x[1, 1, ] # because in the position of [1, 1, ] is time
# x <- x[-1, , ] # compile all matrix into 1 huge df, delete time (position [-1, , ])
data <- readRDS("inputs/pmcmc_data_week_allAge_ser1_test_2agegroups.rds") %>% 
  glimpse()

sir_data <- dplyr::bind_rows(
  data %>% 
    dplyr::transmute(
      replicate = 1,
      # steps = time_start+1,
      weekly = seq_along(replicate),
      value = count_s1_1,
      compartment = "data_count_s1_1"
    )
  ,
  data %>% 
    dplyr::transmute(
      replicate = 1,
      # steps = time_start+1,
      weekly = seq_along(replicate),
      value = count_s1_2,
      compartment = "data_count_s1_2"
    )
) %>%
  tidyr::complete(weekly, compartment,
                  fill = list(value = 0)) %>% 
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
                                   index == 2 ~ "total N",
                                   index == 3 ~ "total S",
                                   index == 4 ~ "total E",
                                   
                                   index == 5 ~ "total I",
                                   index == 6 ~ "total R",
                                   index == 7 ~ "n_EI1_weekly",
                                   index == 8 ~ "n_EI2_weekly",
                                   index == 9 ~ "beta",
                                   
                                   index == 10 ~ "S <14",
                                   index == 11 ~ "S 15+",
                                   index == 12 ~ "E1",
                                   index == 13 ~ "E2",
                                   index == 14 ~ "model_D1",
                                   index == 15 ~ "model_D2",
                                   index == 16 ~ "R <14",
                                   index == 17 ~ "R 15+",
                                   index == 18 ~ "lambda1",
                                   index == 19 ~ "lambda2"
                                   
                                   
                  )) %>% 
  dplyr::select(-index) %>%
  # step --> time (day) --> week adjustment
  dplyr::mutate(
    # day = (steps - 1) %/% 24 + 1
    # day = (steps - 1) %/% 1 + 1
    day = (steps - 365) %/% 1 + 1
  ) %>%
  dplyr::group_by(replicate, day, compartment) %>%
  summarise(
    value = sum(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(day >= 0) %>% 
  dplyr::mutate(weekly = ceiling(day/7)) %>%
  dplyr::group_by(replicate, weekly, compartment) %>% 
  dplyr::summarise(value = sum(value, na.rm = T),
                   # date = max(date),
                   .groups = "drop") %>% 
  dplyr::ungroup() %>% 
  dplyr::bind_rows(sir_data) %>%
  tidyr::complete(weekly,
                  fill = list(value = 0)
  ) %>%
  dplyr::full_join(
    all_dates
    ,
    by = "weekly"
  ) %>%
  # dplyr::filter(date %in% data$yearWeek) %>%
  glimpse()

p1 <- ggplot(incidence_modelled %>% 
               dplyr::filter(
                 compartment %in% c("model_D1",
                                    "n_EI1_weekly",
                                    "data_count_s1_1"),
                 compartment != "Time",
               )
             ,
             aes(x = yearWeek, y = value,
                 group = interaction(compartment,replicate),
                 colour = compartment)) +
  geom_line() +
  geom_vline(aes(xintercept = as.Date("2010-04-01"),
                 colour = "PCV13 (April 2010)"),
             linetype = "dashed") +
  # scale_y_continuous(limits = c(0, 40)) +
  scale_x_date(limits = c(as.Date(min(all_dates$yearWeek)), as.Date(max(all_dates$yearWeek))),
               date_breaks = "year",
               date_labels = "%Y") +
  ggtitle("Cases (Aggregated by Week) for age 0-14") +
  xlab("Time") +
  ylab("Number of People") +
  theme_bw() +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))

p2 <- ggplot(incidence_modelled %>% 
               dplyr::filter(
                 compartment %in% c("model_D2",
                                    "n_EI2_weekly",
                                    "data_count_s1_2"),
                 compartment != "Time",
               )
             ,
             aes(x = yearWeek, y = value,
                 group = interaction(compartment,replicate),
                 colour = compartment)) +
  geom_line() +
  # scale_y_continuous(limits = c(0, 40)) +
  scale_x_date(limits = c(as.Date(min(all_dates$yearWeek)), as.Date(max(all_dates$yearWeek))),
               date_breaks = "year",
               date_labels = "%Y") +
  ggtitle("Cases (Aggregated by Week) for age 15+") +
  xlab("Time") +
  ylab("Number of People") +
  theme_bw() +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))

p_combined <- cowplot::plot_grid(p1, p2,
                                 nrow =2,
                                 labels = c("A", "B"))


print(p_combined)

# lambda & beta figs ###########################################################
# p1 <- ggplot(incidence_modelled %>%
#                dplyr::filter(
#                  compartment %in% c(
#                    # "lambda1",
#                    "beta"
#                  ),
#                  compartment != "Time",
#                )
#              ,
#              aes(x = yearWeek, y = value,
#                  group = interaction(compartment,replicate),
#                  colour = compartment)) +
#   geom_line() +
#   geom_vline(aes(xintercept = as.Date("2010-04-01"),
#                  colour = "PCV13 (April 2010)"),
#              linetype = "dashed") +
#   scale_x_date(limits = c(as.Date(min(all_dates$yearWeek)), as.Date("2018-03-27")),
#                date_breaks = "year",
#                date_labels = "%Y") +
#   # scale_y_continuous(limits = c(0.1225,0.1250)) +
#   ggtitle("Cases (Aggregated by Week) for age 0-14") +
#   xlab("Time") +
#   ylab("Number of People") +
#   theme_bw() +
#   theme(legend.position = c(0.85, 0.85),
#         legend.title = element_blank(),
#         legend.key.size = unit(0.8, "lines"),
#         legend.text = element_text(size = 10),
#         legend.background = element_rect(fill = "transparent", color = "transparent"))
# 
# p2 <- ggplot(incidence_modelled %>%
#                dplyr::filter(
#                  compartment %in% c(
#                    # "lambda2",
#                    "beta"
#                  ),
#                  compartment != "Time",
#                )
#              ,
#              aes(x = yearWeek, y = value,
#                  group = interaction(compartment,replicate),
#                  colour = compartment)) +
#   geom_line() +
#   scale_x_date(limits = c(as.Date(min(all_dates$yearWeek)), as.Date("2018-03-27")),
#                date_breaks = "year",
#                date_labels = "%Y") +
#   # scale_y_continuous(limits = c(0,0.1)) +
#   ggtitle("Cases (Aggregated by Week) for age 15+") +
#   xlab("Time") +
#   ylab("Number of People") +
#   theme_bw() +
#   theme(legend.position = c(0.85, 0.85),
#         legend.title = element_blank(),
#         legend.key.size = unit(0.8, "lines"),
#         legend.text = element_text(size = 10),
#         legend.background = element_rect(fill = "transparent", color = "transparent"))
# 
# p_combined <- cowplot::plot_grid(p1, p2,
#                                  nrow =2,
#                                  labels = c("A", "B"))
# 
# 
# print(p_combined)



