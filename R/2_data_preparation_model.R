rm(list = ls())

library(tidyverse)

if (!dir.exists("inputs")) {
  dir.create("inputs")
}

# 12F all ages
dat_combined_week <- read.csv("raw_data/12F_Jan_2025_combined_cleaned.csv") %>% 
  dplyr::group_by(week_date) %>% 
  dplyr::summarise(counts = sum(counts)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(week_date = as.Date(week_date))

all_date <- data.frame(week_date = seq.Date(from = min(dat_combined_week$week_date),
                                          to = max(dat_combined_week$week_date), 
                                          by = 1)) %>% 
  dplyr::mutate(week_step = 1:nrow(.))

incidence_week_allAge <- dplyr::left_join(all_date, 
                                          dat_combined_week, 
                                          by = "week_date"
                                          ) %>% 
  dplyr::select(week_step, counts) %>% 
  dplyr::rename(time = week_step,
                cases = counts) # Annoying name requirement inputs to monty

write.csv(incidence_week_allAge, "inputs/incidence_week_12F_allAge.csv", row.names = FALSE)


# 12F 3 age groups
dat_combined_week <- read.csv("raw_data/12F_Jan_2025_combined_cleaned.csv") %>% 
  dplyr::group_by(week_date, ageGroup3) %>% 
  dplyr::summarise(counts = sum(counts)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(week_date = as.Date(week_date))

all_date <- data.frame(week_date = seq.Date(from = min(dat_combined_week$week_date),
                                            to = max(dat_combined_week$week_date), 
                                            by = 1)) %>% 
  dplyr::mutate(week_step = 1:nrow(.))

# Annoying R dimension requirement
dat_combined_week_age_ <- list()
incidence_week_3ageG_ <- list()

for (i in 1:3) {
  ageGroups <- c("<2", "2-64", "65+") #, "Unknown")
  
  dat_combined_week_age_[[i]] <- dat_combined_week %>% 
    dplyr::filter(ageGroup3 == ageGroups[i])
  
  incidence_week_3ageG_[[i]] <- dplyr::left_join(all_date, 
                                                 dat_combined_week_age_[[i]], 
                                                 by = "week_date"
  ) %>% 
    dplyr::select(week_step, counts) %>% 
    dplyr::rename(time = week_step,
                  cases = counts) # Annoying name requirement inputs to monty
  
  write.csv(incidence_week_3ageG_[[i]],
            paste0("inputs/incidence_week_12F_3ageG_", i, ".csv"),
            row.names = FALSE)
}

# Test load agegroups
below_2 <- read.csv("inputs/incidence_week_12F_3ageG_1.csv")
adults <- read.csv("inputs/incidence_week_12F_3ageG_2.csv")
elderly <- read.csv("inputs/incidence_week_12F_3ageG_3.csv")

plot(below_2$time, below_2$cases, type = "p", col = "red", ylim = c(0,30))
points(adults$time, adults$cases, col = "green")
points(elderly$time, elderly$cases, col = "blue")


# 12F 6 age groups
dat_combined_week <- read.csv("raw_data/12F_Jan_2025_combined_cleaned.csv") %>% 
  dplyr::group_by(week_date, ageGroup6) %>% 
  dplyr::summarise(counts = sum(counts)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(week_date = as.Date(week_date))

all_date <- data.frame(week_date = seq.Date(from = min(dat_combined_week$week_date),
                                            to = max(dat_combined_week$week_date), 
                                            by = 1)) %>% 
  dplyr::mutate(week_step = 1:nrow(.))

# Annoying R dimension requirement
dat_combined_week_age_ <- list()
incidence_week_6ageG_ <- list()

for (i in 1:6) {
  ageGroups <- c("<2", "2-4", "5-14", "15-44", "45-64", "65+") #, "Unknown")
  
  dat_combined_week_age_[[i]] <- dat_combined_week %>% 
    dplyr::filter(ageGroup6 == ageGroups[i])
  
  incidence_week_6ageG_[[i]] <- dplyr::left_join(all_date, 
                                                 dat_combined_week_age_[[i]], 
                                                 by = "week_date"
  ) %>% 
    dplyr::select(week_step, counts) %>% 
    dplyr::rename(time = week_step,
                  cases = counts) # Annoying name requirement inputs to monty
  
  write.csv(incidence_week_6ageG_[[i]],
            paste0("inputs/incidence_week_12F_6ageG_", i, ".csv"),
            row.names = FALSE)
}

# Test load agegroups




