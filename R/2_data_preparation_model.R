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

# combine data into 1 csv
dat_week_12F_3ageG <- below_2 %>% 
  dplyr::rename(cases_1 = cases) %>% 
  dplyr::left_join(adults %>% 
                     dplyr::rename(cases_2 = cases)
                   , by = "time"
                   ) %>% 
  dplyr::left_join(elderly %>% 
                     dplyr::rename(cases_3 = cases)
                   , by = "time"
  ) %>% 
  dplyr::rename(week = time) # naming error in mcState -_-)

write.csv(dat_week_12F_3ageG, "inputs/incidence_week_12F_3ageG_all.csv", row.names = FALSE)

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
# combine data into 1 csv
dat_week_12F_6ageG <- read.csv("inputs/incidence_week_12F_6ageG_1.csv") %>% 
  dplyr::rename(cases_1 = cases) %>% 
  dplyr::left_join(
    read.csv("inputs/incidence_week_12F_6ageG_2.csv") %>% 
      dplyr::rename(cases_2 = cases)
    , by = "time"
    ) %>% 
  dplyr::left_join(
    read.csv("inputs/incidence_week_12F_6ageG_3.csv") %>% 
      dplyr::rename(cases_3 = cases)
    , by = "time"
  ) %>% 
  dplyr::left_join(
    read.csv("inputs/incidence_week_12F_6ageG_4.csv") %>% 
      dplyr::rename(cases_4 = cases)
    , by = "time"
  ) %>% 
  dplyr::left_join(
    read.csv("inputs/incidence_week_12F_6ageG_5.csv") %>% 
      dplyr::rename(cases_5 = cases)
    , by = "time"
  ) %>% 
  dplyr::left_join(
    read.csv("inputs/incidence_week_12F_6ageG_6.csv") %>% 
      dplyr::rename(cases_6 = cases)
    , by = "time"
  ) %>% 
  dplyr::rename(week = time) # naming error in mcState -_-)

write.csv(dat_week_12F_6ageG, "inputs/incidence_week_12F_6ageG_all.csv", row.names = FALSE)


# mcstate data preparation #####################################################
# load epidata
dat_c <- read.csv("raw_data/12F_Jan_2025_combined_cleaned.csv") %>% 
  dplyr::filter(ageGroup3 != "Unknown")

# load genomic data
gen <- read.csv("raw_data/genomic_data_cleaned.csv") %>% 
  dplyr::filter(!is.na(strain),
                collection_date >= as.Date("2017-08-01")) %>%  # after 2017-08-01
  glimpse()

earlier_ne_df <- read.csv("raw_data/GPSC55_mlesky_cleaned_interpolated_predictedModel_binom.csv") %>% 
  glimpse()

# load ne
# ne_55 <- read.csv("raw_data/GPSC55_mlesky_cleaned.csv") %>% 
#   dplyr::filter(date >= min(dat_c$week_date))
interpolated_ne <- read.csv("raw_data/GPSC55_mlesky_cleaned_interpolated.csv") %>% 
  dplyr::select(change_Ne, yearWeek) %>% 
  dplyr::rename(Ne = change_Ne) %>% 
  dplyr::mutate(yearWeek = as.Date(yearWeek))

# non-heterogeneity (allAges), weekly
allAges_weekly <- dat_c %>% 
  dplyr::mutate(week_date = as.Date(week_date),
                iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                yearWeek =ISOweek::ISOweek2date(iso_week)
  ) %>% 
  dplyr::group_by(yearWeek) %>% 
  dplyr::summarise(count_serotype = sum(counts)) %>% 
  dplyr::ungroup() %>% 
  dplyr::full_join(
    dplyr::bind_rows(
      gen %>% 
        dplyr::filter(strain == "GPSC55") %>% 
        dplyr::mutate(week_date = as.Date(week_date),
                      iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                      yearWeek =ISOweek::ISOweek2date(iso_week)
        ) %>% 
        dplyr::group_by(yearWeek) %>% 
        dplyr::summarise(count_WGS_GPSC55 = n()) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(yearWeek = as.Date(yearWeek))
      ,
      earlier_ne_df %>% 
        dplyr::select(yearWeek, predicted_count_GPSC55) %>% 
        dplyr::mutate(yearWeek = as.Date(yearWeek)) %>% 
        dplyr::filter(yearWeek <= as.Date("2017-08-01")) %>% 
        dplyr::rename(count_WGS_GPSC55 = predicted_count_GPSC55)
    )
  ,
  by = "yearWeek"
) %>% 
  dplyr::full_join(
    gen %>% 
      dplyr::filter(strain == "non55") %>% 
      dplyr::mutate(week_date = as.Date(week_date),
                    iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                    yearWeek =ISOweek::ISOweek2date(iso_week)
      ) %>% 
      dplyr::group_by(yearWeek) %>% 
      dplyr::summarise(count_WGS_non55 = n()) %>% 
      dplyr::ungroup()
    ,
    by = "yearWeek"
  ) %>% 
  dplyr::full_join(
    interpolated_ne
    ,
    by = "yearWeek"
    ,
    relationship = "many-to-many"
  ) %>%
  dplyr::mutate(
    count_serotype = as.numeric(count_serotype),
    count_WGS_GPSC55 = as.numeric(count_WGS_GPSC55),
    count_WGS_non55 = as.numeric(count_WGS_non55),
    Ne = as.numeric(Ne)
  ) %>% 
  # tidyr::pivot_longer(
  #   cols = starts_with(c("count_")), # ignore Ne at the moment
  #   names_to = "type",
  #   values_to = "count"
  # ) %>% 
  dplyr::arrange(yearWeek) %>% 
  dplyr::mutate(yearWeek = as.Date(yearWeek),
                day = as.numeric(round((yearWeek - as.Date("1987-09-14"))))) %>%
  mcstate::particle_filter_data(., time = "day", rate = 1, initial_time = 0) %>%
  glimpse()

saveRDS(allAges_weekly, "inputs/pmcmc_data_week_allAge.rds")


# test viz combined GPSC55
ggplot(allAges_weekly
       , aes(x = yearWeek)) +
  geom_line(aes(y = count_serotype), colour = "maroon") +
  geom_line(aes(y = count_WGS_GPSC55), colour = "black") +
  geom_line(aes(y = count_WGS_non55), colour = "darkgreen", size = 1.5) +
  geom_line(aes(y = Ne), colour = "gold2") +
  geom_vline(xintercept = as.Date("2017-08-01"), color = "steelblue", linetype = "dashed") +
  scale_x_date(limits = c(as.Date("2001-01-01"), as.Date("2022-06-01")), 
               date_breaks = "1 year",
               date_labels = "%Y") +
  theme_bw() +
  labs(
    title = "GPSC55 Counts Prediction + Real Data",
    y = "GPSC55 counts"
  ) +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", colour = "transparent"))


# ageGroup3, weekly
ageGroup3_weekly <- dat_c %>% 
  dplyr::mutate(week_date = as.Date(week_date),
                iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                yearWeek =ISOweek::ISOweek2date(iso_week)
  ) %>% 
  dplyr::group_by(yearWeek, ageGroup3) %>% 
  dplyr::summarise(count_serotype = sum(counts), .groups = "drop") %>% 
  dplyr::ungroup() %>% 
  tidyr::pivot_wider(
    names_from = ageGroup3,
    values_from = count_serotype
  ) %>% 
  dplyr::rename(
    count_serotype_1 = `<2`,
    count_serotype_2 = `2-64`,
    count_serotype_3 = `65+`
  ) %>% 
  dplyr::full_join(
    gen %>% 
      dplyr::filter(strain == "GPSC55") %>% 
      dplyr::mutate(week_date = as.Date(week_date),
                    iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                    yearWeek =ISOweek::ISOweek2date(iso_week)
      ) %>% 
      dplyr::group_by(yearWeek, ageGroup3) %>% 
      dplyr::summarise(count_WGS_GPSC55 = n(), .groups = "drop") %>% 
      dplyr::ungroup() %>% 
      tidyr::pivot_wider(
        names_from = ageGroup3,
        values_from = count_WGS_GPSC55
      ) %>% 
      dplyr::rename(
        count_WGS_GPSC55_1 = `<2`,
        count_WGS_GPSC55_2 = `2-64`,
        count_WGS_GPSC55_3 = `65+`
      )
    ,
    by = c("yearWeek")
    ) %>% 
  dplyr::full_join(
    gen %>% 
      dplyr::filter(strain == "non55") %>% 
      dplyr::mutate(week_date = as.Date(week_date),
                    iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                    yearWeek =ISOweek::ISOweek2date(iso_week)
      ) %>% 
      dplyr::group_by(yearWeek, ageGroup3) %>% 
      dplyr::summarise(count_WGS_non55 = n(), .groups = "drop") %>% 
      dplyr::ungroup() %>% 
      tidyr::pivot_wider(
        names_from = ageGroup3,
        values_from = count_WGS_non55
      ) %>% 
      dplyr::rename(
        count_WGS_non55_1 = `<2`,
        count_WGS_non55_2 = `2-64`,
        count_WGS_non55_3 = `65+`
      )
    ,
    by = c("yearWeek")
  ) %>% 
  dplyr::full_join(
    interpolated_ne
    ,
    by = "yearWeek"
  ) %>%
  # tidyr::pivot_longer(
  #   cols = starts_with("count_"),
  #   names_to = "type",
  #   values_to = "count"
  # ) %>% 
  dplyr::arrange(yearWeek) %>% 
  dplyr::mutate(yearWeek = as.Date(yearWeek),
                day = as.numeric(round((yearWeek - as.Date("1987-09-14"))))) %>%
  mcstate::particle_filter_data(., time = "day", rate = 1, initial_time = 0) %>%
  glimpse()

saveRDS(ageGroup3_weekly, "inputs/pmcmc_data_week_ageGroup3.rds")










