# Total population data by age, year for each region
# SOURCE: https://www.nomisweb.co.uk/

library(tidyverse)

# compartment separation for 12F
age_proportion <- dplyr::left_join(
  read.csv("raw_data/nomis_population_long.csv") %>% 
    # region stratification is not needed
    dplyr::group_by(year) %>% 
    dplyr::summarise(PopSize_year = sum(PopSize)) %>% 
    dplyr::ungroup()
  ,
  read.csv("raw_data/nomis_population_long.csv") %>% 
  # region stratification is not needed
  dplyr::mutate(
    Age = as.numeric(Age),
    ageGroup12F = case_when(
      Age <= 44 ~ "0-44",
      Age > 44 ~ "45+"
    )
  ) %>% 
  dplyr::group_by(ageGroup12F, year) %>% 
  dplyr::summarise(PopSize12F = sum(PopSize)) %>% 
  dplyr::ungroup()
  ,
  by = "year"
  ) %>% 
  dplyr::mutate(
    ageGroup12F = factor(ageGroup12F,
                       levels = c("0-44", "45+")),
    PopProp = round(PopSize12F/PopSize_year, 2)
  ) %>% 
  # view() %>% 
  glimpse()

# younger people (0-44) was consistently 60% of the total population


# load model result
# (based on sir_stochastic_allAge_2post_pmcmc_picts.R, saved in raw_data)
incidence_modelled <- dplyr::bind_rows(
  read.csv("raw_data/incidence_modelled_GPSC55.csv") %>% 
    dplyr::filter(compartment == "model_D1") %>% 
    dplyr::mutate(yearWeek = as.Date(yearWeek),
                  year = year(yearWeek)) %>% 
    dplyr::left_join(
      # age proportion in England
      dplyr::left_join(
        read.csv("raw_data/nomis_population_long.csv") %>% 
          # region stratification is not needed
          dplyr::group_by(year) %>% 
          dplyr::summarise(PopSize_year = sum(PopSize)) %>% 
          dplyr::ungroup()
        ,
        read.csv("raw_data/nomis_population_long.csv") %>% 
          # region stratification is not needed
          dplyr::group_by(ageGroup6, year) %>% 
          dplyr::summarise(PopSize6 = sum(PopSize)) %>% 
          dplyr::ungroup()
        ,
        by = "year"
      ) %>% 
        dplyr::mutate(
          ageGroup6 = factor(ageGroup6,
                             levels = c("<2", "2-4", "5-14", "15-44", "45-64", "65+")),
          PopProp = round(PopSize6/PopSize_year, 1)
        ) %>% 
        glimpse()
      ,
      by = "year",
      relationship = "many-to-many"
    ) %>% 
    dplyr::filter(ageGroup6 %in% c("<2", "2-4", "5-14", "15-44")) %>% # D1 (0-44)
    dplyr::mutate(case_modelled = round(value*PopProp, 1)) %>% 
    dplyr::filter(!is.na(ageGroup6)) %>% 
    dplyr::arrange(yearWeek)
  ,
  read.csv("raw_data/incidence_modelled_GPSC55.csv") %>% 
    dplyr::filter(compartment == "model_D2") %>% 
    dplyr::mutate(yearWeek = as.Date(yearWeek),
                  year = year(yearWeek)) %>% 
    dplyr::left_join(
      # age proportion in England
      dplyr::left_join(
        read.csv("raw_data/nomis_population_long.csv") %>% 
          # region stratification is not needed
          dplyr::group_by(year) %>% 
          dplyr::summarise(PopSize_year = sum(PopSize)) %>% 
          dplyr::ungroup()
        ,
        read.csv("raw_data/nomis_population_long.csv") %>% 
          # region stratification is not needed
          dplyr::group_by(ageGroup6, year) %>% 
          dplyr::summarise(PopSize6 = sum(PopSize)) %>% 
          dplyr::ungroup()
        ,
        by = "year"
      ) %>% 
        dplyr::mutate(
          ageGroup6 = factor(ageGroup6,
                             levels = c("<2", "2-4", "5-14", "15-44", "45-64", "65+")),
          PopProp = round(PopSize6/PopSize_year, 1)
        ) %>% 
        glimpse()
      ,
      by = "year",
      relationship = "many-to-many"
    ) %>% 
    dplyr::filter(ageGroup6 %in% c("45-64", "65+")) %>% # D2 (45-64 & 65+) was 50-50 in proportion
    dplyr::mutate(case_modelled = round(value/2, 1)) %>%  # no need to calculate proportion
    dplyr::filter(!is.na(ageGroup6)) %>% 
    dplyr::arrange(yearWeek)
)

ggplot(incidence_modelled,
       aes(x = yearWeek, y = case_modelled,
           group = ageGroup6,
           colour = ageGroup6)) +
  geom_line() +
  theme_bw()


# comparison with the real case counts data
serotype12F_data <- read.csv("raw_data/12F_Jan_2025_combined_cleaned.csv") %>% 
  glimpse()

GPSC55_data <- read.csv("raw_data/genomic_data_cleaned.csv") %>% 
  dplyr::mutate(week_date = as.Date(week_date),
                iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                yearWeek =ISOweek::ISOweek2date(iso_week)
  ) %>% 
  glimpse()


all_week <- data.frame(week_date = seq.Date(from = min(as.Date(serotype12F_data$week_date)),
                                            to = max(as.Date(serotype12F_data$week_date)), 
                                            by = 1)) %>% 
  dplyr::mutate(week_step = 1:nrow(.),
                iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                yearWeek =ISOweek::ISOweek2date(iso_week)
  ) %>% 
  distinct(yearWeek) %>% 
  glimpse()

# serotype 12F
serotype12F_data_ageGroup6 <- dplyr::left_join(
  all_week
  ,
  serotype12F_data %>% 
    dplyr::mutate(week_date = as.Date(week_date),
                  iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                  yearWeek =ISOweek::ISOweek2date(iso_week)
    ) %>% 
    dplyr::group_by(yearWeek, ageGroup6) %>% 
    dplyr::summarise(case_data_12F = sum(counts)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
      case_data_12F = as.numeric(case_data_12F),
      ageGroup6 = case_when(
        is.na(ageGroup6) ~ "Unknown",
        TRUE ~ ageGroup6
      )
    ) %>% 
    dplyr::arrange(yearWeek)
  ,
  by = "yearWeek"
) %>% 
  dplyr::filter(ageGroup6 != "Unknown") %>% 
  glimpse()

ggplot(serotype12F_data_ageGroup6,
       aes(x = yearWeek, y = case_data_12F,
           group = ageGroup6,
           colour = ageGroup6)) +
  geom_line() +
  theme_bw()

# GPSC55
GPSC55_data_ageGroup6 <- dplyr::left_join(
  all_week
  ,
  GPSC55_data %>% 
    dplyr::filter(strain == "GPSC55") %>% 
    dplyr::mutate(week_date = as.Date(week_date),
                  iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                  yearWeek =ISOweek::ISOweek2date(iso_week)
    ) %>% 
    dplyr::group_by(yearWeek, ageGroup6) %>% 
    dplyr::summarise(case_data_GPSC55 = n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
      case_data_GPSC55 = as.numeric(case_data_GPSC55),
      ageGroup6 = case_when(
        is.na(ageGroup6) ~ "Unknown",
        TRUE ~ ageGroup6
      )
    ) %>% 
    dplyr::arrange(yearWeek)
  ,
  by = "yearWeek"
) %>% 
  dplyr::filter(ageGroup6 != "Unknown") %>% 
  glimpse()

ggplot(GPSC55_data_ageGroup6,
       aes(x = yearWeek, y = case_data_GPSC55,
           group = ageGroup6,
           colour = ageGroup6)) +
  geom_line() +
  theme_bw()

# combine df
png("outputs/genomics/trials_GPSC55/trial_100250/figs/model_vs_data.png",
    width = 24, height = 17, unit = "cm", res = 600)
dplyr::bind_rows(
  incidence_modelled %>% 
    dplyr::transmute(
      yearWeek = as.Date(yearWeek),
      ageGroup6 = ageGroup6,
      version = "model",
      counts = case_modelled
    )
  ,
  serotype12F_data_ageGroup6 %>% 
    dplyr::transmute(
      yearWeek = as.Date(yearWeek),
      ageGroup6 = ageGroup6,
      version = "data serotype12F",
      counts = case_data_12F
    )
  ,
  GPSC55_data_ageGroup6 %>% 
    dplyr::transmute(
      yearWeek = as.Date(yearWeek),
      ageGroup6 = ageGroup6,
      version = "data GPSC55",
      counts = case_data_GPSC55
    )
) %>% 
  dplyr::mutate(
    ageGroup6 = factor(ageGroup6,
                       levels = c("<2", "2-4", "5-14",
                                  "15-44", "45-64", "65+")),
    version = factor(version,
                     levels = c("data serotype12F", "data GPSC55", "model"))
  ) %>% 
  ggplot(aes(x = yearWeek, y = counts,
             color = version, group = version)) +
  geom_line() +
  geom_vline(xintercept = as.Date("2017-09-01"), color = "steelblue", linetype = "dashed") +
  labs(
    title = "Serotype 12F & GPSC55: model vs data"
  ) +
  facet_wrap(~ ageGroup6, 
             # scales = "free_y"
  ) +
  theme_bw()
dev.off()
