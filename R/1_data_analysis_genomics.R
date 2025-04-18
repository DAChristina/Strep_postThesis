library(tidyverse)

# epidata is based on combined new data weekly;
# either non-heterogeneity, 3 or 6 ageG

# for gendata I wanna see what was analysed previously.
# I use floor_date instead of ceiling_date
Ne <- read.csv("raw_data/gen_lw/Ne_over_time.csv") %>% 
  dplyr::mutate(year = round(Time)) %>%
  # dplyr::filter(between(year, 2011, 2022)) %>% 
  # dplyr::mutate(step = round((year - 2010.5) * 365)) %>%
  # dplyr::select(year, step, Ne, Ne_lb, Ne_ub) %>%
  glimpse()

Ne_skygrid_adaptive <- read.csv("raw_data/gen_lw/GPSC55_mlesky.csv") %>% 
  dplyr::mutate(year = floor(Year),
                Model = cumsum(c(0, diff(year) < 0)) + 1) %>% 
  glimpse()

Ne_skygrid_adaptive %>% 
  dplyr::group_by(year) %>% 
  count() %>% 
  glimpse()

Ne_skygrid_adaptive %>% 
  ggplot(aes(x = Year)) +
  geom_ribbon(aes(ymin = Lower_range, ymax = Upper_range), alpha = 0.2) +
  geom_line(aes(y = Ne)) +
  facet_wrap(vars(Model)) +
  coord_cartesian(xlim = c(2010, 2022))

Ne <- Ne_skygrid_adaptive %>% 
  dplyr::filter(between(year, 2011, 2022),
                Model == 3) %>% 
  dplyr::transmute(year = floor(Year),
                   step = round((Year - 2010.5) * 365),
                   Ne = Ne, 
                   Ne_lb = Lower_range,
                   Ne_ub = Upper_range) %>% 
  glimpse()
## Model 3 looks like skygrowth

# tbh I have no idea what is this "invasiveness" file is about.
invasiveness <- read.csv("raw_data/gen_lw/invasiveness_12F.csv") %>%
  dplyr::select(X, Adjusted_invasiveness, Age_group) %>% 
  dplyr::mutate(Age_group = tolower(Age_group)) %>% 
  tidyr::pivot_wider(names_from = Age_group, values_from = Adjusted_invasiveness) %>%
  glimpse()


gen <-  read.csv("raw_data/gen_lw/genomic_epi_12F.csv") %>% 
  dplyr::mutate(strain = if_else(GPSC == 55, "GPSC55", "non55"),
                ageGroup3 = case_when(
                  Age < 2 ~ "<2",
                  Age >= 2 & Age < 65 ~ "2-64",
                  Age >= 65 ~ "65+",
                  is.na(Age) ~ "Unknown"
                ),
                ageGroup6 = case_when(
                  Age < 2 ~ "<2",
                  Age >= 2 & Age < 5 ~ "2-4",
                  Age >= 5 & Age < 15 ~ "5-14",
                  Age >= 15 & Age < 45 ~ "15-44",
                  Age >= 45 & Age < 65 ~ "45-64",
                  Age >= 65 ~ "65+",
                  is.na(Age) ~ "Unknown"
                ),
                collection_date = lubridate::dmy(collection_date),
                month = lubridate::month(collection_date),
                # year = if_else(month < 7, YEAR + 0.5, YEAR + 1.5), # this is epi year
                year = lubridate::year(collection_date),
                month = month(collection_date),
                yearMonth = as.Date(paste0(format(collection_date, "%Y-%m"), "-01")), # year-month as date
                week_date = lubridate::floor_date(collection_date, "week")
                ) %>% 
  dplyr::select(-ngsid, -Age, -YEAR) %>% 
  # dplyr::filter(!is.na(strain), year > 2018) %>% 
  # view() %>% 
  glimpse()

write.csv(gen, "raw_data/genomic_data_cleaned.csv", row.names = FALSE)

# need to be specified by ageGroups
# yearly_cases_by_strain <- gen %>% 
#   # dplyr::filter(!is.na(GPSC)) %>% 
#   dplyr::group_by(year, ageGroup3, strain) %>% 
#   dplyr::count() %>% 
#   glimpse()
# 
# week_dately_cases_by_strain <- gen %>% 
#   dplyr::group_by(week_date, ageGroup3, strain) %>% 
#   dplyr::count() %>% 
#   glimpse()
# 
# monthly_cases_by_strain <- gen %>% 
#   dplyr::group_by(yearMonth, ageGroup3, strain) %>% 
#   dplyr::count() %>% 
#   glimpse()
# 
# yearly_cases_by_strain %>% 
#   ggplot(aes(x = year, y = n)) +
#   geom_line() +
#   facet_grid(vars(ageGroup3), vars(strain), scales = "free_y") +
#   geom_vline(aes(xintercept = 2019.5), lty = 2)
# 
# week_dately_cases_by_strain %>% 
#   ggplot(aes(x = week_date, y = n)) +
#   geom_line() +
#   facet_grid(vars(ageGroup3), vars(strain), scales = "free_y") +
#   geom_vline(aes(xintercept = 2019.5), lty = 2)
# 
# monthly_cases_by_strain %>% 
#   tidyr::pivot_wider(names_from = strain, values_from = n) %>% 
#   dplyr::mutate(non55 = replace_na(non55, 0),
#                 GPSC55 = replace_na(GPSC55, 0),
#                 total = non55 + GPSC55,
#                 p_non55 = non55 / total) %>% 
#   ggplot(aes(x = yearMonth, y = p_non55, colour = ageGroup3)) +
#   geom_line()
# 
# monthly_cases_by_strain %>% 
#   ggplot(aes(x = yearMonth, y = n)) +
#   geom_line() +
#   facet_grid(vars(ageGroup3), vars(strain), scales = "free_y") +
#   geom_vline(aes(xintercept = 2019.5), lty = 2)
# 
# cases_55C <- yearly_cases_by_strain %>% 
#   tidyr::pivot_wider(names_from = strain, values_from = n, values_fill = 0) %>% 
#   dplyr::mutate(`12F` = GPSC55 + non55) %>% 
#   tidyr::pivot_wider(names_from = ageGroup3, values_from = GPSC55:`12F`,
#                      names_glue = "cases_{ageGroup3}_{.value}") %>% 
#   glimpse()
## TODO we have collection dates already so can fit to < yearly data - go for weekly?






