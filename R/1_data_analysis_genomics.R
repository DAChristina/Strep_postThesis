library(tidyverse)

# epidata is based on combined new data weekly;
# either non-heterogeneity, 3 or 6 ageG

# for gendata I wanna see what was analysed previously.
# I use floor_date instead of ceiling_date
Ne <- read.csv("raw_data/gen_lw/Ne_over_time.csv") %>% 
  dplyr::mutate(date = as.Date((Time - 1970) * 365.25, origin = "1970-01-01"),
                year = lubridate::year(date),
                month = lubridate::month(date),
                yearMonth = as.Date(paste0(format(date, "%Y-%m"), "-01")) # year-month as date
  ) %>% 
  # dplyr::filter(between(year, 2011, 2022)) %>% 
  # dplyr::mutate(step = round((year - 2010.5) * 365)) %>%
  # dplyr::select(year, step, Ne, Ne_lb, Ne_ub) %>%
  glimpse()

Ne_skygrid_adaptive <- read.csv("raw_data/gen_lw/GPSC55_mlesky.csv") %>% 
  dplyr::mutate(date = as.Date((Year - 1970) * 365.25, origin = "1970-01-01"),
                year = lubridate::year(date),
                month = lubridate::month(date),
                yearMonth = as.Date(paste0(format(date, "%Y-%m"), "-01")),
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
# Lilith use the Ne_skygrid_adaptive (GPSC55_mlesky.csv), model 3

# I re-open and prepare Ne data per-year, month and week
ne_55 <- read.csv("raw_data/gen_lw/GPSC55_mlesky.csv") %>% 
  dplyr::mutate(date = as.Date((Year - 1970) * 365.25, origin = "1970-01-01"),
                year = lubridate::year(date),
                month = lubridate::month(date),
                yearMonth = as.Date(paste0(format(date, "%Y-%m"), "-01")),
                Model = cumsum(c(0, diff(year) < 0)) + 1) %>% 
  dplyr::filter(Model == 3) %>%
  dplyr::rename(Ne_lb = Lower_range,
                Ne_ub = Upper_range) %>%
  dplyr::select(-c("Strain", "Model")) %>% 
  glimpse()

write.csv(ne_55, "raw_data/GPSC55_mlesky_cleaned.csv", row.names = FALSE)

# use Ne data to predict previous GPSC prevalence ##############################
new_dates <- seq(from = as.Date(min(ne_55$date)), to = as.Date(max(ne_55$date)), by = "1 week") %>% 
  glimpse()

fn_Ne    <- splinefun(x = as.Date(ne_55$date), y = ne_55$Ne, method = "natural")
fn_Ne_lo <- splinefun(x = as.Date(ne_55$date), y = ne_55$Ne_lb, method = "natural")
fn_Ne_up <- splinefun(x = as.Date(ne_55$date), y = ne_55$Ne_ub, method = "natural")

spline_fit_Ne    <- fn_Ne(new_dates)
spline_fit_Ne_lo <- fn_Ne_lo(new_dates)
spline_fit_Ne_up <- fn_Ne_up(new_dates)

time_diffs_weeks <- as.numeric(diff(new_dates)) / 7
change_Ne <- abs(diff(spline_fit_Ne)) / time_diffs_weeks
change_Ne_lo <- abs(diff(spline_fit_Ne_lo)) / time_diffs_weeks
change_Ne_up <- abs(diff(spline_fit_Ne_up)) / time_diffs_weeks

interpolated_df <- tibble(
  date = new_dates[-1],
  itr_Ne = spline_fit_Ne[-1],
  itr_Ne_lo = spline_fit_Ne_lo[-1],
  itr_Ne_up = spline_fit_Ne_up[-1],
  time_diffs_weeks = time_diffs_weeks,
  change_Ne = change_Ne,
  change_Ne_lo = change_Ne_lo,
  change_Ne_up = change_Ne_up
) %>%
  mutate(
    iso_week = paste0(year(date), "-W", sprintf("%02d", week(date)), "-1"),
    yearWeek = ISOweek::ISOweek2date(iso_week)
  ) %>%
  glimpse()

write.csv(interpolated_df, "raw_data/GPSC55_mlesky_cleaned_interpolated.csv", row.names = FALSE)


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
                month = lubridate::month(collection_date),
                yearMonth = as.Date(paste0(format(collection_date, "%Y-%m"), "-01")), # year-month as date
                week_date = lubridate::floor_date(collection_date, "week")
                ) %>% 
  dplyr::select(-ngsid, -Age, -YEAR) %>% 
  # dplyr::filter(!is.na(strain), year > 2018) %>% 
  # view() %>% 
  glimpse()

write.csv(gen, "raw_data/genomic_data_cleaned.csv", row.names = FALSE)


# load gen & interpolated_df first
test <- gen %>% 
  dplyr::filter(strain == "GPSC55") %>% 
  dplyr::mutate(week_date = as.Date(week_date),
                iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                yearWeek =ISOweek::ISOweek2date(iso_week)
  ) %>% 
  dplyr::group_by(yearWeek) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(
    interpolated_df %>%
      dplyr::select(yearWeek, contains("Ne"))
    ,
    by = c("yearWeek")
  ) %>% 
  glimpse()

################################################################################
selected_GPSC55 <- test %>% 
  dplyr::filter(yearWeek >= as.Date("2017-08-01")) %>% 
  dplyr::mutate(sin_week = sin(2*pi*lubridate::isoweek(yearWeek)/52),
                cos_week = cos(2*pi*lubridate::isoweek(yearWeek)/52)) %>% 
  glimpse()

# test other models
model_gam <- mgcv::gam(count ~ s(change_Ne), data = selected_GPSC55)
model_gam_bayesian <- brms::brm(count ~ s(change_Ne), data = selected_GPSC55,
                                family = poisson(), chains = 4, cores = 2,
                                iter = 1000)
model_glm <- stats::glm(count ~ sin_week + cos_week + change_Ne, data = selected_GPSC55,
                        family = poisson)

AIC(model_gam, model_glm) # only works for frequentist models

# save model (use load("model.RData") to load model later)
saveRDS(model_gam, file = "raw_data/model_Ne_gam.rds")
saveRDS(model_gam_bayesian, file = "raw_data/model_Ne_gam_bayes.rds")
saveRDS(model_glm, file = "raw_data/model_Ne_glm.rds")

earlier_ne_df <- interpolated_df %>%
  dplyr::filter(yearWeek < as.Date("2017-08-01")) %>% # midpoint; they started intensively sequenced GPSC55 ("2017-07-01")
  dplyr::mutate(sin_week = sin(2*pi*lubridate::isoweek(yearWeek)/52),
                cos_week = cos(2*pi*lubridate::isoweek(yearWeek)/52)) %>%
  glimpse()

earlier_ne_df <- earlier_ne_df %>% 
  dplyr::mutate(predicted_GPSC55_gam = predict(model_gam, newdata = earlier_ne_df),
                predicted_GPSC55_gam_bayesian = predict(model_gam_bayesian, newdata = earlier_ne_df),
                predicted_GPSC55_gam_bayesian_mean = rowMeans(predicted_GPSC55_gam_bayesian),
                predicted_GPSC55_glm = predict(model_glm, newdata = earlier_ne_df)) %>% 
  glimpse()

# save predicted model result for previous GPSC55 cases pre-August 2017
write.csv(earlier_ne_df, "raw_data/GPSC55_mlesky_cleaned_interpolated_predictedModel.csv", row.names = FALSE)

# test GPSC55 previous WGS
combined <- dplyr::bind_rows(
  selected_GPSC55 %>% 
    select(yearWeek, count) %>%
    mutate(source = "1. WGS data")
  ,
  earlier_ne_df %>%
    select(yearWeek, contains("predicted_GPSC55"),
           -predicted_GPSC55_gam_bayesian) %>% # weird matrix format (4 chains)
    tidyr::pivot_longer(
      cols = contains("predicted_GPSC55"),
      names_to = "source",
      values_to = "count"
    ) %>%
    dplyr::mutate(
      source = case_when(
        source == "predicted_GPSC55_gam" ~ "3.1. Predicted (GAM)",
        source == "predicted_GPSC55_gam_bayesian_mean" ~ "3.2. Predicted mean (Bayesian GAM)",
        source == "predicted_GPSC55_glm" ~ "3.3. Predicted (GLM)",
      )
    ) %>%
    dplyr::filter(source != "predicted_GPSC55_gam_bayesian") %>% 
    # weird array conversion
    # dplyr::mutate(
    #   count = as.data.frame(count)
    # ) %>%
    unnest(cols = count)
  ,
  interpolated_df %>%
    dplyr::select(yearWeek, change_Ne) %>% 
    rename(count = change_Ne) %>% 
    mutate(source = "2. Interpolated Ne")
) %>% 
  # weird array conversion
  dplyr::mutate(
    count = as.data.frame(count)
  ) %>% 
  unnest(cols = count) %>% 
  glimpse()


ggplot(combined, aes(x = yearWeek, y = count, color = source)) +
  geom_line(size = 1) +
  # geom_point(size = 0.5, alpha = 0.6) +
  labs(
    color = "source"
  ) +
  scale_x_date(limits = c(min(as.Date(dat_c$week_date)), max(as.Date(dat_c$week_date))), # 2009 instead of min(as.Date(dat_c$week_date))
               date_breaks = "1 year",
               date_labels = "%Y") +
  theme_bw() +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))




