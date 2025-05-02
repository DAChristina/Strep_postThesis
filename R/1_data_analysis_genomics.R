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


# load gen & interpolated_df first, combine with dat_c to calculate GPSC55/12F proportion
test <- gen %>% 
  dplyr::filter(strain == "GPSC55") %>% 
  dplyr::mutate(week_date = as.Date(week_date),
                iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                yearWeek =ISOweek::ISOweek2date(iso_week)
  ) %>% 
  dplyr::group_by(yearWeek) %>% 
  dplyr::summarise(count_GPSC55 = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(yearWeek = as.Date(yearWeek)) %>% 
  dplyr::left_join(
    interpolated_df %>%
      dplyr::select(yearWeek, contains("Ne")) %>% 
      dplyr::mutate(yearWeek = as.Date(yearWeek))
    ,
    by = c("yearWeek")
  ) %>% 
  dplyr::full_join(
    read.csv("raw_data/12F_Jan_2025_combined_cleaned.csv") %>% 
      dplyr::mutate(week_date = as.Date(week_date),
                    iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                    yearWeek =ISOweek::ISOweek2date(iso_week)
      ) %>% 
      dplyr::group_by(yearWeek) %>% 
      dplyr::summarise(count_12F = sum(counts)) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(yearWeek = as.Date(yearWeek))
    ,
    by = c("yearWeek")
  ) %>% 
  dplyr::mutate(prop_GPSC55 = count_GPSC55/count_12F,
                # prop_GPSC55 = case_when( # test giving prop = 0 for 2001-2010 to catch a grip for the model
                  # yearWeek >= "2001-01-01" & yearWeek <= "2010-01-01" ~ 0,
                  # is.na(prop_GPSC55) ~ 0,
                  # TRUE ~ prop_GPSC55
                # )
  ) %>% 
  glimpse()

################################################################################
selected_GPSC55 <- test %>% 
  dplyr::filter(yearWeek >= as.Date("2001-01-01"), # using proportions, I set up range when the first time 12F data were recorded ("2001-01-01") instead of GPSC data were intensively collected ("2017-08-01")
                !is.na(prop_GPSC55), # no 12F data available after 2020-05-11
                prop_GPSC55 >= 0 & prop_GPSC55 <= 1 # minor correction for missing 12F data (or GPSC counts > 12F counts)
                ) %>% 
  dplyr::mutate(
    count_GPSC55 = case_when(
      yearWeek <= as.Date("2017-12-01") & yearWeek >= as.Date("2016-01-01") ~ NA,
      TRUE ~ count_GPSC55
      ),
    count_12F = case_when(
      yearWeek <= as.Date("2017-12-01") & yearWeek >= as.Date("2016-03-01") ~ NA,
      TRUE ~ count_12F
    )
  ) %>% 
  dplyr::mutate(sin_week = sin(2*pi*lubridate::isoweek(yearWeek)/52),
                cos_week = cos(2*pi*lubridate::isoweek(yearWeek)/52)) %>% 
  glimpse()

dat_model <- ggplot(selected_GPSC55, aes(x = yearWeek)) +
  geom_line(aes(y = count_12F, colour = "12F")) +
  geom_line(aes(y = count_GPSC55, colour = "GPSC55")) +
  geom_vline(xintercept = as.Date("2016-03-01"), color = "steelblue", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2017-12-01"), color = "steelblue", linetype = "dashed") +
  scale_colour_manual(values = c("12F" = "steelblue", "GPSC55" = "maroon")) +
  theme_bw() +
  labs(title = "Data for Model Inference") +
  theme(legend.position = c(0.9, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))
dat_model

# test new model based on GPSC55/12F proportion
model_gam_binom <- mgcv::gam(prop_GPSC55 ~ sin_week + cos_week + yearWeek + s(change_Ne),
                       data = selected_GPSC55,
                       family = binomial("logit"),
                       weights = count_12F)
model_glm_binom <- stats::glm(prop_GPSC55 ~ sin_week + cos_week + yearWeek + change_Ne,
                              data = selected_GPSC55,
                              family = binomial("logit"),
                              weights = count_12F)

AIC(model_gam_binom, model_glm_binom)
BIC(model_gam_binom, model_glm_binom)
mgcv::gam.check(model_gam_binom)
plot(model_gam_binom)

# gam is better than glm
saveRDS(model_gam_binom, file = "raw_data/model_Ne_gam_binom.rds")
saveRDS(model_glm_binom, file = "raw_data/model_Ne_glm_binom.rds")


earlier_ne_df <- interpolated_df %>%
  dplyr::filter(yearWeek <= as.Date("2017-08-01")) %>% # midpoint; they started intensively sequenced GPSC55 ("2017-08-01")
  dplyr::mutate(sin_week = sin(2*pi*lubridate::isoweek(yearWeek)/52),
                cos_week = cos(2*pi*lubridate::isoweek(yearWeek)/52)) %>%
  glimpse()

# new model version; se extraction failed to load within dplyr::mutate command
pred_gam_binom <- predict(model_gam_binom,
                          newdata = earlier_ne_df,
                          se.fit = TRUE,
                          type = "link"
)
pred_glm_binom <- predict(model_glm_binom,
                          newdata = earlier_ne_df,
                          se.fit = TRUE,
                          type = "link"
)

earlier_ne_df <- earlier_ne_df %>%
  dplyr::mutate(
    predicted_prop_GPSC55_gam_binom = plogis(pred_gam_binom$fit),
    predicted_prop_GPSC55_gam_binom_lower = plogis(pred_gam_binom$fit+1.96*pred_gam_binom$se.fit),
    predicted_prop_GPSC55_gam_binom_upper = plogis(pred_gam_binom$fit-1.96*pred_gam_binom$se.fit),
    
    predicted_prop_GPSC55_glm_binom = plogis(pred_glm_binom$fit),
    predicted_prop_GPSC55_glm_binom_lower = plogis(pred_glm_binom$fit+1.96*pred_glm_binom$se.fit),
    predicted_prop_GPSC55_glm_binom_upper = plogis(pred_glm_binom$fit-1.96*pred_glm_binom$se.fit),
  ) %>% 
  glimpse()

write.csv(earlier_ne_df, "raw_data/GPSC55_mlesky_cleaned_interpolated_predictedModel_binom.csv", row.names = FALSE)

# test GPSC55 previous WGS
combined <- dplyr::bind_rows(
  test %>% 
    select(yearWeek, prop_GPSC55) %>%
    mutate(source = "1. Data GPSC55/12F") %>% 
    rename(count = prop_GPSC55)
  ,
  earlier_ne_df %>%
    select(yearWeek, contains("predicted_prop_GPSC55"),
           ) %>%
    tidyr::pivot_longer(
      cols = contains("predicted_prop_"),
      names_to = "source",
      values_to = "count"
    ) %>%
    dplyr::mutate(
      source = case_when(
        source == "predicted_prop_GPSC55_gam_binom" ~ "3.1. Predicted (GAM)",
        source == "predicted_prop_GPSC55_glm_binom" ~ "3.3. Predicted (GLM)",
      )
    ) %>%
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
  dplyr::filter(source != "2. Interpolated Ne") %>% # omit Ne
  glimpse()


dat_fit <- ggplot(combined %>% 
                    dplyr::filter(source != "1. Data GPSC55/12F")
                  , aes(x = yearWeek, y = count, color = source)) +
  geom_line(size = 1) +
  geom_ribbon(data = earlier_ne_df,
              aes(x = yearWeek,
                  ymin = predicted_prop_GPSC55_gam_binom_lower,
                  ymax = predicted_prop_GPSC55_gam_binom_upper),
              inherit.aes = FALSE,
              fill = "darkgreen", alpha = 0.3
  ) +
  geom_ribbon(data = earlier_ne_df,
              aes(x = yearWeek,
                  ymin = predicted_prop_GPSC55_glm_binom_lower,
                  ymax = predicted_prop_GPSC55_glm_binom_upper),
              inherit.aes = FALSE,
              fill = "steelblue", alpha = 0.3
  ) +
  # geom_point(size = 0.5, alpha = 0.6) +
  scale_x_date(# limits = c(min(as.Date(dat_c$week_date)), max(as.Date(dat_c$week_date))), # 2009 instead of min(as.Date(dat_c$week_date))
    date_breaks = "1 year",
    date_labels = "%Y") +
  theme_bw() +
  labs(
    title = "Proportion Result",
    colour = "source",
    y = "proportion (GPSC55/12F)"
  ) +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", colour = "transparent"))
dat_fit


png("report/picts_proportion_modifiedWGS.png",
    width = 24, height = 24, unit = "cm", res = 300)
cowplot::plot_grid(dat_model, dat_fit,
                   ncol = 1,
                   labels = c("A", "B"))
dev.off()



