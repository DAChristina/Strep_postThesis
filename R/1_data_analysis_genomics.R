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

interpolated_df <- tidyr::tibble(
  date = new_dates[-1],
  itr_Ne = spline_fit_Ne[-1],
  itr_Ne_lo = spline_fit_Ne_lo[-1],
  itr_Ne_up = spline_fit_Ne_up[-1],
  time_diffs_weeks = time_diffs_weeks,
) %>%
  dplyr::mutate(
    iso_week = paste0(year(date), "-W", sprintf("%02d", week(date)), "-1"),
    yearWeek = ISOweek::ISOweek2date(iso_week)
  ) %>%
  glimpse()

# write.csv(interpolated_df, "raw_data/GPSC55_mlesky_cleaned_interpolated.csv", row.names = FALSE)

################################################################################
# generate centred Ne
interpolated_df <- interpolated_df %>% 
  dplyr::arrange(yearWeek) %>% 
  dplyr::mutate(
    change_smooth_Ne = abs(itr_Ne - lag(itr_Ne)) / as.numeric(yearWeek - lag(yearWeek)),
    rate_itr_Ne = abs(itr_Ne - lag(itr_Ne)) / as.numeric(yearWeek - lag(yearWeek)),
    # stable_period = rate_itr_Ne < (threshold)
  )
threshold1 <- stats::mad(interpolated_df$rate_itr_Ne, na.rm = TRUE)*1
threshold1

# try threshold before the peak of the epidemic failed, tend to choose 2000 onwards instead
interpolated_df_preEpi <- interpolated_df %>% 
  dplyr::filter(yearWeek >= as.Date("2000-01-01"))
threshold2 <- stats::mad(interpolated_df_preEpi$rate_itr_Ne, na.rm = TRUE)*0.1
threshold2

ggplot(interpolated_df, aes(x = yearWeek, y = change_smooth_Ne)) +
  geom_line() +
  geom_hline(yintercept = threshold1, linetype = "dashed", color = "purple") + 
  geom_hline(yintercept = threshold2, linetype = "dashed", color = "red") + # previously 1.130835e-05
  scale_x_date(breaks = "3 years") +
  scale_y_log10() +
  theme_bw()

plateau_yearWeeks <- interpolated_df$yearWeek[interpolated_df$rate_itr_Ne < threshold1]
plateau_yearWeeks

interpolated_df <- interpolated_df %>%
  dplyr::mutate(
    plateau_flag = ifelse(rate_itr_Ne < threshold1, 1, 0),
    baseline_Ne = ifelse(plateau_flag == 1, itr_Ne, NA)
  ) %>%
  tidyr::fill(baseline_Ne, .direction = "down") %>% 
  dplyr::mutate(
    baseline_Ne = coalesce(baseline_Ne, min(itr_Ne, na.rm = TRUE)),
    centred_Ne = pmax(itr_Ne - baseline_Ne, 0),
    centred_Ne = ifelse(yearWeek <= as.Date("2010-01-01"), 0, centred_Ne)
  ) %>% 
  glimpse()

write.csv(interpolated_df, "raw_data/GPSC55_mlesky_cleaned_interpolated.csv", row.names = FALSE)

png("report/picts_final_centred_Ne.png",
    width = 30, height = 10, unit = "cm", res = 300)
ggplot(interpolated_df, aes(x = yearWeek)) +
  geom_line(aes(y = itr_Ne, color = "Ne")) +
  geom_line(aes(y = centred_Ne, color = "Centred Ne")) +
  scale_colour_manual(values = c("Ne" = "gold2",
                                 "Centred Ne" = "maroon")) +
  theme_bw() +
  scale_x_date(# limits = c(min(as.Date(dat_c$week_date)), max(as.Date(dat_c$week_date))), # 2009 instead of min(as.Date(dat_c$week_date))
    date_breaks = "1 year",
    date_labels = "%Y") +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", colour = "transparent"))
dev.off()


# tbh I have no idea what is this "invasiveness" file is about.
invasiveness <- read.csv("raw_data/gen_lw/invasiveness_12F.csv") %>%
  dplyr::select(X, Adjusted_invasiveness, Age_group) %>% 
  dplyr::mutate(Age_group = tolower(Age_group)) %>% 
  tidyr::pivot_wider(names_from = Age_group, values_from = Adjusted_invasiveness) %>%
  glimpse()

# 28 days of carriage, exp distribution (CI 8.05-103.29 days)
# Lilith's script to specify gamma in terms of mean and variance rather than shape and scale
# convert mean and variance of gamma to shape and scale
mv_to_ss <- function(mean, var) {
  scale <- var / mean
  shape <- mean / scale
  list(shape = shape, scale = scale)
}
# gamma dist functions
qgammamv <- function(p, mean, var) {
  X <- mv_to_ss(mean, var)
  qgamma(p = p, shape = X$shape, scale = X$scale)
}
dgammav <- function(x, mean, var) {
  X <- mv_to_ss(mean, var)
  dgamma(x = x, shape = X$shape, scale = X$scale)
}

# fitting function by least-squares based on mean and CIs
fit_gamma <- function(mean, l, u, max_v = 100, ci = 0.95) {
  
  a <- (1 - ci) / 2
  p <- c(a, 1 - a)
  
  f <- function(v) {
    x <- qgammamv(p = p, mean = mean, var = v)
    sum((x - c(l, u)) ^ 2)
  }
  
  var_D <- optimise(f = f, interval = c(0, max_v), maximum = FALSE)$minimum
  x <- mv_to_ss(mean, var_D)
  qs <- qgamma(p, shape = x$shape, scale = x$scale)
  print(sprintf("fitted qs = (%.3f, %.3f)", qs[1], qs[2]))
  print(sprintf("target qs = (%f, %f)", l, u))
  print(sprintf("fitted var = %.3f", var_D))
  data.frame(dist = "gamma", scale = x$scale, shape = x$shape, mean = mean,
             q2.5 = qs[1], q97.5 = qs[2])
  
}

plot_gamma <- function(shape, scale, xlab, aim = NULL, ci = 0.95, ...) {
  a <- (1 - ci) / 2
  p <- c(a, 1 - a)
  qs <- qgamma(p = p, shape = shape, scale = scale)
  xlim <- c(0, max(qs, aim, na.rm = TRUE) * 1.1)
  curve(dgamma(x, shape = shape, scale = scale),
        xlab = xlab,
        xlim = xlim,
        ylab = "density", ...)
  abline(v = c(shape * scale, qs), lty = c(1, 2, 2), col = "darkred")
  if (!is.null(aim)) {
    abline(v = aim, lty = c(1, 2, 2), col = "darkblue")
  }
}

# test script
qexp(c(0.25, 0.975), 1 / 28)
duration <- fit_gamma(28, 14, 56, max_v = 1e3)
plot_gamma(duration$shape, duration$scale, "Duration of carriage (days)",
           aim = c(28, 14, 56))

# load gen data
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
                # based on non-heterogeneity model fitting,
                # invasiveness differ in <= 44 and 45+.
                ageGroup12F = case_when(
                  Age <= 44 ~ "0-44",
                  Age > 44 ~ "45+"
                ),
                
                collection_date = lubridate::dmy(collection_date),
                month = lubridate::month(collection_date),
                year = lubridate::year(collection_date),
                month = lubridate::month(collection_date),
                yearMonth = as.Date(paste0(format(collection_date, "%Y-%m"), "-01")), # year-month as date
                week_date = lubridate::floor_date(collection_date, "week")
  ) %>% 
  dplyr::select(-ngsid, -Age, -YEAR) %>% 
  glimpse()

write.csv(gen, "raw_data/genomic_data_cleaned.csv", row.names = FALSE)


# I've decided to add data points from GPSC55 sampling period;
# 4 pre-GPSC55 data points sampling period was randomly selected & shared by NC (3 May 2025)
# focused only for GPSC55
rand_ss_counts <- data.frame(
  date = as.Date(c("2001-01-01", 
                   "2001-04-02", "2008-04-01", "2012-12-10", "2014-12-25")), # random date
  child_26 = c(NA, 0, 0, 0, 4),
  child_32 = c(NA, 1, 3, 5, 1),
  child_55 = c(0, 0, 0, 5, 24),
  adult_26 = c(NA, 0, 0, 0, 2),
  adult_32 = c(NA, 3, 3, 2, 1),
  adult_55 = c(0, 0, 0, 0, 10)
) %>% 
  dplyr::transmute(
    iso_week = paste0(year(date), "-W", sprintf("%02d", week(date)), "-1"),
    yearWeek = ISOweek::ISOweek2date(iso_week),
    count_55_1 = child_55,
    count_55_2 = adult_55,
    prop_55_1 = child_55/(child_26 + child_32 + child_55),
    prop_55_2 = adult_55/(adult_26 + adult_32 + adult_55)
  ) %>% 
  glimpse()

# load gen & interpolated_df first, combine with dat_c to calculate GPSC55/12F proportion
test <- gen %>% 
  dplyr::filter(strain == "GPSC55") %>% 
  dplyr::mutate(week_date = as.Date(week_date),
                iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                yearWeek = ISOweek::ISOweek2date(iso_week)
  ) %>% 
  dplyr::group_by(ageGroup12F, yearWeek) %>% 
  dplyr::summarise(count_GPSC55 = n()) %>% 
  dplyr::ungroup() %>% 
  tidyr::pivot_wider(
    .,
    names_from = contains("ageGroup"),
    names_prefix = "count_",
    values_from = "count_GPSC55"
  ) %>% 
  dplyr::rename(
    count_55_1 = "count_0-44",
    count_55_2 = "count_45+"
  ) %>% 
  dplyr::mutate(yearWeek = as.Date(yearWeek)) %>% 
  dplyr::bind_rows(
    rand_ss_counts %>% 
      dplyr::select(yearWeek, count_55_1, count_55_2)
  ) %>% 
  dplyr::full_join(
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
                    yearWeek =ISOweek::ISOweek2date(iso_week),
                    
                    ageGroup12F = case_when(
                      ageGroup6 == "45-64" | ageGroup6 == "65+" ~ "45+",
                      TRUE ~ "0-44"
                    )
      ) %>% 
      dplyr::group_by(yearWeek, ageGroup12F) %>% 
      dplyr::summarise(count_12F = sum(counts)) %>% 
      dplyr::ungroup() %>% 
      tidyr::pivot_wider(
        .,
        names_from = contains("ageGroup"),
        names_prefix = "count_",
        values_from = "count_12F"
      ) %>% 
      dplyr::rename(
        count_12F_1 = "count_0-44",
        count_12F_2 = "count_45+"
      ) %>% 
      dplyr::mutate(yearWeek = as.Date(yearWeek))
    ,
    by = c("yearWeek")
  ) %>% 
  dplyr::mutate(prop_GPSC55_1 = count_55_1/count_12F_1,
                prop_GPSC55_2 = count_55_2/count_12F_2,
                week = week(yearWeek),
  ) %>% 
  dplyr::arrange(yearWeek) %>% 
  glimpse()

################################################################################
# ageGroup 0-44
selected_GPSC55_1 <- test %>% 
  dplyr::mutate(
    year = year(yearWeek),
  ) %>% 
  dplyr::arrange(yearWeek) %>% 
  dplyr::filter(
    yearWeek >= as.Date("2017-09-01") | yearWeek <= as.Date("2016-01-01"), # using proportions, I set up range when the first time 12F data were recorded ("2001-01-01") instead of GPSC data were intensively collected ("2017-08-01")
    prop_GPSC55_1 >= 0 & prop_GPSC55_1 <= 1 # minor correction for missing 12F data (or GPSC counts > 12F counts)
  ) %>% 
  dplyr::mutate(week = week(yearWeek),
                weekAll = as.numeric((as.Date(yearWeek) - as.Date("2001-01-01"))/7),
                sin_week = sin((2*pi*week + 39)/52),
                cos_week = cos((2*pi*week + 39)/52),
                
  ) %>% 
  glimpse()

# quick sinusoidal plot check
ggplot(selected_GPSC55_1, aes(x = yearWeek)) +
  geom_line(aes(y = sin_week)) +
  geom_line(aes(y = cos_week)) +
  theme_bw() +
  # scale_y_log10() +
  scale_x_date(limits = c(as.Date("2017-09-01"), as.Date("2020-01-01")),
               date_breaks = "1 month",
               date_labels = "%m/%Y") +
  labs(title = "Data for Model Inference") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = c(0.2, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))


dat_model_1 <- ggplot(selected_GPSC55_1, aes(x = yearWeek)) +
  geom_line(aes(y = count_12F_1, colour = "12F")) +
  geom_line(aes(y = count_55_1, colour = "GPSC55")) +
  geom_line(aes(y = itr_Ne, colour = "Ne"), size = 1.5) +
  geom_line(aes(y = centred_Ne, colour = "Ne (centred)"), size = 1.5) +
  geom_vline(xintercept = as.Date("2016-01-01"), color = "steelblue", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2017-09-01"), color = "steelblue", linetype = "dashed") +
  scale_colour_manual(values = c("12F" = "steelblue",
                                 "GPSC55" = "maroon",
                                 "Ne" = "gold2",
                                 "Ne (centred)" = "orange")) +
  theme_bw() +
  scale_y_log10() +
  scale_x_date(date_breaks = "1 year",
               date_labels = "%Y") +
  labs(title = "Data for Model Inference (age 0-44)") +
  theme(legend.position = c(0.2, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))
dat_model_1

# test new model based on GPSC55/12F proportion
model_gam_binom_spline <- mgcv::gam(prop_GPSC55_1 ~ sin_week + cos_week + s(centred_Ne) + as.numeric(weekAll),
                                    data = selected_GPSC55_1,
                                    family = binomial("logit"),
                                    weights = count_12F_1)
model_gam_binom_tensor <- mgcv::gam(prop_GPSC55_1 ~ sin_week + cos_week + te(centred_Ne) + as.numeric(weekAll),
                                    data = selected_GPSC55_1,
                                    family = binomial("logit"),
                                    weights = count_12F_1)
model_glm_binom <- stats::glm(prop_GPSC55_1 ~ sin_week + cos_week + centred_Ne + as.numeric(weekAll),
                              data = selected_GPSC55_1,
                              family = binomial("logit"),
                              weights = count_12F_1)

AIC(model_gam_binom_spline, model_gam_binom_tensor, model_glm_binom)
BIC(model_gam_binom_spline, model_gam_binom_tensor, model_glm_binom)
mgcv::gam.check(model_gam_binom_spline)
plot(model_gam_binom_spline)

# test plot for the GLM
par(mfrow = c(2, 2))
plot(model_glm_binom)
par(mfrow = c(1, 1))

# gam is better than glm
# saveRDS(model_gam_binom_spline, file = "raw_data/model_Ne_gam_binom_1.rds")
# saveRDS(model_glm_binom, file = "raw_data/model_Ne_glm_binom_1.rds")


# seasonal effect check
sin_coef <- coef(model_gam_binom_spline)["sin_week"]
cos_coef <- coef(model_gam_binom_spline)["cos_week"]

tibble(
  yearWeek = seq(from = as.Date("2000-01-01"),
                 to = as.Date("2023-01-01"),
                 by = "week")
) %>% 
  dplyr::mutate(
    week = week(yearWeek),
    sin_week = sin((2*pi*week + 39*3)/52), # not sure why time shift should be *3
    cos_week = cos((2*pi*week + 39*3)/52), # not sure why time shift should be *3
    
    seasonal_effect = sin_coef * sin_week + cos_coef * cos_week
    
  ) %>% 
  dplyr::filter(yearWeek >= as.Date("2020-01-01")) %>% 
  ggplot(aes(x = yearWeek, y = seasonal_effect)) +
  geom_line(color = "steelblue", size = 1) +
  labs(
    title = "Estimated Seasonal Wave from the GAM",
    x = "Week",
    y = "Seasonal Effect (sin_coef*sin_week + cos_coef*cos_week)"
  ) +
  scale_x_date(date_breaks = "2 months",
               date_labels = "%m/%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = c(0.2, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))


earlier_ne_df <- interpolated_df %>%
  dplyr::filter(yearWeek >= as.Date("2010-01-01") & yearWeek <= as.Date("2017-08-01")) %>% # midpoint; they started intensively sequenced GPSC55 ("2017-08-01")
  dplyr::mutate(week = week(yearWeek),
                weekAll = as.numeric((as.Date(yearWeek) - as.Date("2001-01-01"))/7),
                sin_week = sin((2*pi*week + 39*3)/52),
                cos_week = cos((2*pi*week + 39*3)/52)) %>%
  glimpse()

# new model version; se extraction failed to load within dplyr::mutate command
pred_gam_binom_spline <- predict(model_gam_binom_spline,
                                 newdata = earlier_ne_df,
                                 se.fit = TRUE,
                                 unconditional = TRUE,
                                 type = "link"
)
pred_gam_binom_tensor <- predict(model_gam_binom_tensor,
                                 newdata = earlier_ne_df,
                                 se.fit = TRUE,
                                 unconditional = TRUE,
                                 type = "link"
)
pred_glm_binom <- predict(model_glm_binom,
                          newdata = earlier_ne_df,
                          se.fit = TRUE,
                          unconditional = TRUE,
                          type = "link"
)

earlier_ne_df <- earlier_ne_df %>%
  dplyr::mutate(
    predicted_prop_GPSC55_1_gam_binom_spline = plogis(pred_gam_binom_spline$fit),
    predicted_prop_GPSC55_1_gam_binom_spline_lower = plogis(pred_gam_binom_spline$fit+1.96*pred_gam_binom_spline$se.fit),
    predicted_prop_GPSC55_1_gam_binom_spline_upper = plogis(pred_gam_binom_spline$fit-1.96*pred_gam_binom_spline$se.fit),
    
    predicted_prop_GPSC55_1_gam_binom_tensor = plogis(pred_gam_binom_tensor$fit),
    predicted_prop_GPSC55_1_gam_binom_tensor_lower = plogis(pred_gam_binom_tensor$fit+1.96*pred_gam_binom_tensor$se.fit),
    predicted_prop_GPSC55_1_gam_binom_tensor_upper = plogis(pred_gam_binom_tensor$fit-1.96*pred_gam_binom_tensor$se.fit),
    
    predicted_prop_GPSC55_1_glm_binom = plogis(pred_glm_binom$fit),
    predicted_prop_GPSC55_1_glm_binom_lower = plogis(pred_glm_binom$fit+1.96*pred_glm_binom$se.fit),
    predicted_prop_GPSC55_1_glm_binom_upper = plogis(pred_glm_binom$fit-1.96*pred_glm_binom$se.fit),
  ) %>% 
  dplyr::left_join(
    read.csv("raw_data/12F_Jan_2025_combined_cleaned.csv") %>% 
      dplyr::mutate(week_date = as.Date(week_date),
                    iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                    yearWeek =ISOweek::ISOweek2date(iso_week)
      ) %>% 
      dplyr::group_by(yearWeek) %>% 
      dplyr::summarise(count_12F_1 = sum(counts)) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(yearWeek = as.Date(yearWeek))
    ,
    by = "yearWeek"
  ) %>% 
  dplyr::mutate(
    predicted_count_55_1 = predicted_prop_GPSC55_1_gam_binom_spline*count_12F_1,
    predicted_count_55_1_lo = predicted_prop_GPSC55_1_gam_binom_spline_lower*count_12F_1,
    predicted_count_55_1_up = predicted_prop_GPSC55_1_gam_binom_spline_upper*count_12F_1
  ) %>%
  glimpse()

write.csv(earlier_ne_df, "raw_data/GPSC55_mlesky_cleaned_interpolated_predictedModel_binom_1.csv", row.names = FALSE)

# test GPSC55 previous WGS
combined <- dplyr::bind_rows(
  test %>% 
    select(yearWeek, prop_GPSC55_1) %>%
    mutate(source = "1. Data GPSC55/12F") %>% 
    rename(count = prop_GPSC55_1)
  ,
  earlier_ne_df %>%
    select(yearWeek, contains("predicted_prop_GPSC55_1"),
    ) %>%
    tidyr::pivot_longer(
      cols = contains("predicted_prop_"),
      names_to = "source",
      values_to = "count"
    ) %>%
    dplyr::mutate(
      source = case_when(
        source == "predicted_prop_GPSC55_1_gam_binom_spline" ~ "3.1. Predicted (GAM spline)",
        source == "predicted_prop_GPSC55_1_gam_binom_tensor" ~ "3.1. Predicted (GAM tensor)",
        source == "predicted_prop_GPSC55_1_glm_binom" ~ "3.3. Predicted (GLM)",
      )
    ) %>%
    unnest(cols = count)
  ,
  interpolated_df %>%
    dplyr::select(yearWeek, itr_Ne) %>% 
    rename(count = itr_Ne) %>% 
    mutate(source = "2. Interpolated Ne")
) %>% 
  # weird array conversion
  dplyr::mutate(
    count = as.data.frame(count)
  ) %>% 
  unnest(cols = count) %>% 
  dplyr::filter(yearWeek >= as.Date("2000-01-01"),
                source != "2. Interpolated Ne") %>% # omit Ne
  glimpse()

dat_fit <- ggplot(combined %>% 
                    dplyr::filter(source != "1. Data GPSC55/12F")
                  , aes(x = yearWeek, y = count, color = source)) +
  geom_line(size = 1) +
  geom_ribbon(data = earlier_ne_df,
              aes(x = yearWeek,
                  ymin = predicted_prop_GPSC55_1_gam_binom_spline_lower,
                  ymax = predicted_prop_GPSC55_1_gam_binom_spline_upper),
              inherit.aes = FALSE,
              fill = "orange", alpha = 0.2
  ) +
  geom_ribbon(data = earlier_ne_df,
              aes(x = yearWeek,
                  ymin = predicted_prop_GPSC55_1_gam_binom_tensor_lower,
                  ymax = predicted_prop_GPSC55_1_gam_binom_tensor_upper),
              inherit.aes = FALSE,
              fill = "green", alpha = 0.2
  ) +
  geom_ribbon(data = earlier_ne_df,
              aes(x = yearWeek,
                  ymin = predicted_prop_GPSC55_1_glm_binom_lower,
                  ymax = predicted_prop_GPSC55_1_glm_binom_upper),
              inherit.aes = FALSE,
              fill = "steelblue", alpha = 0.2
  ) +
  scale_x_date(
    # limits = c(min(as.Date(dat_c$week_date)), max(as.Date(dat_c$week_date))), # 2009 instead of min(as.Date(dat_c$week_date))
    date_breaks = "1 year",
    date_labels = "%m/%Y") +
  theme_bw() +
  labs(
    title = "Proportion Result",
    colour = "source",
    y = "proportion (GPSC55/12F)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", colour = "transparent"))
dat_fit

png("report/picts_proportion_centredNe_allWGS_1.png",
    width = 24, height = 24, unit = "cm", res = 300)
cowplot::plot_grid(dat_model_1, dat_fit,
                   ncol = 1,
                   labels = c("A", "B"))
dev.off()

# test GPSC55 case counts prediction
ggplot(earlier_ne_df
       , aes(x = yearWeek, y = predicted_count_55_1)) +
  geom_line(size = 0.5) +
  geom_ribbon(data = earlier_ne_df,
              aes(x = yearWeek,
                  ymin = predicted_count_55_1_lo,
                  ymax = predicted_count_55_1_up),
              inherit.aes = FALSE,
              fill = "blue", alpha = 0.2
  ) +
  # geom_point(size = 0.5, alpha = 0.6) +
  scale_x_date(limits = c(as.Date("2001-01-01"), as.Date("2018-01-01")), 
               date_breaks = "1 year",
               date_labels = "%Y") +
  theme_bw() +
  labs(
    title = "GPSC55 Counts Prediction Result",
    y = "GPSC55 counts prediction"
  ) +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", colour = "transparent"))


# test combine all model viz
reselected_GPSC55_1 <- test %>% 
  dplyr::arrange(yearWeek) %>% 
  dplyr::filter(
    yearWeek >= as.Date("2010-01-01"), # using proportions, I set up range when the first time 12F data were recorded ("2001-01-01") instead of GPSC data were intensively collected ("2017-08-01")
  ) %>% 
  dplyr::left_join(
    earlier_ne_df %>% 
      dplyr::select(yearWeek, predicted_count_55_1)
    ,
    by = "yearWeek",
    relationship = "many-to-many"
  ) %>% 
  glimpse()

ggplot(reselected_GPSC55_1, aes(x = yearWeek)) +
  geom_line(aes(y = count_12F_1, colour = "12F")) +
  geom_line(aes(y = count_55_1, colour = "GPSC55"), size = 1) +
  geom_line(aes(y = predicted_count_55_1, colour = "GAM prediction")) +
  geom_vline(xintercept = as.Date("2017-09-01"), color = "steelblue", linetype = "dashed") +
  scale_colour_manual(values = c("12F" = "steelblue",
                                 "GPSC55" = "maroon",
                                 "GAM prediction" = "violet",
                                 "Ne" = "gold2",
                                 "Ne (centred)" = "orange")) +
  theme_bw() +
  # scale_y_log10() +
  scale_x_date(date_breaks = "1 year",
               date_labels = "%Y") +
  labs(title = "Data for Model Inference (age 0-44)") +
  theme(legend.position = c(0.1, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))


################################################################################
# ageGroup 45+
selected_GPSC55_2 <- test %>% 
  dplyr::mutate(
    year = year(yearWeek),
  ) %>% 
  dplyr::arrange(yearWeek) %>% 
  dplyr::filter(
    yearWeek >= as.Date("2017-09-01") | yearWeek <= as.Date("2016-01-01"), # using proportions, I set up range when the first time 12F data were recorded ("2001-01-01") instead of GPSC data were intensively collected ("2017-08-01")
    prop_GPSC55_2 >= 0 & prop_GPSC55_2 <= 1 # minor correction for missing 12F data (or GPSC counts > 12F counts)
  ) %>% 
  dplyr::mutate(week = week(yearWeek),
                weekAll = as.numeric((as.Date(yearWeek) - as.Date("2001-01-01"))/7),
                sin_week = sin((2*pi*week + 39)/52),
                cos_week = cos((2*pi*week + 39)/52),
                
  ) %>% 
  glimpse()

# quick sinusoidal plot check
ggplot(selected_GPSC55_2, aes(x = yearWeek)) +
  geom_line(aes(y = sin_week)) +
  geom_line(aes(y = cos_week)) +
  theme_bw() +
  # scale_y_log10() +
  scale_x_date(limits = c(as.Date("2017-09-01"), as.Date("2020-01-01")),
               date_breaks = "1 month",
               date_labels = "%m/%Y") +
  labs(title = "Data for Model Inference") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = c(0.2, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))


dat_model_2 <- ggplot(selected_GPSC55_2, aes(x = yearWeek)) +
  geom_line(aes(y = count_12F_2, colour = "12F")) +
  geom_line(aes(y = count_55_2, colour = "GPSC55")) +
  geom_line(aes(y = itr_Ne, colour = "Ne"), size = 1.5) +
  geom_line(aes(y = centred_Ne, colour = "Ne (centred)"), size = 1.5) +
  geom_vline(xintercept = as.Date("2016-01-01"), color = "steelblue", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2017-09-01"), color = "steelblue", linetype = "dashed") +
  scale_colour_manual(values = c("12F" = "steelblue",
                                 "GPSC55" = "maroon",
                                 "Ne" = "gold2",
                                 "Ne (centred)" = "orange")) +
  theme_bw() +
  scale_y_log10() +
  scale_x_date(date_breaks = "1 year",
               date_labels = "%Y") +
  labs(title = "Data for Model Inference (age 45+)") +
  theme(legend.position = c(0.2, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))
dat_model_2

# test new model based on GPSC55/12F proportion
model_gam_binom_spline <- mgcv::gam(prop_GPSC55_2 ~ sin_week + cos_week + s(centred_Ne) + as.numeric(weekAll),
                                    data = selected_GPSC55_2,
                                    family = binomial("logit"),
                                    weights = count_12F_2)
model_gam_binom_tensor <- mgcv::gam(prop_GPSC55_2 ~ sin_week + cos_week + te(centred_Ne) + as.numeric(weekAll),
                                    data = selected_GPSC55_2,
                                    family = binomial("logit"),
                                    weights = count_12F_2)
model_glm_binom <- stats::glm(prop_GPSC55_2 ~ sin_week + cos_week + centred_Ne + as.numeric(weekAll),
                              data = selected_GPSC55_2,
                              family = binomial("logit"),
                              weights = count_12F_2)

AIC(model_gam_binom_spline, model_gam_binom_tensor, model_glm_binom)
BIC(model_gam_binom_spline, model_gam_binom_tensor, model_glm_binom)
mgcv::gam.check(model_gam_binom_spline)
plot(model_gam_binom_spline)

# test plot for the GLM
par(mfrow = c(2, 2))
plot(model_glm_binom)
par(mfrow = c(1, 1))

# gam is better than glm
# saveRDS(model_gam_binom_spline, file = "raw_data/model_Ne_gam_binom_2.rds")
# saveRDS(model_glm_binom, file = "raw_data/model_Ne_glm_binom_2.rds")


# seasonal effect check
sin_coef <- coef(model_gam_binom_spline)["sin_week"]
cos_coef <- coef(model_gam_binom_spline)["cos_week"]

tibble(
  yearWeek = seq(from = as.Date("2000-01-01"),
                 to = as.Date("2023-01-01"),
                 by = "week")
) %>% 
  dplyr::mutate(
    week = week(yearWeek),
    sin_week = sin((2*pi*week + 39*3)/52), # not sure why time shift should be *3
    cos_week = cos((2*pi*week + 39*3)/52), # not sure why time shift should be *3
    
    seasonal_effect = sin_coef * sin_week + cos_coef * cos_week
    
  ) %>% 
  dplyr::filter(yearWeek >= as.Date("2020-01-01")) %>% 
  ggplot(aes(x = yearWeek, y = seasonal_effect)) +
  geom_line(color = "steelblue", size = 1) +
  labs(
    title = "Estimated Seasonal Wave from the GAM",
    x = "Week",
    y = "Seasonal Effect (sin_coef*sin_week + cos_coef*cos_week)"
  ) +
  scale_x_date(date_breaks = "2 months",
               date_labels = "%m/%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = c(0.2, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))


earlier_ne_df <- interpolated_df %>%
  dplyr::filter(yearWeek >= as.Date("2010-01-01") & yearWeek <= as.Date("2017-08-01")) %>% # midpoint; they started intensively sequenced GPSC55 ("2017-08-01")
  dplyr::mutate(week = week(yearWeek),
                weekAll = as.numeric((as.Date(yearWeek) - as.Date("2001-01-01"))/7),
                sin_week = sin((2*pi*week + 39*3)/52),
                cos_week = cos((2*pi*week + 39*3)/52)) %>%
  glimpse()

# new model version; se extraction failed to load within dplyr::mutate command
pred_gam_binom_spline <- predict(model_gam_binom_spline,
                                 newdata = earlier_ne_df,
                                 se.fit = TRUE,
                                 unconditional = TRUE,
                                 type = "link"
)
pred_gam_binom_tensor <- predict(model_gam_binom_tensor,
                                 newdata = earlier_ne_df,
                                 se.fit = TRUE,
                                 unconditional = TRUE,
                                 type = "link"
)
pred_glm_binom <- predict(model_glm_binom,
                          newdata = earlier_ne_df,
                          se.fit = TRUE,
                          unconditional = TRUE,
                          type = "link"
)

earlier_ne_df <- earlier_ne_df %>%
  dplyr::mutate(
    predicted_prop_GPSC55_2_gam_binom_spline = plogis(pred_gam_binom_spline$fit),
    predicted_prop_GPSC55_2_gam_binom_spline_lower = plogis(pred_gam_binom_spline$fit+1.96*pred_gam_binom_spline$se.fit),
    predicted_prop_GPSC55_2_gam_binom_spline_upper = plogis(pred_gam_binom_spline$fit-1.96*pred_gam_binom_spline$se.fit),
    
    predicted_prop_GPSC55_2_gam_binom_tensor = plogis(pred_gam_binom_tensor$fit),
    predicted_prop_GPSC55_2_gam_binom_tensor_lower = plogis(pred_gam_binom_tensor$fit+1.96*pred_gam_binom_tensor$se.fit),
    predicted_prop_GPSC55_2_gam_binom_tensor_upper = plogis(pred_gam_binom_tensor$fit-1.96*pred_gam_binom_tensor$se.fit),
    
    predicted_prop_GPSC55_2_glm_binom = plogis(pred_glm_binom$fit),
    predicted_prop_GPSC55_2_glm_binom_lower = plogis(pred_glm_binom$fit+1.96*pred_glm_binom$se.fit),
    predicted_prop_GPSC55_2_glm_binom_upper = plogis(pred_glm_binom$fit-1.96*pred_glm_binom$se.fit),
  ) %>% 
  dplyr::left_join(
    read.csv("raw_data/12F_Jan_2025_combined_cleaned.csv") %>% 
      dplyr::mutate(week_date = as.Date(week_date),
                    iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                    yearWeek =ISOweek::ISOweek2date(iso_week)
      ) %>% 
      dplyr::group_by(yearWeek) %>% 
      dplyr::summarise(count_12F_2 = sum(counts)) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(yearWeek = as.Date(yearWeek))
    ,
    by = "yearWeek"
  ) %>% 
  dplyr::mutate(
    predicted_count_55_2 = predicted_prop_GPSC55_2_gam_binom_spline*count_12F_2,
    predicted_count_55_2_lo = predicted_prop_GPSC55_2_gam_binom_spline_lower*count_12F_2,
    predicted_count_55_2_up = predicted_prop_GPSC55_2_gam_binom_spline_upper*count_12F_2
  ) %>%
  glimpse()

write.csv(earlier_ne_df, "raw_data/GPSC55_mlesky_cleaned_interpolated_predictedModel_binom_2.csv", row.names = FALSE)

# test GPSC55 previous WGS
combined <- dplyr::bind_rows(
  test %>% 
    select(yearWeek, prop_GPSC55_2) %>%
    mutate(source = "1. Data GPSC55/12F") %>% 
    rename(count = prop_GPSC55_2)
  ,
  earlier_ne_df %>%
    select(yearWeek, contains("predicted_prop_GPSC55_2"),
    ) %>%
    tidyr::pivot_longer(
      cols = contains("predicted_prop_"),
      names_to = "source",
      values_to = "count"
    ) %>%
    dplyr::mutate(
      source = case_when(
        source == "predicted_prop_GPSC55_2_gam_binom_spline" ~ "3.1. Predicted (GAM spline)",
        source == "predicted_prop_GPSC55_2_gam_binom_tensor" ~ "3.1. Predicted (GAM tensor)",
        source == "predicted_prop_GPSC55_2_glm_binom" ~ "3.3. Predicted (GLM)",
      )
    ) %>%
    unnest(cols = count)
  ,
  interpolated_df %>%
    dplyr::select(yearWeek, itr_Ne) %>% 
    rename(count = itr_Ne) %>% 
    mutate(source = "2. Interpolated Ne")
) %>% 
  # weird array conversion
  dplyr::mutate(
    count = as.data.frame(count)
  ) %>% 
  unnest(cols = count) %>% 
  dplyr::filter(yearWeek >= as.Date("2000-01-01"),
                source != "2. Interpolated Ne") %>% # omit Ne
  glimpse()

dat_fit <- ggplot(combined %>% 
                    dplyr::filter(source != "1. Data GPSC55/12F")
                  , aes(x = yearWeek, y = count, color = source)) +
  geom_line(size = 1) +
  geom_ribbon(data = earlier_ne_df,
              aes(x = yearWeek,
                  ymin = predicted_prop_GPSC55_2_gam_binom_spline_lower,
                  ymax = predicted_prop_GPSC55_2_gam_binom_spline_upper),
              inherit.aes = FALSE,
              fill = "orange", alpha = 0.2
  ) +
  geom_ribbon(data = earlier_ne_df,
              aes(x = yearWeek,
                  ymin = predicted_prop_GPSC55_2_gam_binom_tensor_lower,
                  ymax = predicted_prop_GPSC55_2_gam_binom_tensor_upper),
              inherit.aes = FALSE,
              fill = "green", alpha = 0.2
  ) +
  geom_ribbon(data = earlier_ne_df,
              aes(x = yearWeek,
                  ymin = predicted_prop_GPSC55_2_glm_binom_lower,
                  ymax = predicted_prop_GPSC55_2_glm_binom_upper),
              inherit.aes = FALSE,
              fill = "steelblue", alpha = 0.2
  ) +
  scale_x_date(
    # limits = c(min(as.Date(dat_c$week_date)), max(as.Date(dat_c$week_date))), # 2009 instead of min(as.Date(dat_c$week_date))
    date_breaks = "1 year",
    date_labels = "%m/%Y") +
  theme_bw() +
  labs(
    title = "Proportion Result",
    colour = "source",
    y = "proportion (GPSC55/12F)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", colour = "transparent"))
dat_fit

png("report/picts_proportion_centredNe_allWGS_2.png",
    width = 24, height = 24, unit = "cm", res = 300)
cowplot::plot_grid(dat_model_2, dat_fit,
                   ncol = 1,
                   labels = c("A", "B"))
dev.off()

# test GPSC55 case counts prediction
ggplot(earlier_ne_df
       , aes(x = yearWeek, y = predicted_count_55_2)) +
  geom_line(size = 0.5) +
  geom_ribbon(data = earlier_ne_df,
              aes(x = yearWeek,
                  ymin = predicted_count_55_2_lo,
                  ymax = predicted_count_55_2_up),
              inherit.aes = FALSE,
              fill = "blue", alpha = 0.2
  ) +
  # geom_point(size = 0.5, alpha = 0.6) +
  scale_x_date(limits = c(as.Date("2001-01-01"), as.Date("2018-01-01")), 
               date_breaks = "1 year",
               date_labels = "%Y") +
  theme_bw() +
  labs(
    title = "GPSC55 Counts Prediction Result",
    y = "GPSC55 counts prediction"
  ) +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", colour = "transparent"))


# test combine all model viz
reselected_GPSC55_2 <- test %>% 
  dplyr::arrange(yearWeek) %>% 
  dplyr::filter(
    yearWeek >= as.Date("2010-01-01"), # using proportions, I set up range when the first time 12F data were recorded ("2001-01-01") instead of GPSC data were intensively collected ("2017-08-01")
  ) %>% 
  dplyr::left_join(
    earlier_ne_df %>% 
      dplyr::select(yearWeek, predicted_count_55_2)
    ,
    by = "yearWeek",
    relationship = "many-to-many"
  ) %>% 
  glimpse()

ggplot(reselected_GPSC55_2, aes(x = yearWeek)) +
  geom_line(aes(y = count_12F_2, colour = "12F")) +
  geom_line(aes(y = count_55_2, colour = "GPSC55"), size = 1) +
  geom_line(aes(y = predicted_count_55_2, colour = "GAM prediction")) +
  geom_vline(xintercept = as.Date("2017-09-01"), color = "steelblue", linetype = "dashed") +
  scale_colour_manual(values = c("12F" = "steelblue",
                                 "GPSC55" = "maroon",
                                 "GAM prediction" = "violet",
                                 "Ne" = "gold2",
                                 "Ne (centred)" = "orange")) +
  theme_bw() +
  # scale_y_log10() +
  scale_x_date(date_breaks = "1 year",
               date_labels = "%Y") +
  labs(title = "Data for Model Inference (age 45+)") +
  theme(legend.position = c(0.1, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))


