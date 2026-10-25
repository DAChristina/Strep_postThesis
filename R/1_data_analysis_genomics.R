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
new_dates <- as.numeric(seq(from = as.Date(min(ne_55$date)), to = as.Date(max(ne_55$date)), by = "1 week")) %>% 
  glimpse()

# I change smooth function from natural (cubic spline) to monoH.FC.
fn_Ne    <- splinefun(x = as.Date(ne_55$date), y = ne_55$Ne, method = "monoH.FC")
fn_Ne_lo <- splinefun(x = as.Date(ne_55$date), y = ne_55$Ne_lb, method = "monoH.FC")
fn_Ne_up <- splinefun(x = as.Date(ne_55$date), y = ne_55$Ne_ub, method = "monoH.FC")

spline_fit_Ne    <- fn_Ne(new_dates)
spline_fit_Ne_lo <- fn_Ne_lo(new_dates)
spline_fit_Ne_up <- fn_Ne_up(new_dates)

time_diffs_weeks <- as.numeric(diff(new_dates)) / 7

interpolated_df <- tidyr::tibble(
  date = as.Date(new_dates[-1]),
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
interpolated_df_preEpi <- interpolated_df #%>% 
# dplyr::filter(yearWeek >= as.Date("2000-01-01"))
threshold2 <- stats::mad(interpolated_df_preEpi$rate_itr_Ne, na.rm = TRUE)*0.005
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
    # centred_Ne = ifelse(yearWeek <= as.Date("2010-01-01"), 0, centred_Ne)
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
    count_GPSC55_1 = child_55,
    count_GPSC55_2 = adult_55,
    prop_55_1 = child_55/(child_26 + child_32 + child_55),
    prop_55_2 = adult_55/(adult_26 + adult_32 + adult_55)
  ) %>% 
  glimpse()

# load gen & interpolated_df first, combine with dat_c to calculate GPSC55/12F proportion
age_levels <- c("0-44", "45+")
week_grid_age <- tidyr::expand_grid(
  yearWeek = seq(as.Date("2001-01-01"),
                 max(interpolated_df$yearWeek, na.rm = TRUE), by = "week"), 
  ageGroup12F = age_levels)

gpsc55_age_weekly <- gen %>%
  dplyr::filter(strain == "GPSC55") %>% 
  dplyr::mutate(week_date = as.Date(week_date),
                iso_week = paste0(lubridate::year(week_date), "-W", sprintf("%02d", lubridate::week(week_date)), "-1"),
                yearWeek = ISOweek::ISOweek2date(iso_week)
  ) %>%
  dplyr::count(yearWeek, ageGroup12F, name = "count_GPSC55")

f12_age_weekly <- read.csv("inputs/12F_Jan_2025_combined_cleaned.csv") %>%
  dplyr::mutate(week_date = as.Date(week_date),
                iso_week = paste0(lubridate::year(week_date), "-W", sprintf("%02d", lubridate::week(week_date)), "-1"),
                yearWeek = ISOweek::ISOweek2date(iso_week),
                ageGroup12F = dplyr::case_when(
                  ageGroup6 %in% c("45-64", "65+") ~ "45+",
                  TRUE ~ "0-44")
  ) %>% 
  dplyr::group_by(yearWeek, ageGroup12F) %>%
  dplyr::summarise(count_12F = sum(counts, na.rm = TRUE), .groups = "drop")

test_long <- week_grid_age %>%
  dplyr::left_join(f12_age_weekly
                   ,
                   by = c("yearWeek", "ageGroup12F")
  ) %>%
  dplyr::left_join(gpsc55_age_weekly
                   ,
                   by = c("yearWeek", "ageGroup12F")) %>%
  dplyr::mutate(count_12F = dplyr::coalesce(count_12F, 0),
                count_GPSC55 = dplyr::coalesce(count_GPSC55, 0L),
                week = lubridate::week(yearWeek),
  ) %>%
  # dplyr::select(-count_GPSC55_rand) %>%
  dplyr::left_join(
    interpolated_df %>% 
      dplyr::select(yearWeek, contains("Ne"))
    ,
    by = "yearWeek") %>%
  dplyr::mutate(prop_GPSC55 = ifelse(count_12F > 0, count_GPSC55 / count_12F, 
                                     NA_real_
  ),
  # weekAll = as.numeric((as.Date(yearWeek) - as.Date("2001-01-01")) / 7)
  ) %>%
  dplyr::arrange(ageGroup12F, yearWeek) %>% 
  glimpse()

################################################################################
# the GAM helper function
fit_age_gam <- function(df_age){
  mgcv::gam(cbind(count_GPSC55, count_12F - count_GPSC55) ~
              s(week, bs = "cc", k = 15) +
              s(centred_Ne, bs = "tp", k = 5) +
              s(as.numeric(weekAll), bs = "tp", k = 15)
            ,
            data = df_age %>% 
              dplyr::filter(
                yearWeek >= as.Date("2017-09-01") | yearWeek <= as.Date("2010-01-01"), # using proportions, I set up range when the first time 12F data were recorded ("2001-01-01") instead of GPSC data were intensively collected ("2017-08-01")
                # prop_GPSC55 >= 0 & prop_GPSC55 <= 1 # minor correction for missing 12F data (or GPSC counts > 12F counts)
                prop_GPSC55 <= 1,
                count_12F > 0
              ),
            family = binomial("logit"),
            # weights = count_12F,
            method  = "REML", # ML provides good result, NCV degrades the seasonality, not stable
            knots = list(week = c(0.5, 52.5))
  )
}



# report function
# compile_all_df <- function(df1, df2, df3){
#   
#   
# }
# 
# produce_picts <- function(df1, df2, df3){
#   
#   
# }


################################################################################
dat_age1 <- test_long %>% 
  dplyr::filter(ageGroup12F == "0-44") %>% 
  dplyr::mutate(weekAll = as.numeric((as.Date(yearWeek) - as.Date("2001-01-01"))/7))
dat_age2 <- test_long %>% 
  dplyr::filter(ageGroup12F == "45+") %>% 
  dplyr::mutate(weekAll = as.numeric((as.Date(yearWeek) - as.Date("2001-01-01"))/7))
model_age1 <- fit_age_gam(dat_age1)
model_age2 <- fit_age_gam(dat_age2)

pred_base <- test_long %>%
  dplyr::filter(yearWeek >= as.Date("2010-01-01") | yearWeek <= as.Date("2017-08-01")) %>% 
  dplyr::mutate(weekAll = as.numeric((as.Date(yearWeek) - as.Date("2001-01-01"))/7)) %>% 
  dplyr::select(yearWeek, ageGroup12F, week, weekAll, centred_Ne, count_12F)

pred_age1 <- pred_base %>% 
  dplyr::filter(ageGroup12F == "0-44")
pred_age2 <- pred_base %>% 
  dplyr::filter(ageGroup12F == "45+")
p1 <- predict(model_age1,
              newdata = pred_age1,
              se.fit = TRUE,
              type = "link")
p2 <- predict(model_age2,
              newdata = pred_age2,
              se.fit = TRUE,
              type = "link")
pred_age1 <- pred_age1 %>%
  dplyr::mutate(
    p_hat = as.vector(plogis(p1$fit)),
    p_lo  = as.vector(plogis(p1$fit - 1.96 * p1$se.fit)),
    p_hi  = as.vector(plogis(p1$fit + 1.96 * p1$se.fit)),
    predicted_count_55 = p_hat * count_12F,
    predicted_count_55_lo = p_lo * count_12F,
    predicted_count_55_hi = p_hi * count_12F
  ) %>% 
  glimpse()
pred_age2 <- pred_age2 %>%
  dplyr::mutate(    
    p_hat = as.vector(plogis(p2$fit)),    
    p_lo  = as.vector(plogis(p2$fit - 1.96 * p2$se.fit)),    
    p_hi  = as.vector(plogis(p2$fit + 1.96 * p2$se.fit)),    
    predicted_count_55 = p_hat * count_12F,
    predicted_count_55_lo = p_lo * count_12F,
    predicted_count_55_hi = p_hi * count_12F
  ) %>% 
  glimpse()
pred_total <- dplyr::bind_rows(pred_age1, pred_age2) %>%  
  dplyr::group_by(yearWeek) %>%  
  dplyr::summarise(
    predicted_count_55_total = sum(predicted_count_55, na.rm = TRUE),
    predicted_count_55_total_lo = sum(predicted_count_55_lo, na.rm = TRUE),
    predicted_count_55_total_hi = sum(predicted_count_55_hi, na.rm = TRUE),
    .groups = "drop") %>% 
  glimpse()

# save report
write.csv(pred_age1,
          "raw_data/GPSC55_mlesky_cleaned_interpolated_predictedModel_binom_1.csv",
          row.names = FALSE)
write.csv(pred_age2,
          "raw_data/GPSC55_mlesky_cleaned_interpolated_predictedModel_binom_2.csv",
          row.names = FALSE)
write.csv(pred_total,
          "raw_data/GPSC55_mlesky_cleaned_interpolated_predictedModel_binom.csv",
          row.names = FALSE)


# produce picts
# 1. proportion
dat_fit1 <- ggplot(pred_age1
                   , aes(x = yearWeek, y = p_hat, color = p_hat)) +
  geom_line(size = 1) +
  geom_ribbon(data = pred_age1,
              aes(x = yearWeek,
                  ymin = p_lo,
                  ymax = p_hi),
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
dat_fit1

dat_fit2 <- ggplot(pred_age2
                   , aes(x = yearWeek, y = p_hat, color = p_hat)) +
  geom_line(size = 1) +
  geom_ribbon(data = pred_age2,
              aes(x = yearWeek,
                  ymin = p_lo,
                  ymax = p_hi),
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
dat_fit2

# case count prediction
pred_ribbon1 <- ggplot(pred_age1
                       , aes(x = yearWeek, y = predicted_count_55)) +
  geom_line(size = 0.5) +
  geom_ribbon(data = pred_age1,
              aes(x = yearWeek,
                  ymin = predicted_count_55_lo,
                  ymax = predicted_count_55_hi),
              inherit.aes = FALSE,
              fill = "blue", alpha = 0.2
  ) +
  # geom_point(size = 0.5, alpha = 0.6) +
  scale_x_date(limits = c(as.Date("2010-01-01"), as.Date("2017-11-01")), 
               date_breaks = "1 year",
               date_labels = "%Y") +
  theme_bw() +
  labs(
    title = "GPSC55 Counts Prediction Result (age 0-44)",
    y = "GPSC55 counts prediction"
  ) +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", colour = "transparent"))
pred_ribbon1

pred_ribbon2 <- ggplot(pred_age2
                       , aes(x = yearWeek, y = predicted_count_55)) +
  geom_line(size = 0.5) +
  geom_ribbon(data = pred_age2,
              aes(x = yearWeek,
                  ymin = predicted_count_55_lo,
                  ymax = predicted_count_55_hi),
              inherit.aes = FALSE,
              fill = "blue", alpha = 0.2
  ) +
  # geom_point(size = 0.5, alpha = 0.6) +
  scale_x_date(limits = c(as.Date("2010-01-01"), as.Date("2017-11-01")), 
               date_breaks = "1 year",
               date_labels = "%Y") +
  theme_bw() +
  labs(
    title = "GPSC55 Counts Prediction Result (age 45+)",
    y = "GPSC55 counts prediction"
  ) +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", colour = "transparent"))
pred_ribbon2

# test combine all model viz
reselected_GPSC55_1 <- test_long %>% 
  dplyr::arrange(yearWeek) %>% 
  dplyr::filter(
    ageGroup12F == "0-44",
    yearWeek >= as.Date("2010-01-01"), # using proportions, I set up range when the first time 12F data were recorded ("2001-01-01") instead of GPSC data were intensively collected ("2017-08-01")
  ) %>%
  dplyr::left_join(
    pred_age1 %>% 
      dplyr::select(yearWeek, predicted_count_55)
    ,
    by = "yearWeek",
    relationship = "many-to-many"
  ) %>% 
  dplyr::mutate(
    count_GPSC55 = ifelse(yearWeek <= min(gen$collection_date), NA,
                          count_GPSC55
    ),
  ) %>% 
  dplyr::distinct(yearWeek, .keep_all = TRUE) %>% 
  glimpse()

reselected_1 <- ggplot(reselected_GPSC55_1, aes(x = yearWeek)) +
  geom_line(aes(y = count_12F, colour = "12F")) +
  geom_line(aes(y = count_GPSC55, colour = "GPSC55"), size = 1) +
  geom_line(aes(y = predicted_count_55, colour = "GAM prediction")) +
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
reselected_1

initial_reselected_1 <- ggplot(reselected_GPSC55_1, aes(x = yearWeek)) +
  geom_line(aes(y = count_12F, colour = "12F")) +
  geom_line(aes(y = count_GPSC55, colour = "GPSC55"), size = 1) +
  # geom_line(aes(y = predicted_count_55, colour = "GAM prediction")) +
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
initial_reselected_1

reselected_GPSC55_2 <- test_long %>% 
  dplyr::arrange(yearWeek) %>% 
  dplyr::filter(
    ageGroup12F == "45+",
    yearWeek >= as.Date("2010-01-01"), # using proportions, I set up range when the first time 12F data were recorded ("2001-01-01") instead of GPSC data were intensively collected ("2017-08-01")
  ) %>%
  dplyr::left_join(
    pred_age2 %>% 
      dplyr::select(yearWeek, predicted_count_55)
    ,
    by = "yearWeek",
    relationship = "many-to-many"
  ) %>% 
  dplyr::mutate(
    count_GPSC55 = ifelse(yearWeek <= min(gen$collection_date), NA,
                          count_GPSC55
    ),
  ) %>% 
  dplyr::distinct(yearWeek, .keep_all = TRUE) %>% 
  glimpse()

reselected_2 <- ggplot(reselected_GPSC55_2, aes(x = yearWeek)) +
  geom_line(aes(y = count_12F, colour = "12F")) +
  geom_line(aes(y = count_GPSC55, colour = "GPSC55"), size = 1) +
  geom_line(aes(y = predicted_count_55, colour = "GAM prediction")) +
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
reselected_2

initial_reselected_2 <- ggplot(reselected_GPSC55_2, aes(x = yearWeek)) +
  geom_line(aes(y = count_12F, colour = "12F")) +
  geom_line(aes(y = count_GPSC55, colour = "GPSC55"), size = 1) +
  # geom_line(aes(y = predicted_count_55, colour = "GAM prediction")) +
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
initial_reselected_2

# combined count
reselected_GPSC55_allAge <- test_long %>% 
  dplyr::arrange(yearWeek) %>% 
  dplyr::group_by(yearWeek) %>% 
  dplyr::summarise(
    # count_12F = sum(count_12F, na.rm = TRUE),
    count_GPSC55 = sum(count_GPSC55, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  # dplyr::left_join(
  #   gen %>%
  #     dplyr::filter(strain == "GPSC55") %>% 
  #     dplyr::mutate(week_date = as.Date(week_date),
  #                   iso_week = paste0(lubridate::year(week_date), "-W", sprintf("%02d", lubridate::week(week_date)), "-1"),
  #                   yearWeek = ISOweek::ISOweek2date(iso_week)
  #     ) %>%
  #     dplyr::count(yearWeek, name = "count_GPSC55")
  #   ,
  #   by = "yearWeek"
  # ) %>% 
  dplyr::left_join(
    read.csv("inputs/12F_Jan_2025_combined_cleaned.csv") %>%
      dplyr::mutate(week_date = as.Date(week_date),
                    iso_week = paste0(lubridate::year(week_date), "-W", sprintf("%02d", lubridate::week(week_date)), "-1"),
                    yearWeek = ISOweek::ISOweek2date(iso_week),
                    ageGroup12F = dplyr::case_when(
                      ageGroup6 %in% c("45-64", "65+") ~ "45+",
                      TRUE ~ "0-44")
      ) %>% 
      dplyr::group_by(yearWeek) %>%
      dplyr::summarise(count_12F = sum(counts, na.rm = TRUE), .groups = "drop") %>% 
      dplyr::select(yearWeek, count_12F)
    ,
    by = "yearWeek"
    
  ) %>% 
  dplyr::filter(
    yearWeek >= as.Date("2010-01-01"), # using proportions, I set up range when the first time 12F data were recorded ("2001-01-01") instead of GPSC data were intensively collected ("2017-08-01")
  ) %>%
  dplyr::left_join(
    pred_total %>% 
      dplyr::select(yearWeek, predicted_count_55_total) %>% 
      dplyr::rename(
        predicted_count_55 = predicted_count_55_total,
        # predicted_count_55_lo = predicted_count_55_total_lo,
        # predicted_count_55_hi = predicted_count_55_total_hi
      )
    ,
    by = "yearWeek",
    relationship = "many-to-many"
  ) %>% 
  dplyr::mutate(
    count_GPSC55 = ifelse(yearWeek <= min(gen$collection_date), NA,
                          count_GPSC55
    ),
    count_12F = coalesce(count_12F, 0)
  ) %>% 
  dplyr::distinct(yearWeek, .keep_all = TRUE) %>% 
  glimpse()

reselected_all <- ggplot(reselected_GPSC55_allAge, aes(x = yearWeek)) +
  geom_line(aes(y = count_12F, colour = "12F")) +
  geom_line(aes(y = count_GPSC55, colour = "GPSC55"), size = 1) +
  geom_line(aes(y = predicted_count_55, colour = "GAM prediction")) +
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
  labs(title = "Data for Model Inference (all age group combined)") +
  theme(legend.position = c(0.1, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))
reselected_all


################################################################################
# recollect all picts
for (i in 1:2) {
  
  save_plot <- cowplot::plot_grid(get(paste0("initial_reselected_", i)),
                                  get(paste0("dat_fit", i)),
                                  ncol = 1,
                                  labels = c("A", "B")
  )
  
  ribbon_picts <- cowplot::plot_grid(get(paste0("pred_ribbon", i)),
                                     get(paste0("reselected_", i)),
                                     ncol = 1,
                                     labels = c("C", "D")
  )
  
  combined_plot <- cowplot::plot_grid(save_plot, ribbon_picts,
                                      ncol = 2,
                                      labels = c("", "")
  )
  
  png(filename = paste0("report/picts_proportion_centredNe_allWGS_combined_",
                        i, ".png"),
      width = 40, height = 20, units = "cm", res = 300)
  
  print(combined_plot)
  dev.off()
  
  # additional pics for model adjustment
  
  # par(mfrow = c(2, 2))
  # mgcv::gam.check(get(paste0("model_age",i)))
  # par(mfrow = c(1, 1))
  sim <- DHARMa::simulateResiduals(get(paste0("model_age",i)),
                                   n = 1000)
  
  png(filename = paste0("report/picts_proportion_centredNe_allWGS_combined_1gam.check_",
                        i, ".png"),
      width = 20, height = 20, units = "cm", res = 300)
  par(mfrow = c(2, 2))
  
  DHARMa::plotQQunif(sim)
  DHARMa::plotResiduals(sim)
  disp <- DHARMa::testDispersion(sim)
  zero <- DHARMa::testZeroInflation(sim)  # draws zero-inflation histogram
  
  dev.off()
  
  png(filename = paste0("report/picts_proportion_centredNe_allWGS_combined_2plot_",
                        i, ".png"),
      width = 20, height = 20, units = "cm", res = 300)
  par(mfrow = c(2, 2))
  plot(get(paste0("model_age",i)))
  par(mfrow = c(1, 1))
  dev.off()
  
  # save model result
  check_output <- capture.output(
    mgcv::gam.check(get(paste0("model_age",i)))
  )
  
  writeLines(check_output,
             paste0("report/report_gam_check_",
                    i, ".txt"
                    )
  )
}


################################################################################
# test gratia to compare age groups
# model_by <- mgcv::gam(cbind(count_GPSC55, count_12F - count_GPSC55) ~
#                         ageGroup12F +
#                         s(week, bs = "cc", by = ageGroup12F, k = 15) +
#                         s(centred_Ne, bs = "tp", by = ageGroup12F, k = 5) +
#                         s(weekAll, bs = "tp", by = ageGroup12F, k = 15)
#                       ,
#                       data = test_long %>% 
#                         dplyr::filter(
#                           yearWeek >= as.Date("2017-09-01") | yearWeek <= as.Date("2010-01-01"), # using proportions, I set up range when the first time 12F data were recorded ("2001-01-01") instead of GPSC data were intensively collected ("2017-08-01")
#                           # prop_GPSC55 >= 0 & prop_GPSC55 <= 1 # minor correction for missing 12F data (or GPSC counts > 12F counts)
#                           prop_GPSC55 <= 1,
#                           count_12F > 0
#                         ) %>% 
#                         dplyr::mutate(
#                           weekAll = as.numeric((as.Date(yearWeek) - as.Date("2001-01-01"))/7),
#                           ageGroup12F = factor(ageGroup12F,
#                                                levels = c("0-44", "45+"))
#                           )
#                       ,
#                       family = binomial("logit"),
#                       # weights = count_12F,
#                       method  = "REML", # ML provides good result, NCV degrades the seasonality, not stable
#                       knots = list(week = c(0.5, 52.5))
# )
# gratia::difference_smooths(model_by, smooth = "s(weekAll)")
# gratia::draw(model_by)
# 
# gam.check(model_by)
# mgcv::concurvity(model_by, full = TRUE)

