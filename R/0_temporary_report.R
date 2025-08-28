
################################################################################
selected_GPSC55 <- test %>% 
  dplyr::arrange(yearWeek) %>% 
  dplyr::filter(
    yearWeek >= as.Date("2017-08-01"), # using proportions, I set up range when the first time 12F data were recorded ("2001-01-01") instead of GPSC data were intensively collected ("2017-08-01")
    # !is.na(prop_GPSC55), # no 12F data available after 2020-05-11
    prop_GPSC55 >= 0 & prop_GPSC55 <= 1 # minor correction for missing 12F data (or GPSC counts > 12F counts)
  ) %>% 
  # dplyr::mutate(
  #   count_GPSC55 = case_when(
  #     yearWeek <= as.Date("2017-12-01") & yearWeek >= as.Date("2016-01-01") ~ NA,
  #     TRUE ~ count_GPSC55
  #     ),
  #   count_12F = case_when(
  #     yearWeek <= as.Date("2017-12-01") & yearWeek >= as.Date("2016-03-01") ~ NA,
  #     TRUE ~ count_12F
  #   )
  # ) %>% 
  dplyr::mutate(sin_week = sin(2*pi*lubridate::isoweek(yearWeek)/52),
                cos_week = cos(2*pi*lubridate::isoweek(yearWeek)/52)) %>% 
  glimpse()

dat_model <- ggplot(selected_GPSC55, aes(x = yearWeek)) +
  geom_line(aes(y = count_12F, colour = "12F")) +
  geom_line(aes(y = count_GPSC55, colour = "GPSC55")) +
  geom_line(aes(y = itr_Ne, colour = "Ne"), size = 1.5) +
  geom_line(aes(y = centred_Ne, colour = "Ne (centred)"), size = 1.5) +
  geom_vline(xintercept = as.Date("2016-03-01"), color = "steelblue", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2017-08-01"), color = "steelblue", linetype = "dashed") +
  scale_colour_manual(values = c("12F" = "steelblue",
                                 "GPSC55" = "maroon",
                                 "Ne" = "gold2",
                                 "Ne (centred)" = "orange")) +
  theme_bw() +
  scale_y_log10() +
  scale_x_date() +
  labs(title = "Data for Model Inference") +
  theme(legend.position = c(0.9, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))
dat_model

# test new model based on GPSC55/12F proportion
model_gam_binom_spline <- mgcv::gam(prop_GPSC55 ~ sin_week + cos_week + s(centred_Ne) + as.numeric(yearWeek),
                                    data = selected_GPSC55,
                                    family = binomial("logit"),
                                    weights = count_12F)
model_gam_binom_tensor <- mgcv::gam(prop_GPSC55 ~ sin_week + cos_week + te(centred_Ne) + as.numeric(yearWeek),
                                    data = selected_GPSC55,
                                    family = binomial("logit"),
                                    weights = count_12F)
model_gam_binom_spline_ml <- mgcv::gam(prop_GPSC55 ~ sin_week + cos_week + s(centred_Ne) + as.numeric(yearWeek),
                                    data = selected_GPSC55,
                                    family = binomial("logit"), method = "ML",
                                    weights = count_12F)
model_gam_binom_tensor_ml <- mgcv::gam(prop_GPSC55 ~ sin_week + cos_week + te(centred_Ne) + as.numeric(yearWeek),
                                    data = selected_GPSC55,
                                    family = binomial("logit"), method = "ML",
                                    weights = count_12F)
model_glm_binom <- stats::glm(prop_GPSC55 ~ sin_week + cos_week + centred_Ne + as.numeric(yearWeek),
                              data = selected_GPSC55,
                              family = binomial("logit"),
                              weights = count_12F)

AIC(model_gam_binom_spline, model_gam_binom_tensor,
    model_gam_binom_spline_ml, model_gam_binom_tensor_ml,
    model_glm_binom)
BIC(model_gam_binom_spline, model_gam_binom_tensor,
    model_gam_binom_spline_ml, model_gam_binom_tensor_ml,
    model_glm_binom)
mgcv::gam.check(model_gam_binom_spline)
plot(model_gam_binom_spline)

# test plot for the GLM
par(mfrow = c(2, 2))
plot(model_glm_binom)
par(mfrow = c(1, 1))

# gam is better than glm
saveRDS(model_gam_binom_spline, file = "raw_data/model_Ne_gam_binom.rds")
saveRDS(model_glm_binom, file = "raw_data/model_Ne_glm_binom.rds")


earlier_ne_df <- interpolated_df %>%
  dplyr::filter(yearWeek <= as.Date("2017-08-01")) %>% # midpoint; they started intensively sequenced GPSC55 ("2017-08-01")
  dplyr::mutate(sin_week = sin(2*pi*lubridate::isoweek(yearWeek)/52),
                cos_week = cos(2*pi*lubridate::isoweek(yearWeek)/52)) %>%
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
pred_gam_binom_spline_ml <- predict(model_gam_binom_spline_ml,
                                 newdata = earlier_ne_df,
                                 se.fit = TRUE,
                                 unconditional = TRUE,
                                 type = "link"
)
pred_gam_binom_tensor_ml <- predict(model_gam_binom_tensor_ml,
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
    predicted_prop_GPSC55_gam_binom_spline = plogis(pred_gam_binom_spline$fit),
    predicted_prop_GPSC55_gam_binom_spline_lower = plogis(pred_gam_binom_spline$fit+1.96*pred_gam_binom_spline$se.fit),
    predicted_prop_GPSC55_gam_binom_spline_upper = plogis(pred_gam_binom_spline$fit-1.96*pred_gam_binom_spline$se.fit),
    
    predicted_prop_GPSC55_gam_binom_tensor = plogis(pred_gam_binom_tensor$fit),
    predicted_prop_GPSC55_gam_binom_tensor_lower = plogis(pred_gam_binom_tensor$fit+1.96*pred_gam_binom_tensor$se.fit),
    predicted_prop_GPSC55_gam_binom_tensor_upper = plogis(pred_gam_binom_tensor$fit-1.96*pred_gam_binom_tensor$se.fit),
    
    predicted_prop_GPSC55_gam_binom_spline_ml = plogis(pred_gam_binom_spline_ml$fit),
    predicted_prop_GPSC55_gam_binom_spline_ml_lower = plogis(pred_gam_binom_spline_ml$fit+1.96*pred_gam_binom_spline_ml$se.fit),
    predicted_prop_GPSC55_gam_binom_spline_ml_upper = plogis(pred_gam_binom_spline_ml$fit-1.96*pred_gam_binom_spline_ml$se.fit),
    
    predicted_prop_GPSC55_gam_binom_tensor_ml = plogis(pred_gam_binom_tensor_ml$fit),
    predicted_prop_GPSC55_gam_binom_tensor_ml_lower = plogis(pred_gam_binom_tensor_ml$fit+1.96*pred_gam_binom_tensor_ml$se.fit),
    predicted_prop_GPSC55_gam_binom_tensor_ml_upper = plogis(pred_gam_binom_tensor_ml$fit-1.96*pred_gam_binom_tensor_ml$se.fit),
    
    predicted_prop_GPSC55_glm_binom = plogis(pred_glm_binom$fit),
    predicted_prop_GPSC55_glm_binom_lower = plogis(pred_glm_binom$fit+1.96*pred_glm_binom$se.fit),
    predicted_prop_GPSC55_glm_binom_upper = plogis(pred_glm_binom$fit-1.96*pred_glm_binom$se.fit),
  ) %>% 
  dplyr::left_join(
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
    by = "yearWeek"
  ) %>% 
  dplyr::mutate(
    predicted_count_GPSC55 = predicted_prop_GPSC55_gam_binom_spline*count_12F,
    predicted_count_GPSC55_lo = predicted_prop_GPSC55_gam_binom_spline_lower*count_12F,
    predicted_count_GPSC55_up = predicted_prop_GPSC55_gam_binom_spline_upper*count_12F
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
        source == "predicted_prop_GPSC55_gam_binom_spline" ~ "3.1. Predicted (GAM spline)",
        source == "predicted_prop_GPSC55_gam_binom_tensor" ~ "3.2. Predicted (GAM tensor)",
        source == "predicted_prop_GPSC55_gam_binom_spline_ml" ~ "3.3. Predicted (GAM spline ML method)",
        source == "predicted_prop_GPSC55_gam_binom_tensor_ml" ~ "3.4. Predicted (GAM tensor ML method)",
        source == "predicted_prop_GPSC55_glm_binom" ~ "3.5. Predicted (GLM)",
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
  dplyr::filter(source != "2. Interpolated Ne") %>% # omit Ne
  glimpse()


dat_fit <- ggplot(combined %>% 
                    dplyr::filter(source != "1. Data GPSC55/12F")
                  , aes(x = yearWeek, y = count, color = source)) +
  geom_line(size = 1) +
  geom_ribbon(data = earlier_ne_df,
              aes(x = yearWeek,
                  ymin = predicted_prop_GPSC55_gam_binom_spline_lower,
                  ymax = predicted_prop_GPSC55_gam_binom_spline_upper),
              inherit.aes = FALSE,
              fill = "orange", alpha = 0.2
  ) +
  geom_ribbon(data = earlier_ne_df,
              aes(x = yearWeek,
                  ymin = predicted_prop_GPSC55_gam_binom_tensor_lower,
                  ymax = predicted_prop_GPSC55_gam_binom_tensor_upper),
              inherit.aes = FALSE,
              fill = "darkgreen", alpha = 0.2
  ) +
  geom_ribbon(data = earlier_ne_df,
              aes(x = yearWeek,
                  ymin = predicted_prop_GPSC55_gam_binom_spline_ml_lower,
                  ymax = predicted_prop_GPSC55_gam_binom_spline_ml_upper),
              inherit.aes = FALSE,
              fill = "green", alpha = 0.2
  ) +
  geom_ribbon(data = earlier_ne_df,
              aes(x = yearWeek,
                  ymin = predicted_prop_GPSC55_gam_binom_tensor_ml_lower,
                  ymax = predicted_prop_GPSC55_gam_binom_tensor_ml_upper),
              inherit.aes = FALSE,
              fill = "steelblue", alpha = 0.2
  ) +
  geom_ribbon(data = earlier_ne_df,
              aes(x = yearWeek,
                  ymin = predicted_prop_GPSC55_glm_binom_lower,
                  ymax = predicted_prop_GPSC55_glm_binom_upper),
              inherit.aes = FALSE,
              fill = "purple", alpha = 0.2
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
        legend.background = element_rect(fill = "transparent", colour = "transparent")) #+
  # facet_wrap(~ source)
dat_fit


png("report/picts_proportion_centredNe_allWGS.png",
    width = 24, height = 24, unit = "cm", res = 300)
cowplot::plot_grid(dat_model, dat_fit,
                   ncol = 1,
                   labels = c("A", "B"))
dev.off()

# test GPSC55 case counts prediction
ggplot(earlier_ne_df
       , aes(x = yearWeek, y = predicted_count_GPSC55)) +
  geom_line(size = 0.5) +
  geom_ribbon(data = earlier_ne_df,
              aes(x = yearWeek,
                  ymin = predicted_count_GPSC55_lo,
                  ymax = predicted_count_GPSC55_up),
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
