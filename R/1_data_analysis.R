# Data preparation
library(tidyverse)

## 1. Data Viz and Analysis! ###################################################
# New updated data with meningitis (25.04.2024)
# All df are stored in raw_data
dat <- readxl::read_excel("raw_data/12F_for_NickC_MSc_with_Sample_info_for_DC.xlsx") %>% 
  dplyr::mutate(collection_date = as.Date(collection_date),
                year = year(collection_date),
                month = month(collection_date),
                yearMonth = as.Date(paste0(format(collection_date, "%Y-%m"), "-01")), # year-month as date
                vacc = case_when(
                  year < 2006 ~ "Pre-PCV7",
                  year >= 2006 & year < 2011 ~ "PCV7",
                  year >= 2011 ~ "PCV13",
                  TRUE ~ NA_character_
                ),
                ageGroup2 = case_when(
                  Age < 15 ~ "children",
                  Age >= 15 ~ "adults",
                  is.na(Age) ~ "Unknown"
                ),
                ageGroup3 = case_when(
                  Age < 2 ~ "<2",
                  Age >= 2 & Age < 65 ~ "2-64",
                  Age >= 65 ~ "65+",
                  is.na(Age) ~ "Unknown"
                ),
                ageGroup5 = case_when( # edit 5 age bands
                  Age < 5 ~ "<5",
                  Age >= 5 & Age < 19 ~ "5-18",
                  Age >= 19 & Age < 31 ~ "19-30",
                  Age >= 31 & Age < 65 ~ "31-64",
                  Age >= 65 ~ "65+",
                  is.na(Age) ~ "Unknown"
                ),
                ageGroup6 = case_when( # UKHSA
                  Age < 2 ~ "<2",
                  Age >= 2 & Age < 5 ~ "2-4",
                  Age >= 5 & Age < 15 ~ "5-14",
                  Age >= 15 & Age < 45 ~ "15-44",
                  Age >= 45 & Age < 65 ~ "45-64",
                  Age >= 65 ~ "65+",
                  is.na(Age) ~ "Unknown"
                ),
                ageGroup7 = case_when(
                  Age < 2 ~ "<2",
                  Age >= 2 & Age < 5 ~ "2-4",
                  Age >= 5 & Age < 15 ~ "5-14",
                  Age >= 15 & Age < 31 ~ "15-30", # Edit the Age-band into 15-30 & 31-44
                  Age >= 31 & Age < 45 ~ "31-44", # Edit the Age-band into 15-30 & 31-44
                  Age >= 45 & Age < 65 ~ "45-64",
                  Age >= 65 ~ "65+",
                  is.na(Age) ~ "Unknown"
                )
  )

# Save data to inputs
write.csv(dat, "raw_data/12F_cleaned.csv", row.names = FALSE)

dat_v2 <- readxl::read_excel("raw_data/12F_v2_DC_edit.xlsx") %>% 
  dplyr::mutate(year = as.numeric(paste0(substr(epiyear, 1, 4), ".5")),
                ageGroup6 = gsub("to", "-", NewGraphAgegroup))

write.csv(dat_v2, "raw_data/12F_v2_cleaned.csv", row.names = FALSE)

pop <- readxl::read_excel("raw_data/nomis_2024_10_17_DCedit.xlsx") %>%  # ver.2 2001-2023
  dplyr::mutate(`2024` = `2023`) # Temporary for 2024 population; ONS hasn't released the data yet!


pop_l <- pop %>% 
  tidyr::pivot_longer(cols = `2001`:`2024`,
                      names_to = "year",
                      values_to = "PopSize") %>% 
  dplyr::mutate(Age = gsub("Age ", "", Age),
                Age = ifelse(Age == "Aged 90+", 90, as.numeric(Age) # For incidence calculation, data grouped for people aged 90+
                             ),
                ageGroup2 = case_when(
                  Age < 15 ~ "children",
                  Age >= 15 ~ "adults",
                  is.na(Age) ~ "Unknown" # 16 IDs have no AGEYR
                ),
                ageGroup3 = case_when(
                  Age < 2 ~ "<2",
                  Age >= 2 & Age < 65 ~ "2-64",
                  Age >= 65 ~ "65+",
                  is.na(Age) ~ "Unknown"
                ),
                ageGroup5 = case_when( # edit 5 age bands
                  Age < 5 ~ "<5",
                  Age >= 5 & Age < 19 ~ "5-18",
                  Age >= 19 & Age < 31 ~ "19-30",
                  Age >= 31 & Age < 65 ~ "31-64",
                  Age >= 65 ~ "65+",
                  is.na(Age) ~ "Unknown" # 16 IDs have no Age
                  # TRUE ~ "Unknown" 
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
                ageGroup7 = case_when(
                  Age < 2 ~ "<2",
                  Age >= 2 & Age < 5 ~ "2-4",
                  Age >= 5 & Age < 15 ~ "5-14",
                  Age >= 15 & Age < 31 ~ "15-30", # Edit the Age-band into 15-30 & 31-44
                  Age >= 31 & Age < 45 ~ "31-44", # Edit the Age-band into 15-30 & 31-44
                  Age >= 45 & Age < 65 ~ "45-64",
                  Age >= 65 ~ "65+",
                  is.na(Age) ~ "Unknown" # 16 IDs have no AGEYR
                ),
                year = as.numeric(year)) %>% 
  glimpse()

write.csv(pop_l, "raw_data/nomis_population_long.csv", row.names = FALSE)

# Vaccination programme:
# SOURCE: https://www.gov.uk/government/publications/pneumococcal-the-green-book-chapter-25
vaccine_UK <- data.frame(
  year = c(2006, 2011),
  vaccine = c("PCV7", "PCV13")
)

col_map <- c(
  # Vaccination
  "Pre-PCV7" = "lightblue",
  "PCV7" = "gray70",
  "PCV13" = "gray20",
  # 5 age bands 
  "<5" = "indianred4",
  "5-18" = "orange",
  "19-30" = "seagreen4",
  "31-64" = "steelblue",
  "65+" = "purple3",
  "Unknown" = "black",
  # 7 age bands
  "<2" = "indianred4", 
  "2-4" = "indianred2", 
  "5-14" = "orange",
  "15-30" = "seagreen1", # Edit the Age-band into 15-30 & 31-44
  "31-44" = "seagreen4", # Edit the Age-band into 15-30 & 31-44
  "45-64" = "steelblue",
  "65+" = "purple3",
  # 2 age bands
  "children" = "darkred",
  "adults" = "darkblue",
  # Regions (From north to south)
  "North West" = "indianred4",
  "North East" = "steelblue",
  "Yorkshire and The Humber" = "seagreen4",
  "East Midlands" = "purple3",
  "West Midlands" = "orange",
  "East of England" = "indianred2",
  "London" = "seagreen1",
  "South East" = "deepskyblue",
  "South West" = "gold1",
  # Cases vs sequenced
  "Serotype 1 Case" = "gray75",
  "Sequenced" = "deepskyblue3",
  "Meningitis" = "green",
  "30 Day Death" = "maroon"
)

png("pictures/counts_hist_ages.png", width = 17, height = 12, unit = "cm", res = 1200)
hist(dat$Age)
dev.off()

# Seasonality #################################################################
all_year <- dat %>% 
  dplyr::group_by(year) %>% 
  dplyr::summarise(counts = n()) %>% 
  dplyr::ungroup()

all_month <- dat %>% 
  dplyr::group_by(month) %>% 
  dplyr::summarise(counts = n()) %>% 
  dplyr::ungroup()

# Generate new dataframe for yearMonth because no data in some yearMonth
all_yearMonth <- dplyr::left_join(
  data.frame(yearMonth = seq(as.Date(paste0(min(dat$year), "-01-01")),
                             as.Date(paste0(max(dat$year), "-12-01")),
                             by = "1 month"))
  ,
  dat %>% 
    dplyr::group_by(yearMonth) %>% 
    dplyr::summarise(counts = n()) %>% 
    dplyr::ungroup(),
  by = "yearMonth"
)

# Data viz!
viz_year <- ggplot(all_year, aes(x = year, y = counts)) +
  geom_line(size = 1.5) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  theme_bw()

viz_month <- ggplot(all_month, aes(x = month, y = counts)) +
  geom_line(size = 1.5) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  theme_bw()

# Additional data to viz seasonality pattern in yearMonth
december_lines <- seq(from = as.Date(paste0(min(dat$year)-1, "-12-01")),
                      to = as.Date(paste0(max(dat$year)-1, "-12-01")),
                      by = "1 year")


blocks <- data.frame(
  xmin = as.Date(paste0((min(dat$year) - 1):(max(dat$year) - 1), "-12-01")),  # Start December
  xmax = as.Date(paste0((min(dat$year)):(max(dat$year)), "-03-01")),  # End March
  year = (min(dat$year) - 1):(max(dat$year) - 1)
)

viz_yearMonth <- ggplot(all_yearMonth, aes(x = yearMonth, y = counts)) +
  geom_line(size = 1.5) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year", limits = c(as.Date("2015-12-01"), as.Date("2022-04-01"))) +
  geom_vline(xintercept = december_lines, linetype = "dashed", colour = "darkred") +
  geom_rect(data = blocks, aes(xmin = xmin, xmax = xmax,
                               ymin = -Inf, ymax = Inf),
            fill = "darkred",
            inherit.aes = FALSE, alpha = 0.3
  ) +
  ggtitle("The Counts of Serotype 12F by Months Specific to Every Year") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(size = 0.2, color = "grey80"),
    panel.grid.minor = element_line(size = 0.2, color = "grey80"),
    axis.line = element_line(color = "black")
  )

p1 <- cowplot::plot_grid(viz_year, viz_month, labels = c("A", "B"))

# Annoying plot for label
cow_yearMonth <- cowplot::plot_grid(NULL, viz_yearMonth,
                                    ncol = 1,
                                    rel_heights = c(0.1, 2),
                                    labels = c("", "C"))

png("pictures/counts_seasonality_allAges.png", width = 20, height = 18, unit = "cm", res = 1200)
p2 <- cowplot::plot_grid(p1, cow_yearMonth,
                         ncol = 1,
                         rel_heights = c(0.75, 1))
p2
dev.off()


# Age Group ####################################################################
ageGroup2_y <- dat %>% 
  dplyr::group_by(year, ageGroup2) %>% 
  dplyr::summarise(counts = n()) %>% 
  dplyr::ungroup()

ageGroup2_yM <- dat %>% 
  dplyr::group_by(yearMonth, ageGroup2) %>% 
  dplyr::summarise(counts = n()) %>% 
  dplyr::ungroup()

ageGroup5_y <- dat %>% 
  dplyr::group_by(year, ageGroup5) %>% 
  dplyr::summarise(counts = n()) %>% 
  dplyr::ungroup()

ageGroup5_yM <- dat %>% 
  dplyr::group_by(yearMonth, ageGroup5) %>% 
  dplyr::summarise(counts = n()) %>% 
  dplyr::ungroup()

ageGroup7_y <- dat %>% 
  dplyr::group_by(year, ageGroup7) %>% 
  dplyr::summarise(counts = n()) %>% 
  dplyr::ungroup()

ageGroup7_yM <- dat %>% 
  dplyr::group_by(yearMonth, ageGroup7) %>% 
  dplyr::summarise(counts = n()) %>% 
  dplyr::ungroup()

# Data viz!
viz_agegroup2_y <- ggplot(ageGroup2_y, aes(x = year, y = counts, group = ageGroup2,
                                           color = ageGroup2)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = c(col_map),
                     name = "Demographic",
                     breaks = c("children", "adults"),
                     labels = c("Children (< 15)", "Adults")
  ) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  theme_bw() +
  theme(legend.position = "none")

viz_agegroup2_yM <- ggplot(ageGroup2_yM, aes(x = yearMonth, y = counts, group = ageGroup2,
                                             color = ageGroup2)) +
  geom_line(size = 1) +
  scale_color_manual(values = c(col_map),
                     name = "Demographic",
                     breaks = c("children", "adults"),
                     labels = c("Children (< 15)", "Adults")
  ) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year", limits = c(as.Date("2015-12-01"), as.Date("2022-04-01"))) +
  theme_bw() +
  theme(
    legend.text = element_text(size=10),
    legend.title = element_text(size=11),
    legend.key.size = unit(0.25, "cm")
  )

viz_agegroup5_y <- ggplot(ageGroup5_y, aes(x = year, y = counts, group = ageGroup5,
                                           color = ageGroup5)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = c(col_map),
                     name = "Demographic",
                     breaks = c("<5", "5-18", "19-30", "31-64", "65+", "Unknown"),
                     labels = c("<5", "5-18", "19-30", "31-64", "65+", "Unknown")
  ) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  theme_bw() +
  theme(legend.position = "none")

viz_agegroup5_yM <- ggplot(ageGroup5_yM, aes(x = yearMonth, y = counts, group = ageGroup5,
                                             color = ageGroup5)) +
  geom_line(size = 1) +
  scale_color_manual(values = c(col_map),
                     name = "Demographic",
                     breaks = c("<5", "5-18", "19-30", "31-64", "65+", "Unknown"),
                     labels = c("<5", "5-18", "19-30", "31-64", "65+", "Unknown")
  ) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year", limits = c(as.Date("2015-12-01"), as.Date("2022-04-01"))) +
  theme_bw() +
  theme(
    legend.text = element_text(size=10),
    legend.title = element_text(size=11),
    legend.key.size = unit(0.25, "cm")
  )

viz_agegroup7_y <- ggplot(ageGroup7_y, aes(x = year, y = counts, group = ageGroup7,
                                           color = ageGroup7)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = c(col_map),
                     name = "Demographic",
                     breaks = c("<2", "2-4", "5-14", "15-30", "31-44", "45-64", "65+", "Unknown"),
                     labels = c("<2", "2-4", "5-14", "15-30", "31-44", "45-64", "65+", "Unknown")
  ) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  theme_bw() +
  theme(legend.position = "none")

viz_agegroup7_yM <- ggplot(ageGroup7_yM, aes(x = yearMonth, y = counts, group = ageGroup7,
                                             color = ageGroup7)) +
  geom_line(size = 1) +
  scale_color_manual(values = c(col_map),
                     name = "Demographic",
                     breaks = c("<2", "2-4", "5-14", "15-30", "31-44", "45-64", "65+", "Unknown"),
                     labels = c("<2", "2-4", "5-14", "15-30", "31-44", "45-64", "65+", "Unknown")
  ) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year", limits = c(as.Date("2015-12-01"), as.Date("2022-04-01"))) +
  theme_bw() +
  theme(
    legend.text = element_text(size=10),
    legend.title = element_text(size=11),
    legend.key.size = unit(0.25, "cm")
  )

p3_y <- cowplot::plot_grid(viz_agegroup2_y, viz_agegroup5_y, viz_agegroup7_y,
                           ncol = 1,
                           labels = c("A", "B", "C"))

p3_yM <- cowplot::plot_grid(viz_agegroup2_yM, viz_agegroup5_yM, viz_agegroup7_yM,
                            ncol = 1)

png("pictures/counts_ageHeterogeneity_allAges.png", width = 20, height = 18, unit = "cm", res = 1200)
p3 <- cowplot::plot_grid(p3_y, p3_yM,
                         ncol = 2,
                         rel_widths = c(0.65, 1.25))
p3
dev.off()
