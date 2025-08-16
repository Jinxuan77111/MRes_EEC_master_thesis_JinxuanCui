## Master thesis 
setwd()

# 1. calculatation: habitat loss and habitat gain 
# Habitat loss: forest cover & distance to forest edge
# forest cover

# Load packages
library(terra)
library(sf)
library(landscapemetrics)
library(vegan)
library(ggplot2)
library(dplyr)

# Load data 
# Sampling sites from three farms
sites <- read.csv("sampling_points.csv")
sites <- st_as_sf(sites, coords = c("Longitude", "Latitude"), crs = 4326)
plot(sites[c('Site')], key.pos=4, axes=TRUE)

# Load 2023 Braizl land coverage and use data (updated version)
landcover <- rast("brasil_coverage_2023.tif")
print(landcover)
res(landcover)
crs(landcover)

# Crop the landcover raster to the interested region - easy for calculation
# Extract buffer area around sampling region 
sites_region <- st_as_sfc(st_bbox(sites))
sites_region <- st_buffer(sites_region, 1000)

# Crop the landcover map to interested region
sites_landcover <- crop(landcover, sites_region)

# Convert landcover map to binary forest map (1 = forest, 0 = non-forest) 
# Based on the legend from MapBiomas, value = c(1, 3, 4, 5, 6, 49) considered as forest
# In interested region, only value = 3 exists.
sites_forest <- sites_landcover == 3
print(sites_forest)
plot(sites_forest)


# Convert to the UTM 23S (unit:meters) 
sites_utm23S <- st_transform(sites, 32723)
sites_forest_utm23S <- project(sites_forest, "epsg:32723", res=30, method='near')
plot(st_geometry(sites_utm23S), add=TRUE, col = "red")


# Calculation of forest cover with varying buffer size
# radius = 60m
lsm_60 <- sample_lsm(sites_forest_utm23S, sites_utm23S, 
                     shape = "circle", size = 60,  
                     plot_id = sites_utm23S$Site, 
                     what = c('lsm_c_pland'))
# radius = 100m
lsm_100 <- sample_lsm(sites_forest_utm23S, sites_utm23S, 
                      shape = "circle", size = 100,  
                      plot_id = sites_utm23S$Site, 
                      what = c('lsm_c_pland'))
# radius = 200m
lsm_200 <- sample_lsm(sites_forest_utm23S, sites_utm23S, 
                      shape = "circle", size = 200,  
                      plot_id = sites_utm23S$Site, 
                      what = c('lsm_c_pland'))
# radius = 300m
lsm_300 <- sample_lsm(sites_forest_utm23S, sites_utm23S, 
                      shape = "circle", size = 300,  
                      plot_id = sites_utm23S$Site, 
                      what = c('lsm_c_pland'))
# radius = 400m
lsm_400 <- sample_lsm(sites_forest_utm23S, sites_utm23S, 
                      shape = "circle", size = 400,  
                      plot_id = sites_utm23S$Site, 
                      what = c('lsm_c_pland'))
# radius = 500m
lsm_500 <- sample_lsm(sites_forest_utm23S, sites_utm23S, 
                      shape = "circle", size = 500,  
                      plot_id = sites_utm23S$Site, 
                      what = c('lsm_c_pland'))


# Save the result to csv file
lsm_60_df <- as.data.frame(lsm_50)
write.csv(lsm_60_df, "lsm_60_forest_cover.csv", row.names = FALSE)


# Calculation of distance to forest edge
# Convert the binary forest raster to polygons
sites_landcover_vect <- sites_forest_utm23S
sites_landcover_vect <- as.polygons(sites_landcover_vect)

# convert to the simple feature and define the forest edge
forest_sf <- st_as_sf(sites_landcover_vect)
forest_edge <- st_boundary(forest_sf)
plot(forest_edge)

# Check the location of sampling points: interior, exterior the forest 
sites_utm23S$inside_forest <- st_within(sites_utm23S, forest_sf)

# measure the distance 
sites_utm23S$distance_to_edge <- st_distance(sites_utm23S, forest_edge)

# Save the data frame to a CSV file
sites_utm23S_df <- as.data.frame(sites_utm23S)
write.csv(sites_utm23S_df, "distance_to_edge.csv", row.names = FALSE)



# 2.Invesitigate the scale of effect for forest cover
# Load the packages 
library(dplyr)
library(lme4)
library(ggplot2)
library(lubridate)
library(sjPlot)
library(MuMIn)

# Load the data 
landscape <- read.csv("Landscape_result.csv")
baseline <- read.csv("Baseline_output.csv")
baseline <- left_join(baseline, landscape, by = c("Site", "Year"))

# scale the predictor
baseline$Forest_cover_60m_scaled <- scale(baseline$Forest_cover_60m)
baseline$Forest_cover_100m_scaled <- scale(baseline$Forest_cover_100m)
baseline$Forest_cover_200m_scaled <- scale(baseline$Forest_cover_200m)
baseline$Forest_cover_300m_scaled <- scale(baseline$Forest_cover_300m)
baseline$Forest_cover_400m_scaled <- scale(baseline$Forest_cover_400m)
baseline$Forest_cover_500m_scaled <- scale(baseline$Forest_cover_500m)
baseline$Distance_to_edge_scaled <- scale(baseline$Distance_to_edge.m.)


## Fit the models for varying buffer size
# 60m
fc60m <- glmer(Validation ~ Forest_cover_60m_scaled 
               + Distance_to_edge_scaled + (1 | Site), 
               data = baseline, 
               family = binomial)

# 100m
fc100m <- glmer(Validation ~ Forest_cover_100m_scaled 
                + Distance_to_edge_scaled + (1 | Site), 
                data = baseline, 
                family = binomial)

# 200m
fc200m <- glmer(Validation ~ Forest_cover_200m_scaled 
                + Distance_to_edge_scaled + (1 | Site), 
                data = baseline, 
                family = binomial)

# 300m
fc300m <- glmer(Validation ~ Forest_cover_300m_scaled 
                + Distance_to_edge_scaled + (1 | Site), 
                data = baseline, 
                family = binomial)

# 400m
fc400m <- glmer(Validation ~ Forest_cover_400m_scaled 
                + Distance_to_edge_scaled + (1 | Site), 
                data = baseline, 
                family = binomial)

# 500m
fc500m <- glmer(Validation ~ Forest_cover_500m_scaled 
                + Distance_to_edge_scaled + (1 | Site), 
                data = baseline, 
                family = binomial)

# Fit the model without forest cover for calculation of R^2
fc_reduced <- glmer(Validation ~ Distance_to_edge_scaled + (1 | Site), 
                       data = baseline, 
                       family = binomial)


#### 3. Testing for Hypothesis
# Load the packages 
library(dplyr)
library(lme4)
library(ggplot2)
library(lubridate)
library(sjPlot)
library(MuMIn)


# Hypothesis 1: habitat loss vs probaility of frog presence
# Load the data 
landscape <- read.csv("Landscape_result.csv")
baseline <- read.csv("Baseline_output.csv")
baseline <- left_join(baseline, landscape, by = c("Site", "Year"))

# scale the predictors
baseline$Forest_cover_100m_scaled <- scale(baseline$Forest_cover_100m)
baseline$Distance_to_edge_scaled <- scale(baseline$Distance_to_edge.m.)

# check coreleation for predictors
cor(baseline$Forest_cover_100m_scaled, baseline$Distance_to_edge_scaled)

# fit the glmm model
model1 <- glmer(Validation ~ Forest_cover_100m_scaled 
                + Distance_to_edge_scaled + (1 | Site), 
                data = baseline, 
                family = binomial)

summary(model1)

# better model selection
model2 <- glmer(Validation ~ Forest_cover_100m_scaled + (1 | Site), 
                data = baseline, 
                family = binomial)

summary(model2)


# Hypothesis 2: habitat gain vs probaility of frog presence
## Fit lm model 
prop_difference <- read.csv("prop_data_time_23_24.csv")
lm_model <- lm(Difference ~ time_after_reforestation_days, data = prop_difference)
summary(lm_model)

# Check diagnoistic plots
plot(lm_model$fitted.values, lm_model$residuals,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)
qqnorm(lm_model$residuals)
qqline(lm_model$residuals, col = "red")
plot(lm_model$fitted.values, sqrt(abs(scale(lm_model$residuals))),
     xlab = "Fitted values",
     ylab = "Square root of standardized residuals",
     main = "Scale-Location Plot")



### 4.Figure plotting
# Load the packages 
library(dplyr)
library(lme4)
library(ggplot2)
library(lubridate)
library(sjPlot)
library(MuMIn)
library(tidyr)
library(magick)
library(cowplot)

# Hypothesis 1
# Calculation of mean and sd for forest cover 
mean_fc <- mean(baseline$Forest_cover_100m, na.rm = TRUE)
sd_fc <- sd(baseline$Forest_cover_100m, na.rm = TRUE)
baseline <- baseline %>%
  mutate(
    Forest_cover_100m_unscaled = Forest_cover_100m_scaled * sd_fc + mean_fc
  )

# Summarize frog presence proportion per site fro better visualise 
prop_data <- baseline %>%
  group_by(Site) %>%
  summarise(
    forest_cover = mean(Forest_cover_100m_unscaled, na.rm = TRUE),
    frog_presence_rate = mean(Validation, na.rm = TRUE)
  )


# Create new data frame for prediction based on GLMM
new_data <- data.frame(
  Forest_cover_100m_scaled = seq(min(baseline$Forest_cover_100m_scaled, na.rm = TRUE),
                                 max(baseline$Forest_cover_100m_scaled, na.rm = TRUE),
                                 length.out = 100)
)
# Remove random effect in plot
new_data$Site <- NA  

# Prediction based on GLMM and back to percentage of forest cover instread of scaled value
pred_link <- predict(model4, newdata = new_data, re.form = NA, type = "link", se.fit = TRUE)

new_data <- new_data %>%
  mutate(
    fit_link = pred_link$fit,
    se_link = pred_link$se.fit,
    lower_link = fit_link - 1.96 * se_link,
    upper_link = fit_link + 1.96 * se_link,
    pred = 1 / (1 + exp(-fit_link)),          
    lower_CI = 1 / (1 + exp(-lower_link)),    
    upper_CI = 1 / (1 + exp(-upper_link)),    
    Forest_cover_100m_unscaled = Forest_cover_100m_scaled * sd_fc + mean_fc
  )

# Plot the figure using ggplot2
ggplot() +
  geom_ribbon(data = new_data, aes(x = Forest_cover_100m_unscaled, 
                                   ymin = lower_CI, ymax = upper_CI),
              fill = "#3B6E52", alpha = 0.15, color = "gray80", size = 0.1) +
  geom_line(data = new_data, aes(x = Forest_cover_100m_unscaled, y = pred), 
            color = "#3B6E52", size = 1.2) +
  geom_point(data = prop_data, aes(x = forest_cover, y = frog_presence_rate), 
             size = 2.7, 
             color = "#3B6E52",    
             shape = 17,
             alpha = 0.9) +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  labs(
    x = "Forest Cover (%)",
    y = "Proportion of Frog Presence",
    title = "",
    subtitle = ""
  ) +
  annotate("text", x = 93, y = 0.36, label = "", 
           size = 5, fontface = "italic", hjust = 1) +
  theme_classic(base_family = "Times") +
  theme(
    axis.title.x = element_text(size = 16, margin = margin(t = 11)),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    axis.text.x = element_text(size = 13),   
    axis.text.y = element_text(size = 13),
    panel.grid.minor = element_blank()
  )


# Hypothsis 2
# Create new data frame for prediction
new_data_difference <- data.frame(
  time_after_reforestation_days = seq(
    min(prop_difference$time_after_reforestation_days, na.rm = TRUE),
    max(prop_difference$time_after_reforestation_days, na.rm = TRUE),
    length.out = 100
  )
)

prediction_difference <- predict(lm_model, newdata = new_data_difference, se.fit = TRUE)

new_data_difference <- new_data_difference %>%
  mutate(
    fit = prediction_difference$fit,
    se = prediction_difference$se.fit,
    lower_CI = fit - 1.96 * se,
    upper_CI = fit + 1.96 * se
  )

# Plot the figure using ggplot2
ggplot() +
  geom_ribbon(data = new_data_difference,
              aes(x = time_after_reforestation_days, ymin = lower_CI, ymax = upper_CI),
              fill = "#e7c66b", alpha = 0.2, color = "gray80", size = 0.1) +
  geom_line(data = new_data_difference,
            aes(x = time_after_reforestation_days, y = fit),
            color = "#f3a361", size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray30", size = 0.6) +
  geom_point(data = prop_difference,
             aes(x = time_after_reforestation_days, y = Difference),
             shape = 17,
             size = 3,
             color = "#3B6E52",
             alpha = 0.9) +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.25), limits = c(-1, 1)) +
  labs(
    x = "Time Since Reforestation (days)",
    y = "Difference in Proportion of Frog Presence",
    title = ""
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, margin = margin(t = 11)),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    panel.grid.minor = element_blank()
  ) +
  annotate("text", x = max(new_data_difference$time_after_reforestation_days), y = 0.55,
           label = "", size = 4, fontface = "italic", hjust = 1)














