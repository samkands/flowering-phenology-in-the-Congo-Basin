

# INDEPENDENT FOREST PLOTS FOR LUKI AND YANGAMBI PHENOLOGY



# 1. SET WORKING DIRECTORY AND LOAD PACKAGES


setwd("E:/PAPER 2")

library(ggplot2)
library(dplyr)



# 2. CREATE MASTER DATAFRAME WITH MODEL ESTIMATES


climate_effects_master <- data.frame(
  
  Site = c(
    "Luki", "Luki", "Luki",
    "Luki", "Luki", "Luki",
    "Luki", "Luki",
    "Yangambi", "Yangambi",
    "Yangambi", "Yangambi", "Yangambi",
    "Yangambi", "Yangambi"
  ),
  
  Metric = c(
    "Onset", "Onset", "Onset",
    "End", "End", "End",
    "Duration", "Duration",
    "Onset", "Onset",
    "End", "End", "End",
    "Duration", "Duration"
  ),
  
  Term = c(
    "Temperature", "Precipitation", "Temp x Precip",
    "Temperature", "Precipitation", "Temp x Precip",
    "Temperature", "Precipitation",
    "Temperature", "Precipitation",
    "Temperature", "Precipitation", "Temp x Precip",
    "Temperature", "Precipitation"
  ),
  
  Estimate = c(
    # Luki Onset
    21.4188, 4.6270, 30.3277,
    
    # Luki End
    3.0707, 5.6206, 21.3541,
    
    # Luki Duration
    -17.0588, 1.1005,
    
    # Yangambi Onset
    -6.3817, 1.4552,
    
    # Yangambi End
    0.13951, -0.12480, 0.68032,
    
    # Yangambi Duration
    6.5289, -2.0809
  ),
  
  StdError = c(
    # Luki Onset
    9.2826, 8.9034, 17.7961,
    
    # Luki End
    6.1094, 5.8689, 11.7156,
    
    # Luki Duration
    4.9075, 4.8468,
    
    # Yangambi Onset
    6.3386, 6.9058,
    
    # Yangambi End
    4.75712, 4.21850, 4.43438,
    
    # Yangambi Duration
    4.0963, 4.4248
  ),
  
  PValue = c(
    # Luki Onset
    0.02103, 0.60328, 0.08835,
    
    # Luki End
    0.61524, 0.33822, 0.06835,
    
    # Luki Duration
    0.0005088, 0.8203771,
    
    # Yangambi Onset
    0.3140, 0.8331,
    
    # Yangambi End
    0.9766, 0.9764, 0.8781,
    
    # Yangambi Duration
    0.1110, 0.6382
  ),
  
  stringsAsFactors = FALSE
)



# 3. CALCULATE CONFIDENCE INTERVALS AND SIGNIFICANCE CATEGORIES


climate_effects_master <- climate_effects_master %>%
  
  mutate(
    
    # Calculate 95% confidence intervals
    Lower = Estimate - 1.96 * StdError,
    Upper = Estimate + 1.96 * StdError,
    
    # Three significance categories
    Significance = case_when(
      
      PValue < 0.05 ~ "Significant",
      
      PValue < 0.10 ~ "Marginally significant",
      
      TRUE ~ "Non-significant"
    ),
    
    Significance = factor(
      Significance,
      levels = c(
        "Significant",
        "Marginally significant",
        "Non-significant"
      )
    ),
    
    # Add significance symbols
    SigLabel = case_when(
      
      PValue < 0.001 ~ "***",
      
      PValue < 0.01 ~ "**",
      
      PValue < 0.05 ~ "*",
      
      PValue < 0.10 ~ "\u2020",
      
      TRUE ~ ""
    ),
    
    # Format coefficient labels
    LabelText = sprintf(
      "\u03B2 = %.2f%s",
      Estimate,
      SigLabel
    ),
    
    # Order phenological metrics
    Metric = factor(
      Metric,
      levels = c("Onset", "End", "Duration"),
      labels = c(
        "Flowering Onset",
        "Flowering End Date",
        "Flowering Duration"
      )
    ),
    
    # Order climate predictors
    Term = factor(
      Term,
      levels = c(
        "Temp x Precip",
        "Precipitation",
        "Temperature"
      ),
      labels = c(
        "Temperature x Precipitation",
        "Precipitation (Z)",
        "Temperature (Z)"
      )
    )
  )


# Inspect the statistical classifications
print(
  climate_effects_master %>%
    select(Site, Metric, Term, Estimate, StdError,
           Lower, Upper, PValue, Significance, LabelText)
)



# 4. DEFINE COMMON COLOUR SCALE


significance_colours <- c(
  
  "Significant" = "firebrick3",
  
  "Marginally significant" = "darkorange2",
  
  "Non-significant" = "black"
)



# 5. LUKI FOREST PLOT


luki_data <- climate_effects_master %>%
  filter(Site == "Luki")


p_luki <- ggplot(
  luki_data,
  aes(
    x = Estimate,
    y = Term,
    color = Significance
  )
) +
  
  # Reference line at zero
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.7
  ) +
  
  # 95% confidence intervals
  geom_errorbar(
    aes(
      xmin = Lower,
      xmax = Upper
    ),
    width = 0.15,
    linewidth = 0.8
  ) +
  
  # Effect estimates
  geom_point(
    size = 3.5,
    shape = 16
  ) +
  
  # Coefficient labels
  geom_text(
    aes(label = LabelText),
    vjust = -1.3,
    color = "black",
    size = 3,
    fontface = "bold"
  ) +
  
  # Separate phenological metrics
  facet_wrap(
    ~ Metric,
    ncol = 3,
    scales = "free_x"
  ) +
  
  # Significance colours
  scale_color_manual(
    values = significance_colours,
    drop = FALSE
  ) +
  
  # Axis labels and title
  labs(
    title = "Climate Drivers of Flowering Phenology: Luki Forest Site",
    
    x = "Effect Size (Days per SD Change in Climate Predictor)",
    
    y = "Climate Predictor",
    
    color = "Statistical Significance"
  ) +
  
  # Plot theme
  theme_bw(base_size = 12) +
  
  theme(
    
    plot.title = element_text(
      size = 13,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    
    strip.background = element_rect(
      fill = "grey92",
      color = "black",
      linewidth = 0.6
    ),
    
    strip.text = element_text(
      size = 10,
      face = "bold"
    ),
    
    axis.title.x = element_text(
      face = "bold",
      margin = margin(t = 10)
    ),
    
    axis.title.y = element_text(
      face = "bold",
      margin = margin(r = 10)
    ),
    
    panel.grid.major.y = element_line(
      color = "grey85",
      linetype = "solid"
    ),
    
    panel.grid.minor = element_blank(),
    
    legend.position = "bottom"
  )


# Display Luki plot
print(p_luki)


# Save Luki plot
ggsave(
  filename = "Luki_Phenology_Forest_Plot.tiff",
  plot = p_luki,
  device = "tiff",
  width = 10,
  height = 4.5,
  dpi = 300,
  compression = "lzw"
)



# 6. YANGAMBI FOREST PLOT


yangambi_data <- climate_effects_master %>%
  filter(Site == "Yangambi")


p_yangambi <- ggplot(
  yangambi_data,
  aes(
    x = Estimate,
    y = Term,
    color = Significance
  )
) +
  
  # Reference line at zero
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.7
  ) +
  
  # 95% confidence intervals
  geom_errorbar(
    aes(
      xmin = Lower,
      xmax = Upper
    ),
    width = 0.15,
    linewidth = 0.8
  ) +
  
  # Effect estimates
  geom_point(
    size = 3.5,
    shape = 17
  ) +
  
  # Coefficient labels
  geom_text(
    aes(label = LabelText),
    vjust = -1.3,
    color = "black",
    size = 3,
    fontface = "bold"
  ) +
  
  # Separate phenological metrics
  facet_wrap(
    ~ Metric,
    ncol = 3,
    scales = "free_x"
  ) +
  
  # Significance colours
  scale_color_manual(
    values = significance_colours,
    drop = FALSE
  ) +
  
  # Axis labels and title
  labs(
    title = "Climate Drivers of Flowering Phenology: Yangambi Forest Site",
    
    x = "Effect Size (Days per SD Change in Climate Predictor)",
    
    y = "Climate Predictor",
    
    color = "Statistical Significance"
  ) +
  
  # Plot theme
  theme_bw(base_size = 12) +
  
  theme(
    
    plot.title = element_text(
      size = 13,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    
    strip.background = element_rect(
      fill = "grey92",
      color = "black",
      linewidth = 0.6
    ),
    
    strip.text = element_text(
      size = 10,
      face = "bold"
    ),
    
    axis.title.x = element_text(
      face = "bold",
      margin = margin(t = 10)
    ),
    
    axis.title.y = element_text(
      face = "bold",
      margin = margin(r = 10)
    ),
    
    panel.grid.major.y = element_line(
      color = "grey85",
      linetype = "solid"
    ),
    
    panel.grid.minor = element_blank(),
    
    legend.position = "bottom"
  )


# Display Yangambi plot
print(p_yangambi)


# Save Yangambi plot
ggsave(
  filename = "Yangambi_Phenology_Forest_Plot.tiff",
  plot = p_yangambi,
  device = "tiff",
  width = 10,
  height = 4.5,
  dpi = 300,
  compression = "lzw"
)



# 7. EXPORT STATISTICAL RESULTS


write.csv(
  climate_effects_master,
  file = "Climate_Effects_Significance_Summary.csv",
  row.names = FALSE
)

