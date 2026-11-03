

# SPECIES-SPECIFIC FLOWERING TEMPORAL TRENDS: DENSITY PLOTS
# Luki and Yangambi forests


# 1. LOAD PACKAGES

library(ggplot2)
library(dplyr)


# 2. PREPARE SPECIES-SPECIFIC SLOPE DATA

species_slopes_plot <- species_slopes %>%
  mutate(
    
    Forest = factor(
      forest,
      levels = c("Luki", "Yangambi")
    ),
    
    Variable = factor(
      measure,
      levels = c("Onset", "End", "Duration"),
      labels = c(
        "Flowering Onset",
        "Flowering End Date",
        "Flowering Duration"
      )
    ),
    
    Slope = total_slope
    
  ) %>%
  filter(
    is.finite(Slope),
    !is.na(Forest),
    !is.na(Variable)
  )


# 3. CHECK NUMBER OF SPECIES

cat("\n--- NUMBER OF SPECIES IN DENSITY PLOTS ---\n")

print(
  species_slopes_plot %>%
    group_by(Forest, Variable) %>%
    summarise(
      n_species = n_distinct(species),
      .groups = "drop"
    )
)


# 4. CREATE DENSITY PLOT

p_density <- ggplot(
  species_slopes_plot,
  aes(
    x = Slope,
    fill = Forest,
    colour = Forest
  )
) +
  
  # Density curves
  
  geom_density(
    alpha = 0.38,
    linewidth = 1,
    adjust = 1,
    na.rm = TRUE
  ) +
  
  # Zero reference line
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "black",
    linewidth = 0.8
  ) +
  
  # Separate panels for phenological metrics
  
  facet_wrap(
    ~ Variable,
    nrow = 1,
    scales = "free_y"
  ) +
  
  # Forest colours
  
  scale_fill_manual(
    values = c(
      "Luki" = "#D94B45",
      "Yangambi" = "#58758D"
    )
  ) +
  
  scale_colour_manual(
    values = c(
      "Luki" = "#D94B45",
      "Yangambi" = "#58758D"
    )
  ) +
  
  # Axis labels and title
  
  labs(
    
    title = "Species-specific temporal trends in flowering phenology",
    
    x = expression(
      "Species-specific phenological trend (days year"^{-1}*")"
    ),
    
    y = "Density",
    
    fill = "Forest",
    
    colour = "Forest"
    
  ) +
  
  # Theme
  
  theme_classic(
    base_size = 15
  ) +
  
  theme(
    
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    
    axis.title.x = element_text(
      size = 16,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 16,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 13,
      colour = "black"
    ),
    
    strip.text = element_text(
      size = 15,
      face = "bold"
    ),
    
    legend.position = "bottom",
    
    legend.title = element_text(
      face = "bold",
      size = 13
    ),
    
    legend.text = element_text(
      size = 12
    ),
    
    panel.spacing = grid::unit(
      1.2,
      "lines"
    )
    
  )


# 5. DISPLAY PLOT

print(p_density)


# ==============================================================================
# 6. SAVE HIGH-RESOLUTION PNG
# ==============================================================================

ggsave(
  
  filename = "Species_specific_flowering_trends_density.png",
  
  plot = p_density,
  
  width = 14,
  height = 6,
  
  units = "in",
  
  dpi = 600,
  
  bg = "white"
  
)


# ==============================================================================
# 7. SAVE PDF
# ==============================================================================

ggsave(
  
  filename = "Species_specific_flowering_trends_density.pdf",
  
  plot = p_density,
  
  width = 14,
  height = 6,
  
  units = "in",
  
  bg = "white"
  
)


# ==============================================================================
# 8. EXPORT DATA USED FOR THE FIGURE
# ==============================================================================

write.csv(
  
  species_slopes_plot,
  
  file = "Species_specific_flowering_trends_density_data.csv",
  
  row.names = FALSE
  
)


cat("\nDensity plots successfully saved.\n")
