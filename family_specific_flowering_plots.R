

# COMPLETE FLOWERING SPECIES-SPECIFIC AND FAMILY-LEVEL ANALYSIS
# LUKI AND YANGAMBI FORESTS
#
# Based on six existing phyr::pglmm models
#
# 



# 1. LOAD PACKAGES


packages <- c(
  "phyr",
  "dplyr",
  "tidyr",
  "purrr",
  "stringr",
  "ggplot2",
  "Matrix"
)

missing_packages <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

invisible(
  lapply(packages, library, character.only = TRUE)
)



# 2. CHECK REQUIRED OBJECTS


required_objects <- c(
  "luki_onset_no_year2",
  "luki_end_no_year2",
  "luki_duration_no_year2",
  "yangambi_onset_no_year2",
  "yangambi_end_no_year2",
  "yangambi_duration_no_year2",
  "luki_raw",
  "sp.list_yangambi"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    envir = .GlobalEnv
  )
]

if (length(missing_objects) > 0) {
  stop(
    "Missing objects: ",
    paste(missing_objects, collapse = ", ")
  )
}



# 3. EXTRACT SPECIES-SPECIFIC TEMPORAL SLOPES


# This function extracts:
# - Fixed temporal slope
# - Non-phylogenetic species random-slope contribution
# - Phylogenetic species random-slope contribution
# - Total species-specific slope
#
# Random-effect components:
# 3 = 0 + year_c | species
# 4 = 0 + year_c | species__
#
# The extraction uses the internal PGLMM matrices.

extract_species_slopes <- function(model, forest_name, measure_name) {
  
  # Check model components
  required_components <- c(
    "B", "St", "Zt", "s2r", "iV", "H", "random.effects"
  )
  
  missing_components <- setdiff(
    required_components,
    names(model)
  )
  
  if (length(missing_components) > 0) {
    stop(
      forest_name, " ", measure_name,
      ": missing model components: ",
      paste(missing_components, collapse = ", ")
    )
  }
  
  # Fixed year coefficient
  fixed_slope <- as.numeric(
    model$B["year_c", 1]
  )
  
  # Identify the two random-slope variance components
  variance_names <- names(model$s2r)
  
  identity_index <- grep(
    "^0 \\+ year_c\\|species$",
    variance_names
  )
  
  phylogenetic_index <- grep(
    "^0 \\+ year_c\\|species__$",
    variance_names
  )
  
  if (length(identity_index) != 1 ||
      length(phylogenetic_index) != 1) {
    stop(
      "Could not uniquely identify the two random-slope ",
      "variance components in ", forest_name, " ", measure_name,
      ". Inspect names(model$s2r)."
    )
  }
  
  identity_variance <- as.numeric(
    model$s2r[identity_index]
  )
  
  phylogenetic_variance <- as.numeric(
    model$s2r[phylogenetic_index]
  )
  
  # Determine random-effect block sizes
  block_sizes <- apply(
    as.matrix(model$St != 0),
    1,
    sum
  )
  
  block_start <- c(
    1,
    cumsum(block_sizes)[-length(block_sizes)] + 1
  )
  
  block_end <- cumsum(block_sizes)
  
  if (length(block_sizes) < 4) {
    stop(
      "Expected four random-effect blocks in ",
      forest_name, " ", measure_name
    )
  }
  
  # Extract species identity random-slope design matrix
  Zt_identity <- model$Zt[
    block_start[3]:block_end[3],
    ,
    drop = FALSE
  ]
  
  # Extract phylogenetic random-slope design matrix
  Zt_phylogenetic <- model$Zt[
    block_start[4]:block_end[4],
    ,
    drop = FALSE
  ]
  
  # Species names in the model's species random-effect order
  species_names <- levels(
    model$random.effects[[3]]$species
  )
  
  if (is.null(species_names)) {
    stop(
      "Could not retrieve species names from the third ",
      "random-effect component in ",
      forest_name, " ", measure_name
    )
  }
  
  # Extract non-phylogenetic species random-slope deviations
  identity_slope <- identity_variance *
    as.numeric(
      Zt_identity %*%
        model$iV %*%
        model$H
    )
  
  # Retrieve phylogenetic covariance matrix
  phy_covariance <- model$random.effects[[4]]$covar
  
  # Cholesky transformation of phylogenetic covariance
  phy_chol <- t(
    chol(as.matrix(phy_covariance))
  )
  
  # Extract phylogenetic species random-slope deviations
  phylogenetic_slope <- phylogenetic_variance *
    as.numeric(
      phy_chol %*%
        Zt_phylogenetic %*%
        model$iV %*%
        model$H
    )
  
  # Check dimensions
  if (length(identity_slope) != length(species_names)) {
    stop(
      "Identity slope length does not match species names in ",
      forest_name, " ", measure_name
    )
  }
  
  if (length(phylogenetic_slope) != length(species_names)) {
    stop(
      "Phylogenetic slope length does not match species names in ",
      forest_name, " ", measure_name
    )
  }
  
  # Total species-specific temporal slope
  total_slope <- fixed_slope +
    identity_slope +
    phylogenetic_slope
  
  # Directional interpretation
  direction <- ifelse(
    total_slope < 0,
    ifelse(
      measure_name == "Onset",
      "Earlier onset",
      ifelse(
        measure_name == "End",
        "Earlier termination",
        "Shorter duration"
      )
    ),
    ifelse(
      total_slope > 0,
      ifelse(
        measure_name == "Onset",
        "Later onset",
        ifelse(
          measure_name == "End",
          "Later termination",
          "Longer duration"
        )
      ),
      "No change"
    )
  )
  
  # Return species-specific results
  tibble(
    forest = forest_name,
    measure = measure_name,
    species = species_names,
    fixed_slope = fixed_slope,
    identity_slope = identity_slope,
    phylogenetic_slope = phylogenetic_slope,
    total_slope = total_slope,
    direction = direction
  )
}



# 4. APPLY EXTRACTION TO THE SIX EXISTING MODELS


species_luki_onset <- extract_species_slopes(
  luki_onset_no_year2, "Luki", "Onset"
)

species_luki_end <- extract_species_slopes(
  luki_end_no_year2, "Luki", "End"
)

species_luki_duration <- extract_species_slopes(
  luki_duration_no_year2, "Luki", "Duration"
)


species_yangambi_onset <- extract_species_slopes(
  yangambi_onset_no_year2, "Yangambi", "Onset"
)

species_yangambi_end <- extract_species_slopes(
  yangambi_end_no_year2, "Yangambi", "End"
)

species_yangambi_duration <- extract_species_slopes(
  yangambi_duration_no_year2, "Yangambi", "Duration"
)


# Combine the six extraction results
all_species_slopes <- bind_rows(
  species_luki_onset,
  species_luki_end,
  species_luki_duration,
  species_yangambi_onset,
  species_yangambi_end,
  species_yangambi_duration
)



# 5. CHECK SPECIES COUNTS AND SLOPE CALCULATIONS


print(
  all_species_slopes %>%
    group_by(forest, measure) %>%
    summarise(
      species = n_distinct(species),
      records = n(),
      .groups = "drop"
    ),
  n = Inf
)

# Check total slope arithmetic
slope_check <- all_species_slopes %>%
  mutate(
    calculated_total =
      fixed_slope +
      identity_slope +
      phylogenetic_slope,
    difference = total_slope - calculated_total
  ) %>%
  summarise(
    maximum_difference = max(
      abs(difference),
      na.rm = TRUE
    )
  )

print(slope_check)



# 6. PREPARE SPECIES-TO-FAMILY LOOKUP


# Luki taxonomy
luki_family_lookup <- luki_raw %>%
  transmute(
    forest = "Luki",
    species = as.character(species),
    family = as.character(family)
  )

# Corrected Yangambi taxonomy
yangambi_family_lookup <- sp.list_yangambi %>%
  transmute(
    forest = "Yangambi",
    species = as.character(species),
    family = as.character(family)
  )

# Combine family lookup tables
family_lookup <- bind_rows(
  luki_family_lookup,
  yangambi_family_lookup
) %>%
  mutate(
    species = str_replace_all(
      str_squish(species),
      " ",
      "_"
    ),
    family = str_to_sentence(
      str_to_lower(str_squish(family))
    )
  ) %>%
  filter(
    !is.na(species),
    species != "",
    !is.na(family),
    family != ""
  ) %>%
  distinct(forest, species, .keep_all = TRUE)


# Check duplicate species-family assignments
duplicate_family_assignments <- family_lookup %>%
  count(forest, species) %>%
  filter(n > 1)

if (nrow(duplicate_family_assignments) > 0) {
  print(duplicate_family_assignments)
  stop("Resolve duplicate species-family assignments.")
}



# 7. ATTACH BOTANICAL FAMILIES TO SPECIES SLOPES


species_slopes <- all_species_slopes %>%
  mutate(
    species = str_replace_all(
      str_squish(species),
      " ",
      "_"
    )
  ) %>%
  left_join(
    family_lookup,
    by = c("forest", "species")
  )

# Check missing family assignments
missing_families <- species_slopes %>%
  filter(is.na(family)) %>%
  distinct(forest, species)

if (nrow(missing_families) > 0) {
  print(missing_families)
  warning(
    "Some species have no family assignment. ",
    "Check the species names and taxonomy lookup."
  )
}

# Keep species with family assignments for family summaries
species_slopes_family <- species_slopes %>%
  filter(
    !is.na(family),
    family != ""
  )

# Check final species-level dataset
cat(
  "Total species-measure records:",
  nrow(species_slopes), "\n"
)

cat(
  "Unique species in Luki:",
  n_distinct(
    species_slopes %>%
      filter(forest == "Luki") %>%
      pull(species)
  ), "\n"
)

cat(
  "Unique species in Yangambi:",
  n_distinct(
    species_slopes %>%
      filter(forest == "Yangambi") %>%
      pull(species)
  ), "\n"
)



# 8. OVERALL SPECIES-LEVEL SLOPE STATISTICS


species_slope_statistics <- species_slopes %>%
  group_by(forest, measure) %>%
  summarise(
    n_species = n_distinct(species),
    
    mean_slope = mean(total_slope, na.rm = TRUE),
    median_slope = median(total_slope, na.rm = TRUE),
    sd_slope = sd(total_slope, na.rm = TRUE),
    
    minimum_slope = min(total_slope, na.rm = TRUE),
    maximum_slope = max(total_slope, na.rm = TRUE),
    
    negative_n = sum(total_slope < 0, na.rm = TRUE),
    positive_n = sum(total_slope > 0, na.rm = TRUE),
    
    .groups = "drop"
  )

print(species_slope_statistics, n = Inf)



# 9. SPECIES DIRECTIONAL CONSISTENCY


species_direction_summary <- species_slopes %>%
  group_by(forest, measure) %>%
  summarise(
    total_species = n_distinct(species),
    
    negative_n = sum(total_slope < 0, na.rm = TRUE),
    positive_n = sum(total_slope > 0, na.rm = TRUE),
    
    negative_percent =
      100 * negative_n / total_species,
    
    positive_percent =
      100 * positive_n / total_species,
    
    .groups = "drop"
  )

print(species_direction_summary, n = Inf)



# 10. FAMILY-LEVEL SLOPE SUMMARY


family_slope_summary <- species_slopes_family %>%
  group_by(forest, measure, family) %>%
  summarise(
    n_species = n_distinct(species),
    
    mean_slope = mean(total_slope, na.rm = TRUE),
    median_slope = median(total_slope, na.rm = TRUE),
    
    Q1 = unname(
      quantile(total_slope, 0.25, na.rm = TRUE)
    ),
    
    Q3 = unname(
      quantile(total_slope, 0.75, na.rm = TRUE)
    ),
    
    IQR = IQR(total_slope, na.rm = TRUE),
    SD = sd(total_slope, na.rm = TRUE),
    
    negative_n = sum(total_slope < 0, na.rm = TRUE),
    positive_n = sum(total_slope > 0, na.rm = TRUE),
    
    negative_percent =
      100 * negative_n / n_species,
    
    positive_percent =
      100 * positive_n / n_species,
    
    .groups = "drop"
  )

print(family_slope_summary, n = Inf)



# 11. FAMILY RESPONSE CONSISTENCY


# 80% is a descriptive threshold, not a significance test.

family_response_consistency <- family_slope_summary %>%
  mutate(
    dominant_direction = case_when(
      negative_percent >= 80 ~ "Predominantly negative",
      positive_percent >= 80 ~ "Predominantly positive",
      TRUE ~ "Mixed response"
    )
  ) %>%
  select(
    forest,
    measure,
    family,
    n_species,
    negative_n,
    positive_n,
    negative_percent,
    positive_percent,
    dominant_direction
  )

print(family_response_consistency, n = Inf)



# 12. FAMILY COMPARISON BETWEEN FORESTS


family_comparison <- family_slope_summary %>%
  select(
    forest,
    measure,
    family,
    n_species,
    mean_slope,
    median_slope,
    Q1,
    Q3,
    negative_percent,
    positive_percent
  ) %>%
  pivot_wider(
    names_from = forest,
    values_from = c(
      n_species,
      mean_slope,
      median_slope,
      Q1,
      Q3,
      negative_percent,
      positive_percent
    )
  ) %>%
  arrange(measure, family)

print(family_comparison, n = Inf)



# 13. FAMILY-SPECIFIC PLOT FUNCTION


plot_family_forest <- function(
    data,
    forest_name,
    measure_name
) {
  
  plot_data <- data %>%
    filter(
      forest == forest_name,
      measure == measure_name
    ) %>%
    filter(
      is.finite(mean_slope),
      is.finite(Q1),
      is.finite(Q3)
    ) %>%
    arrange(mean_slope)
  
  if (nrow(plot_data) == 0) {
    stop(
      "No data for ",
      forest_name, " - ", measure_name
    )
  }
  
  if (measure_name == "Onset") {
    negative_label <- "Earlier flowering onset"
    positive_label <- "Later flowering onset"
    
  } else if (measure_name == "End") {
    negative_label <- "Earlier flowering termination"
    positive_label <- "Later flowering termination"
    
  } else if (measure_name == "Duration") {
    negative_label <- "Shorter flowering duration"
    positive_label <- "Longer flowering duration"
    
  } else {
    stop("Unknown measure.")
  }
  
  plot_data <- plot_data %>%
    mutate(
      phenological_shift = case_when(
        mean_slope < 0 ~ negative_label,
        mean_slope > 0 ~ positive_label,
        TRUE ~ "Zero slope"
      ),
      Plant_Family = factor(
        family,
        levels = family
      )
    )
  
  n_families <- nrow(plot_data)
  n_species <- sum(plot_data$n_species)
  
  x_range <- range(
    c(
      plot_data$Q1,
      plot_data$Q3,
      plot_data$mean_slope
    ),
    na.rm = TRUE
  )
  
  padding <- diff(x_range) * 0.12
  
  if (!is.finite(padding) || padding == 0) {
    padding <- 0.5
  }
  
  x_min <- x_range[1] - padding
  x_max <- x_range[2] + padding
  
  ggplot(
    plot_data,
    aes(
      x = mean_slope,
      y = Plant_Family,
      colour = phenological_shift
    )
  ) +
    
    geom_vline(
      xintercept = 0,
      colour = "grey45",
      linetype = "dashed",
      linewidth = 0.7
    ) +
    
    # Q1-Q3: interquartile range
    geom_segment(
      aes(
        x = Q1,
        xend = Q3,
        yend = Plant_Family
      ),
      linewidth = 0.8
    ) +
    
    # Family mean
    geom_point(size = 2.8) +
    
    scale_colour_manual(
      name = "Phenological shift",
      values = setNames(
        c("#2C7FB8", "#E66101", "grey50"),
        c(
          negative_label,
          positive_label,
          "Zero slope"
        )
      ),
      breaks = c(
        negative_label,
        positive_label
      )
    ) +
    
    scale_x_continuous(
      limits = c(x_min, x_max),
      breaks = pretty(
        c(x_min, x_max),
        n = 6
      ),
      expand = c(0, 0)
    ) +
    
    labs(
      x = expression(
        "Mean species-specific temporal slope (days year"^{-1}*")"
      ),
      y = "Plant family",
      
      title = paste(
        forest_name,
        "- Flowering",
        tolower(measure_name),
        "temporal trends across families"
      ),
      
      subtitle = paste0(
        "Point = mean species slope; ",
        "line = interquartile range (Q1-Q3)\n",
        n_species, " species across ",
        n_families, " families"
      )
    ) +
    
    theme_bw(base_size = 11) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 12,
        hjust = 0.5
      ),
      
      plot.subtitle = element_text(
        size = 9,
        hjust = 0.5
      ),
      
      axis.text.y = element_text(
        colour = "black",
        size = ifelse(n_families > 50, 6.5, 8.5),
        face = "italic"
      ),
      
      axis.text.x = element_text(
        colour = "black"
      ),
      
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      
      legend.position = "right",
      
      plot.margin = margin(12, 15, 12, 15)
    )
}



# 14. CREATE SIX FAMILY PLOTS


family_luki_onset_plot <- plot_family_forest(
  family_slope_summary, "Luki", "Onset"
)

family_luki_end_plot <- plot_family_forest(
  family_slope_summary, "Luki", "End"
)

family_luki_duration_plot <- plot_family_forest(
  family_slope_summary, "Luki", "Duration"
)

family_yangambi_onset_plot <- plot_family_forest(
  family_slope_summary, "Yangambi", "Onset"
)

family_yangambi_end_plot <- plot_family_forest(
  family_slope_summary, "Yangambi", "End"
)

family_yangambi_duration_plot <- plot_family_forest(
  family_slope_summary, "Yangambi", "Duration"
)


# Display plots
print(family_luki_onset_plot)
print(family_luki_end_plot)
print(family_luki_duration_plot)

print(family_yangambi_onset_plot)
print(family_yangambi_end_plot)
print(family_yangambi_duration_plot)



# 15. EXPORT RESULTS


output_folder <- file.path(
  getwd(),
  "Flowering_Species_and_Family_Results"
)

dir.create(
  output_folder,
  showWarnings = FALSE,
  recursive = TRUE
)


# Save species-specific slopes
write.csv(
  species_slopes,
  file.path(
    output_folder,
    "Species_specific_flowering_slopes.csv"
  ),
  row.names = FALSE
)

# Save overall species statistics
write.csv(
  species_slope_statistics,
  file.path(
    output_folder,
    "Overall_species_slope_statistics.csv"
  ),
  row.names = FALSE
)

# Save species directional consistency
write.csv(
  species_direction_summary,
  file.path(
    output_folder,
    "Species_directional_consistency.csv"
  ),
  row.names = FALSE
)

# Save family slope summary
write.csv(
  family_slope_summary,
  file.path(
    output_folder,
    "Family_slope_summary.csv"
  ),
  row.names = FALSE
)

# Save family response consistency
write.csv(
  family_response_consistency,
  file.path(
    output_folder,
    "Family_response_consistency.csv"
  ),
  row.names = FALSE
)

# Save forest comparison
write.csv(
  family_comparison,
  file.path(
    output_folder,
    "Family_comparison_Luki_Yangambi.csv"
  ),
  row.names = FALSE
)


# Save all six plots
plot_list <- list(
  Luki_Onset = family_luki_onset_plot,
  Luki_End = family_luki_end_plot,
  Luki_Duration = family_luki_duration_plot,
  Yangambi_Onset = family_yangambi_onset_plot,
  Yangambi_End = family_yangambi_end_plot,
  Yangambi_Duration = family_yangambi_duration_plot
)

for (nm in names(plot_list)) {
  
  forest_name <- ifelse(
    grepl("^Luki_", nm),
    "Luki",
    "Yangambi"
  )
  
  measure_name <- sub(
    "^(Luki|Yangambi)_",
    "",
    nm
  )
  
  n_families <- family_slope_summary %>%
    filter(
      forest == forest_name,
      measure == measure_name
    ) %>%
    nrow()
  
  plot_height <- max(
    7,
    min(24, 0.20 * n_families + 2.5)
  )
  
  ggsave(
    file.path(
      output_folder,
      paste0("Figure_Family_", nm, ".png")
    ),
    plot = plot_list[[nm]],
    width = 9,
    height = plot_height,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    file.path(
      output_folder,
      paste0("Figure_Family_", nm, ".pdf")
    ),
    plot = plot_list[[nm]],
    width = 9,
    height = plot_height,
    units = "in",
    bg = "white"
  )
}


# Save the combined species-level extraction too
write.csv(
  all_species_slopes,
  file.path(
    output_folder,
    "All_species_slope_components.csv"
  ),
  row.names = FALSE
)


message(
  "Analysis completed. Results saved in: ",
  normalizePath(output_folder)
)



# Check species numbers for every forest and measure
species_slopes %>%
  group_by(forest, measure) %>%
  summarise(
    species = n_distinct(species),
    records = n(),
    .groups = "drop"
  )

# Check family assignment
species_slopes %>%
  group_by(forest) %>%
  summarise(
    total_species = n_distinct(species),
    species_with_family = n_distinct(
      species[!is.na(family)]
    ),
    missing_family = sum(is.na(family)),
    .groups = "drop"
  )

# Check the number of families
family_slope_summary %>%
  group_by(forest) %>%
  summarise(
    number_of_families = n_distinct(family),
    .groups = "drop"
  )

# Check that all three measures are represented
table(
  species_slopes$forest,
  species_slopes$measure
)

# Check species slope arithmetic
species_slopes %>%
  mutate(
    recalculated_slope =
      fixed_slope + identity_slope + phylogenetic_slope,
    difference = total_slope - recalculated_slope
  ) %>%
  summarise(
    maximum_difference = max(abs(difference), na.rm = TRUE)
  )