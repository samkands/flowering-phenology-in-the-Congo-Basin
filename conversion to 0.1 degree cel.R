# Standard deviation used in the original standardization
temp_sd_luki <- sd(luki_raw$mean_temp, na.rm = TRUE)

# Flowering onset
beta_temp_onset <- 21.4188
onset_effect_01 <- beta_temp_onset * (0.1 / temp_sd_luki)

# Flowering duration
beta_temp_duration <- -17.0588
duration_effect_01 <- beta_temp_duration * (0.1 / temp_sd_luki)

# Display results
cat("Temperature SD:", round(temp_sd_luki, 6), "°C\n\n")

cat("Flowering onset effect:",
    round(onset_effect_01, 2), "days\n")

cat("Flowering duration effect:",
    round(duration_effect_01, 2), "days\n")