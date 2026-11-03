


# SECTION 1: ENVIRONMENT SETUP & LIBRARY INITIALIZATION

setwd("E:/PAPER 2/PGLMM")

# Required package checks and installation
required_packages <- c("dplyr", "glmmTMB", "ggplot2", "patchwork", "VennDiagram", 
                       "V.PhyloMaker", "ape", "phytools", "raster", "stringr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(dplyr)
library(stringr)
library(raster)
library(glmmTMB)
library(ggplot2)
library(patchwork)
library(VennDiagram)
library(V.PhyloMaker)
library(ape)
library(phytools)

# Conflict Resolution Rule
# Force R to use dplyr versions for overlapping commands

select <- dplyr::select
rename <- dplyr::rename

# Load raw phenology data

full_data <- readRDS("PhenoDRC.rds")


# SECTION 2: RAW DATA FILTERING & TAXONOMIC CLEANING

flowering <- full_data %>%
  filter(
    phenophase == "flowers",     # Isolate the flowers phenophase
    value == 1,                  # Keep only active presence records
    !is.na(species),             # Remove rows missing a species name
    species != "sp.",            # Remove unresolved species labels
    species != "NA"
  ) %>%
  mutate(
    # Clean taxonomic authorities or noise (e.g., "puberula H." -> "puberula")
    species = str_extract(species, "^[a-z0-9_]+"),
    # Clean historical binomial format if column exists
    historical_genus_species = str_extract(historical_genus_species, "^[A-Za-z0-9_]+\\s+[a-z0-9_]+")
  ) %>%
  # Filter out any rows where the regex extraction failed to match
  filter(!is.na(species))



# Clean strings and force tree-safe text connections (convert dashes to underscores)

library(dplyr)
library(stringr)
library(tidyr)

flowering <- flowering %>%
  select(
    -any_of(c(
      "latitude", "longitude", "elevation", "Id", "id",
      "specie_qa", "species_qa", "historical_common_name",
      "species"
    ))
  ) %>%
  relocate(
    historical_genus_species, family, genus,
    .after = coord_type
  ) %>%
  filter(!is.na(historical_genus_species)) %>%
  mutate(
    historical_genus_species =
      str_squish(historical_genus_species)
  ) %>%
  filter(
    !str_detect(
      historical_genus_species,
      regex("\\b(sp|spp|cf|cfr)\\b", ignore_case = TRUE)
    ),
    str_detect(historical_genus_species, "\\s+"),
    !is.na(site),
    !is.na(year),
    !is.na(week),
    !is.na(doy),
    !is.na(genus),
    !is.na(family)
  )

# --- VERIFICATION CHECKS ---
print("Unique sites remaining in data:")
print(unique(flowering$site))

print("Total number of raw rows isolated:")
print(nrow(flowering))

## Rename historical_genus_species to species and put an _

library(dplyr)
library(stringr)

flowering <- flowering %>%
  rename(
    species = historical_genus_species
  ) %>%
  mutate(
    species = str_replace_all(
      str_trim(species),
      "\\s+",
      "_"
    )
  )

# Check the updated column names
names(flowering)

# Check the species names
flowering %>%
  select(species) %>%
  distinct() %>%
  print(n = 20)

#Save the main dataset

write.csv(
  flowering,
  file = "flowering.csv",
  row.names = FALSE
)


# Extract clean spatial operational subsets

luki_flowering <- flowering %>%
  filter(site == "Luki")


print(paste("Luki isolated observations:", nrow(luki_flowering)))


luki_flowering <- luki_flowering %>%
  mutate(
    across(c(family, genus, species), str_trim),
    species = str_replace_all(species, "-", "_"),
    species_binomial = paste(genus, species, sep = "_")
  )


# Post-cleaning diagnostic verification scan for non-alphanumeric text remains
luki_flowering <- luki_flowering %>% filter(str_detect(species, "[^A-Za-z0-9_]"))

print(paste("Lingering anomalies in Luki:", nrow(luki_anomalies)))


# Climate data:  NETCDF CLIMATE EXTRACTION & PROJECTION MATCHING

library(raster)
library(dplyr)

setwd("E:/CLIMATE DATA")

# Load CRU TS monthly climate grids
tmp_nc <- raster::brick(
  "cru_ts4.09.1901.2024.tmp.dat.nc",
  varname = "tmp"
)

pre_nc <- raster::brick(
  "cru_ts4.09.1901.2024.pre.dat.nc",
  varname = "pre"
)

# Luki coordinates


sites <- data.frame(
  site      = "Luki",
  latitude  = -5.62,  
  longitude = 13.17
)

# Study period
study_years <- 1947:1957

# Verify the number of monthly layers
nlayers(tmp_nc)
nlayers(pre_nc)

# CRU TS should have 12 monthly layers per year
stopifnot(nlayers(tmp_nc) == nlayers(pre_nc))

# Build the time index for every layer
full_time <- expand.grid(
  month = 1:12,
  year = 1901:2024
)

# expand.grid varies year fastest; reorder chronologically
full_time <- full_time[order(full_time$year, full_time$month), ]
rownames(full_time) <- NULL

# Check that the time index matches the raster layers
stopifnot(nrow(full_time) == nlayers(tmp_nc))

# Identify the layers for the study years
index <- which(full_time$year %in% study_years)

stopifnot(length(index) == length(study_years) * 12)

time_subset <- full_time[index, ]

# Extract all monthly values for Luki
tmp_all <- raster::extract(
  tmp_nc,
  sites[, c("longitude", "latitude")]
)

pre_all <- raster::extract(
  pre_nc,
  sites[, c("longitude", "latitude")]
)

# Ensure the extraction results are matrices
tmp_all <- as.matrix(tmp_all)
pre_all <- as.matrix(pre_all)

# Select the correct months
tmp_vals <- tmp_all[, index, drop = FALSE]
pre_vals <- pre_all[, index, drop = FALSE]

# Calculate annual climate summaries

climate_summary <- data.frame(
  site = "Luki",
  year = study_years,
  mean_temp = vapply(study_years, function(y) {
    mean(tmp_vals[1, time_subset$year == y], na.rm = TRUE)
  }, numeric(1)),
  total_precip = vapply(study_years, function(y) {
    sum(pre_vals[1, time_subset$year == y], na.rm = TRUE)
  }, numeric(1))
)

print(climate_summary)

# 8. MERGE DATA


library(dplyr)

luki_raw <- luki_flowering %>%
  group_by(species,  genus, family, year) %>%
  filter(sum(!is.na(doy)) >= 5) %>%
  summarise(
    n_obs = sum(!is.na(doy)),
    onset_doy = quantile(doy, 0.10, na.rm = TRUE, type = 7),
    end_doy = quantile(doy, 0.90, na.rm = TRUE, type = 7),
    .groups = "drop"
  ) %>%
  mutate(
    site = "Luki",
    duration_doy = end_doy - onset_doy
  ) %>%
  left_join(climate_summary, by = c("site", "year")) %>%
  mutate(
    year_c = year - 1947,
    precip_100 = total_precip / 100
  ) %>%
  select(
    site, family, genus, species,year,
    n_obs, onset_doy, end_doy, duration_doy,
    mean_temp, total_precip, year_c, precip_100
  )

library(dplyr)

species_year_coverage <- luki_raw %>%
  group_by(species, genus, family) %>%
  summarise(
    n_years = n_distinct(year),
    first_year = min(year),
    last_year = max(year),
    .groups = "drop"
  ) %>%
  arrange(n_years)

summary(species_year_coverage$n_years)

table(species_year_coverage$n_years)


#SECTION 7: INDEPENDENT TEMPORAL BOUNDS & DISTRIBUTIONS


# Apply specific historical subset ranges per site
luki_raw <- luki_raw %>%
  filter(year >= 1947 & year <= 1957)

# --- FINAL QUALITY SANITY CONTROLS ---
print(paste("Total Luki analysis rows:", nrow(luki_raw)))

print("Confirming DOY fields are structured as continuous numeric variables:")
print(is.numeric(luki_raw$onset_doy))
print(paste("Unique value diversity:", length(unique(luki_raw$onset_doy))))

# Generate visual distribution layout curve
plot(density(luki_raw$onset_doy, na.rm = TRUE), 
     main = "Density of Onset Days (Luki)", xlab = "Day of Year")


# reconstruct the phylogeny using phylomaker

#Put an _ between species


sp.list <- luki_raw %>%
  dplyr::select(species, genus, family) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    family = dplyr::case_when(
      genus == "Erythrophleum" ~ "Fabaceae",
      genus == "Hexalobus"     ~ "Annonaceae",
      genus == "Mammea"        ~ "Calophyllaceae",
      genus == "Musanga"       ~ "Urticaceae",
      genus == "Schrebera"     ~ "Oleaceae",
      genus == "Staudtia"      ~ "Myristicaceae",
      TRUE ~ as.character(family)
    )
  ) %>%
  as.data.frame()


# Replace spaces with underscores in species names

sp.list$species <- gsub(" ", "_", sp.list$species)

rownames(sp.list) <- NULL


nodes.info.1 <- nodes.info.1 %>%
  as.data.frame()


#  LOAD V.PHYLOMAKER DATA & UPDATE TAXONOMY

data(GBOTB.extended)
data(nodes.info.1)

phylo_result <- V.PhyloMaker::phylo.maker(
  sp.list = sp.list,
  tree = GBOTB.extended,
  nodes = nodes.info.1,
  scenarios = "S3"
)


# Extract the phylogenetic tree

phy <- phylo_result$scenario.3


# Save the phylogenetic tree

ape::write.tree(
  phy,
  file = "phylogenetic_tree_130_species.newick"
)

# Save the cleaned species list
write.csv(
  sp.list,
  file = "cleaned_species_list_130.csv",
  row.names = FALSE
)

# Plot tree
plot(
  phy,
  cex = 0.5,
  no.margin = TRUE
)

# Match clean dataset to the new tree

luki_raw$species_tree <- gsub(" ", "_", luki_raw$species)

df_phy <- luki_raw %>%
  filter(species_tree %in% phy$tip.label)


# CREATE PHYLOGENETIC COVARIANCE MATRIX 

phy_cov <- vcv(phy, corr = TRUE)



# PHYLOGENETIC SIGNAL TEST 

library(dplyr)
library(ape)
library(phytools)


# 1. Calculate mean flowering onset for each species

trait_df <- luki_raw %>%
  group_by(species) %>%
  summarise(
    trait_value = mean(end_doy, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(trait_value))

# 2. Create a named trait vector

trait <- trait_df$trait_value
names(trait) <- trait_df$species

# 3. Identify species shared by the trait data and tree
common_sp <- intersect(names(trait), phy$tip.label)

cat("Species in trait data:", length(trait), "\n")
cat("Species in phylogeny:", length(phy$tip.label), "\n")
cat("Matching species:", length(common_sp), "\n")

# 4. Prune the tree to matching species
phy_sub <- ape::keep.tip(phy, common_sp)

# 5. Align trait values to the tree's tip order
trait_aligned <- trait[phy_sub$tip.label]

# 6. Verify alignment
stopifnot(
  identical(names(trait_aligned), phy_sub$tip.label),
  !anyNA(trait_aligned)
)

cat("Aligned trait values:", length(trait_aligned), "\n")
cat("Tree tips:", length(phy_sub$tip.label), "\n")

# 7. Calculate Blomberg's K with a randomization test

K_result <- phytools::phylosig(
  tree = phy_sub,
  x = trait_aligned,
  method = "K",
  test = TRUE,
  nsim = 999
)

print(K_result)

# trait         K   P_value
# 1    Flowering onset 0.177034  0.34034
# 2      Flowering end 0.125825  0.692693
# 3 Flowering duration 0.109788  0.826827 

##PART 2. Climate sensitivity 

##Standardization 

luki_raw <- luki_raw %>%
  mutate(
    temp_z   = as.numeric(scale(mean_temp)),
    precip_z = as.numeric(scale(total_precip))
  )


install.packages("phyr")
library(phyr)

#Luki onset with interaction temp_z * precip_z

# Remove Scottellia_klaineana

luki_raw <- luki_raw %>%
  filter(species != "Scottellia_klaineana") %>%
  mutate(species = droplevels(factor(species)))

# Confirm 130 species remain
n_distinct(luki_raw$species)

# Fit the phylogenetic mixed model

#Climate sensitivity

#Climate Onset

luki_onset <- pglmm(
  onset_doy ~ temp_z * precip_z +
    (1 | species__) +
    (1 | year),
  data = luki_raw,
  family = "gaussian",
  cov_ranef = list(species = phy)
)

summary(luki_onset)

# Linear mixed model fit by restricted maximum likelihood
# 
# Call:onset_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
#  -6063  12142  12173 
# 
# Random effects:
#              Variance Std.Dev
# 1|species   2.421e+03   49.20
# 1|species__ 3.957e+02   19.89
# 1|year      8.409e+07 9170.00
# residual    6.010e+03   77.52
# 
# Fixed effects:
#                     Value Std.Error Zscore Pvalue
# (Intercept)      101.9303 3221.9828 0.0316 0.9748
# temp_z            21.4116 2959.3190 0.0072 0.9942
# precip_z           4.7439 2831.3200 0.0017 0.9987
# temp_z:precip_z   30.2427 5671.3816 0.0053 0.9957



luki_onset_additive <- pglmm(
  onset_doy ~ temp_z + precip_z +
    (1 | species__) +
    (1 | year),
  data = luki_raw,
  family = "gaussian",
  cov_ranef = list(species = phy)
)

summary(luki_onset_additive)

# Linear mixed model fit by restricted maximum likelihood
# 
# Call:onset_doy ~ temp_z + precip_z
# 
# logLik    AIC    BIC 
#  -6076  12165  12192 
# 
# Random effects:
#              Variance Std.Dev
# 1|species   2.378e+03   48.77
# 1|species__ 3.951e+02   19.88
# 1|year      8.646e+07 9298.48
# residual    6.013e+03   77.54
# 
# Fixed effects:
#                 Value Std.Error Zscore Pvalue
# (Intercept)  110.6396 2810.2962 0.0394 0.9686
# temp_z        17.9110 2926.1312 0.0061 0.9951
# precip_z       3.2212 2856.2534 0.0011 0.9991


luki_onset_with_bobyqa <- phyr::pglmm(
  onset_doy ~ temp_z * precip_z +
    (1 | species__) +
    (1 | year),
  data = luki_raw,
  family = "gaussian",
  cov_ranef = list(species = phy),
  REML = FALSE,
  optimizer = "bobyqa",
  maxit = 1000,
  reltol = 1e-8,
  verbose = FALSE
)

summary(luki_onset_year_bobyqa)

# Linear mixed model fit by maximum likelihood
# 
# Call:onset_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
#  -6055  12125  12156 
# 
# Random effects:
#             Variance Std.Dev
# 1|species     2471.1   49.71
# 1|species__    351.6   18.75
# 1|year         760.4   27.57
# residual      6055.6   77.82
# 
# Fixed effects:
#                    Value Std.Error Zscore    Pvalue    
# (Intercept)     102.4101   15.5492 6.5862 4.512e-11 ***
# temp_z           21.4188    9.2826 2.3074   0.02103 *  
# precip_z          4.6270    8.9034 0.5197   0.60328    
# temp_z:precip_z  30.3277   17.7961 1.7042   0.08835 .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


luki_onset_year_boby <- phyr::pglmm(
  onset_doy ~ temp_z + precip_z +
    (1 | species__) +
    (1 | year),
  data = luki_raw,
  family = "gaussian",
  cov_ranef = list(species = phy),
  REML = FALSE,
  optimizer = "bobyqa",
  maxit = 1000,
  reltol = 1e-8,
  verbose = FALSE
)

summary(luki_onset_year_boby)

# summary(luki_onset_year_boby)
# Linear mixed model fit by maximum likelihood
# 
# Call:onset_doy ~ temp_z + precip_z
# 
# logLik    AIC    BIC 
#  -6111  12236  12262 
# 
# Random effects:
#              Variance Std.Dev
# 1|species   2.413e+03   49.12
# 1|species__ 3.902e+02   19.75
# 1|year      6.338e+07 7960.87
# residual    5.988e+03   77.38
# 
# Fixed effects:
#                 Value Std.Error Zscore Pvalue
# (Intercept)  110.7215 2406.0365 0.0460 0.9633
# temp_z        17.9136 2505.2001 0.0072 0.9943
# precip_z       3.2169 2445.3744 0.0013 0.9990

#luki_onset_with_bobyqa is the BEST MODEL

##END PHENOLOGY

#Luki end with interaction temp_z * precip_z

luki_end <- pglmm(
  end_doy ~ temp_z * precip_z +
    (1 | species__) +
    (1 | year),
  data = luki_raw,
  family = "gaussian",
  cov_ranef = list(species = phy)
)


summary(luki_end)

# Linear mixed model fit by restricted maximum likelihood
# 
# Call:end_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
#  -5922  11860  11891 
# 
# Random effects:
#              Variance Std.Dev
# 1|species   2.985e+03   54.63
# 1|species__ 3.987e+02   19.97
# 1|year      5.177e+07 7195.12
# residual    4.349e+03   65.94
# 
# Fixed effects:
#                     Value Std.Error Zscore Pvalue
# (Intercept)      262.5060 2528.1002 0.1038 0.9173
# temp_z             2.9530 2321.9916 0.0013 0.9990
# precip_z           5.9419 2221.5589 0.0027 0.9979
# temp_z:precip_z   21.4493 4449.9767 0.0048 0.9962


# No Interaction Term

luki_end_additive <- pglmm(
  end_doy ~ temp_z + precip_z +
    (1 | species__) +
    (1 | year),
  data = luki_raw,
  family = "gaussian",
  cov_ranef = list(species = phy)
)

summary(luki_end_additive)
# Linear mixed model fit by restricted maximum likelihood
# 
# Call:end_doy ~ temp_z + precip_z
# 
# logLik    AIC    BIC 
#  -5893  11800  11826 
# 
# Random effects:
#             Variance Std.Dev
# 1|species     3175.1   56.35
# 1|species__    303.0   17.41
# 1|year         574.1   23.96
# residual      4377.6   66.16
# 
# Fixed effects:
#                 Value Std.Error  Zscore Pvalue    
# (Intercept) 269.87959  13.69764 19.7026 <2e-16 ***
# temp_z        0.55382   7.85216  0.0705 0.9438    
# precip_z      4.68612   7.68444  0.6098 0.5420    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Model for flowering End (end_doy) with bobyqa

luki_end_with_bobyqa <- phyr::pglmm(
  end_doy ~ temp_z * precip_z +
    (1 | species__) +
    (1 | year),
  data = luki_raw,
  family = "gaussian",
  cov_ranef = list(species = phy),
  REML = FALSE,
  optimizer = "bobyqa",
  maxit = 1000,
  reltol = 1e-8,
  verbose = FALSE
)

summary(luki_end_with_bobyqa)

# Linear mixed model fit by maximum likelihood
# 
# Call:end_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
#  -5911  11838  11869 
# 
# Random effects:
#             Variance Std.Dev
# 1|species     3034.7   55.09
# 1|species__    338.3   18.39
# 1|year         309.7   17.60
# residual      4385.9   66.23
# 
# Fixed effects:
#                    Value Std.Error  Zscore  Pvalue    
# (Intercept)     263.3678   13.5924 19.3761 < 2e-16 ***
# temp_z            3.0707    6.1094  0.5026 0.61524    
# precip_z          5.6206    5.8689  0.9577 0.33822    
# temp_z:precip_z  21.3541   11.7156  1.8227 0.06835 .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Luki additive with Bobyqa was the best model

luki_end_bobyqa <- phyr::pglmm(
  end_doy ~ temp_z + precip_z +
    (1 | species__) +
    (1 | year),
  data = luki_raw,
  family = "gaussian",
  cov_ranef = list(species = phy),
  REML = FALSE,
  optimizer = "bobyqa",
  maxit = 1000,
  reltol = 1e-8,
  verbose = FALSE
)

summary(luki_end_bobyqa)
# Linear mixed model fit by maximum likelihood
# 
# Call:end_doy ~ temp_z + precip_z
# 
# logLik    AIC    BIC 
#  -5955  11924  11951 
# 
# Random effects:
#              Variance Std.Dev
# 1|species      2952.5   54.34
# 1|species__     415.5   20.38
# 1|year      3222773.2 1795.21
# residual       4330.6   65.81
# 
# Fixed effects:
#                Value Std.Error Zscore Pvalue
# (Intercept) 268.5466  542.7191 0.4948 0.6207
# temp_z        0.4699  564.9368 0.0008 0.9993
# precip_z      4.8574  551.4460 0.0088 0.9930

#Luki_end_

##LUKI DURATION

luki_duration <- pglmm(
  duration_doy ~ temp_z * precip_z +
    (1 | species__) +
    (1 | year),
  data = luki_raw,
  family = "gaussian",
  cov_ranef = list(species = phy)
)

summary(luki_duration)

# summary(luki_duration)
# Linear mixed model fit by restricted maximum likelihood
# 
# Call:duration_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
#  -6238  12492  12522 
# 
# Random effects:
#              Variance   Std.Dev
# 1|species   3.100e+03 5.567e+01
# 1|species__ 3.710e-01 6.091e-01
# 1|year      1.452e+08 1.205e+04
# residual    8.748e+03 9.353e+01
# 
# Fixed effects:
#                      Value  Std.Error  Zscore Pvalue
# (Intercept)      161.91048 4233.58175  0.0382 0.9695
# temp_z           -18.43465 3888.47588 -0.0047 0.9962
# precip_z           0.99255 3720.28816  0.0003 0.9998
# temp_z:precip_z   -8.47562 7452.06262 -0.0011 0.9991

luki_duration_additive <- pglmm(
  duration_doy ~ temp_z + precip_z +
    (1 | species__) +
    (1 | year),
  data = luki_raw,
  family = "gaussian",
  cov_ranef = list(species = phy)
)

summary(luki_duration_additive)

# summary(luki_duration_additive)
# Linear mixed model fit by restricted maximum likelihood
# 
# Call:duration_doy ~ temp_z + precip_z
# 
# logLik    AIC    BIC 
#  -6250  12515  12541 
# 
# Random effects:
#              Variance   Std.Dev
# 1|species   3.102e+03 5.569e+01
# 1|species__ 3.712e-01 6.093e-01
# 1|year      1.453e+08 1.205e+04
# residual    8.738e+03 9.348e+01
# 
# Fixed effects:
#                 Value Std.Error  Zscore Pvalue
# (Intercept)  159.4479 3642.7874  0.0438 0.9651
# temp_z       -17.4543 3792.9688 -0.0046 0.9963
# precip_z       1.4216 3702.3904  0.0004 0.9997

luki_duration_year_bobyqa <- phyr::pglmm(
  duration_doy ~ temp_z * precip_z +
    (1 | species__) +
    (1 | year),
  data = luki_raw,
  family = "gaussian",
  cov_ranef = list(species = phy),
  REML = FALSE,
  optimizer = "bobyqa",
  maxit = 1000,
  reltol = 1e-8,
  verbose = FALSE
)

summary(luki_duration_year_bobyqa)

# summary(luki_duration_year_bobyqa)
# Linear mixed model fit by maximum likelihood
# 
# Call:duration_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
#  -6291  12597  12628 
# 
# Random effects:
#              Variance   Std.Dev
# 1|species   3.102e+03 5.570e+01
# 1|species__ 9.975e-01 9.988e-01
# 1|year      2.269e+08 1.506e+04
# residual    8.710e+03 9.333e+01
# 
# Fixed effects:
#                      Value  Std.Error  Zscore Pvalue
# (Intercept)      161.99783 5292.05225  0.0306 0.9756
# temp_z           -18.43476 4860.66502 -0.0038 0.9970
# precip_z           0.99553 4650.42730  0.0002 0.9998
# temp_z:precip_z   -8.47637 9315.21275 -0.0009 0.9993


#Duration additive with 

luki_duration_additive_bobyqa <- phyr::pglmm(
     duration_doy ~ temp_z + precip_z +
         (1 | species__) +
         (1 | year),
     data = luki_raw,
     family = "gaussian",
     cov_ranef = list(species = phy),
     REML = FALSE,
     optimizer = "bobyqa",
     maxit = 1000,
     reltol = 1e-8,
     verbose = FALSE
 )
 
# summary(luki_duration_additive_bobyqa)
# Linear mixed model fit by maximum likelihood
# 
# Call:duration_doy ~ temp_z + precip_z
# 
# logLik    AIC    BIC 
#  -6290  12595  12621 
# 
# Random effects:
#              Variance   Std.Dev
# 1|species   3.076e+03 5.546e+01
# 1|species__ 2.399e-02 1.549e-01
# 1|year      2.201e+08 1.484e+04
# residual    8.718e+03 9.337e+01
# 
# Fixed effects:
#                 Value Std.Error  Zscore Pvalue
# (Intercept)  159.4664 4484.2519  0.0356 0.9716
# temp_z       -17.4543 4669.1261 -0.0037 0.9970
# precip_z       1.4185 4557.6245  0.0003 0.9998

####YANGAMBI


# SECTION 1: ENVIRONMENT SETUP & LIBRARY INITIALIZATION

setwd("E:/PAPER 2/PGLMM")

# Required package checks and installation
required_packages <- c("dplyr", "glmmTMB", "ggplot2", "patchwork", "VennDiagram", 
                       "V.PhyloMaker", "ape", "phytools", "raster", "stringr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(dplyr)
library(stringr)
library(raster)
library(glmmTMB)
library(ggplot2)
library(patchwork)
library(VennDiagram)
library(V.PhyloMaker)
library(ape)
library(phytools)

# Conflict Resolution Rule
# Force R to use dplyr versions for overlapping commands

select <- dplyr::select
rename <- dplyr::rename

# Load raw phenology data

full_data <- readRDS("PhenoDRC.rds")

# SECTION 2: RAW DATA FILTERING & TAXONOMIC CLEANING

flowering <- full_data %>%
  filter(
    phenophase == "flowers",     # Isolate the flowers phenophase
    value == 1,                  # Keep only active presence records
    !is.na(species),             # Remove rows missing a species name
    species != "sp.",            # Remove unresolved species labels
    species != "NA"
  ) %>%
  mutate(
    # Clean taxonomic authorities or noise (e.g., "puberula H." -> "puberula")
    species = str_extract(species, "^[a-z0-9_]+"),
    # Clean historical binomial format if column exists
    historical_genus_species = str_extract(historical_genus_species, "^[A-Za-z0-9_]+\\s+[a-z0-9_]+")
  ) %>%
  # Filter out any rows where the regex extraction failed to match
  filter(!is.na(species))

# Clean strings and force tree-safe text connections (convert dashes to underscores)

library(dplyr)
library(stringr)
library(tidyr)

flowering <- flowering %>%
  select(
    -any_of(c(
      "latitude", "longitude", "elevation", "Id", "id",
      "specie_qa", "species_qa", "historical_common_name",
      "species"
    ))
  ) %>%
  relocate(
    historical_genus_species, family, genus,
    .after = coord_type
  ) %>%
  filter(!is.na(historical_genus_species)) %>%
  mutate(
    historical_genus_species =
      str_squish(historical_genus_species)
  ) %>%
  filter(
    !str_detect(
      historical_genus_species,
      regex("\\b(sp|spp|cf|cfr)\\b", ignore_case = TRUE)
    ),
    str_detect(historical_genus_species, "\\s+"),
    !is.na(site),
    !is.na(year),
    !is.na(week),
    !is.na(doy),
    !is.na(genus),
    !is.na(family)
  )

# --- VERIFICATION CHECKS ---

print("Unique sites remaining in data:")
print(unique(flowering$site))

print("Total number of raw rows isolated:")
print(nrow(flowering))

## Rename historical_genus_species to species and put an _

library(dplyr)
library(stringr)

flowering <- flowering %>%
  rename(
    species = historical_genus_species
  ) %>%
  mutate(
    species = str_replace_all(
      str_trim(species),
      "\\s+",
      "_"
    )
  )

# Check the updated column names
names(flowering)

# Check the species names
flowering %>%
  select(species) %>%
  distinct() %>%
  print(n = 20)

#Save the main dataset

write.csv(
  flowering,
  file = "flowering.csv",
  row.names = FALSE
)

# Extract clean spatial operational subsets

yangambi_flowering <- flowering %>%
  filter(site == "Yangambi")

print(paste("Yangambi isolated observations:", nrow(yangambi_flowering)))

yangambi_flowering <- yangambi_flowering %>%
  mutate(
    across(c(family, genus, species), str_trim),
    species = str_replace_all(species, "-", "_"),
    species_binomial = paste(genus, species, sep = "_")
  )

# Post-cleaning diagnostic verification scan for non-alphanumeric text remains
yangambi_flowering <- yangambi_flowering %>% filter(str_detect(species, "[^A-Za-z0-9_]"))

print(paste("Lingering anomalies in Yangambi:", nrow(yangambi_flowering)))


# Statistically stable subset

yangambi_df <- yangambi_flowering %>%
  filter(year >= 1938 & year <= 1956)


# Climate data:  NETCDF CLIMATE EXTRACTION & PROJECTION MATCHING

library(raster)
library(dplyr)

setwd("E:/CLIMATE DATA")

# Load CRU TS monthly climate grids
tmp_nc <- raster::brick(
  "cru_ts4.09.1901.2024.tmp.dat.nc",
  varname = "tmp"
)

pre_nc <- raster::brick(
  "cru_ts4.09.1901.2024.pre.dat.nc",
  varname = "pre"
)

sites <- data.frame(
  site      = "Yangambi",
  latitude  = 0.7623637,  
  longitude = 24.46295
)


# View the single-site object
print(sites)

# Study period
study_years <- 1938:1956

# Verify the number of monthly layers
nlayers(tmp_nc)
nlayers(pre_nc)

# CRU TS should have 12 monthly layers per year
stopifnot(nlayers(tmp_nc) == nlayers(pre_nc))

# Build the time index for every layer
full_time <- expand.grid(
  month = 1:12,
  year = 1901:2024
)

# expand.grid varies year fastest; reorder chronologically
full_time <- full_time[order(full_time$year, full_time$month), ]
rownames(full_time) <- NULL

# Check that the time index matches the raster layers
stopifnot(nrow(full_time) == nlayers(tmp_nc))

# Identify the layers for the study years
index <- which(full_time$year %in% study_years)

stopifnot(length(index) == length(study_years) * 12)

time_subset <- full_time[index, ]

# Extract all monthly values for Yangambi
tmp_all <- raster::extract(
  tmp_nc,
  sites[, c("longitude", "latitude")]
)

pre_all <- raster::extract(
  pre_nc,
  sites[, c("longitude", "latitude")]
)

# Ensure the extraction results are matrices
tmp_all <- as.matrix(tmp_all)
pre_all <- as.matrix(pre_all)

# Select the correct months
tmp_vals <- tmp_all[, index, drop = FALSE]
pre_vals <- pre_all[, index, drop = FALSE]

# Calculate annual climate summaries

climate_summary <- data.frame(
  site = "Yangambi",
  year = study_years,
  mean_temp = vapply(study_years, function(y) {
    mean(tmp_vals[1, time_subset$year == y], na.rm = TRUE)
  }, numeric(1)),
  total_precip = vapply(study_years, function(y) {
    sum(pre_vals[1, time_subset$year == y], na.rm = TRUE)
  }, numeric(1))
)

print(climate_summary)

# 8. MERGE DATA

library(dplyr)

yangambi_df <- yangambi_df %>%
  group_by(species, genus, family, year) %>%
  filter(sum(!is.na(doy)) >= 5) %>%
  summarise(
    n_obs = sum(!is.na(doy)),
    onset_doy = quantile(doy, 0.10, na.rm = TRUE, type = 7),
    end_doy = quantile(doy, 0.90, na.rm = TRUE, type = 7),
    .groups = "drop"
  ) %>%
  mutate(
    site = "Yangambi",
    duration_doy = end_doy - onset_doy
  ) %>%
  left_join(climate_summary, by = c("site", "year")) %>%
  mutate(
    year_c = year - 1938,
    precip_100 = total_precip / 100
  ) %>%
  select(
    site, family, genus, species,year,
    n_obs, onset_doy, end_doy, duration_doy,
    mean_temp, total_precip, year_c, precip_100
  )


species_year_coverage <- yangambi_data %>%
  group_by(species, genus, family) %>%
  summarise(
    n_years = n_distinct(year),
    first_year = min(year),
    last_year = max(year),
    .groups = "drop"
  ) %>%
  arrange(n_years)

summary(species_year_coverage$n_years)

table(species_year_coverage$n_years)


# INDEPENDENT TEMPORAL BOUNDS & DISTRIBUTIONS


# Apply specific historical subset ranges 
yangambi_df <- yangambi_df %>%
  filter(year >= 1938 & year <= 1956)

# --- FINAL QUALITY SANITY CONTROLS ---
print(paste("Total Yangambi analysis rows:", nrow(yangambi_data)))

print("Confirming DOY fields are structured as continuous numeric variables:")
print(is.numeric(yangambi_df$onset_doy))
print(paste("Unique value diversity:", length(unique(yangambi_data$onset_doy))))

# Generate visual distribution layout curve
plot(density(yangambi_df$onset_doy, na.rm = TRUE), 
     main = "Density of Onset Days (Yangambi)", xlab = "Day of Year")


# reconstruct the phylogeny using phylomaker
#  Remove duplicate species entries

cat("Species after removing duplicate entries:",
    nrow(sp.list_clean), "\n")

# Standardize family names to match the tree

sp.list_clean <- yangambi_df %>%
  mutate(
    family = case_when(
      genus %in% fabaceae_genera ~ "Fabaceae",
      genus == "Bersama" ~ "Francoaceae",
      genus == "Oncoba" ~ "Salicaceae",
      genus == "Rhaph iostylis" ~ family,
      genus == "Rhaphiostylis" ~ "Metteniusaceae",
      genus == "Thomandersia" ~ "Thomandersiaceae",
      tolower(family) == "leguminosae" ~ "Fabaceae",
      TRUE ~ family
    )
  )

# Inspect the updated classifications

sp.list_clean %>%
  distinct(genus, family) %>%
  arrange(family, genus) %>%
  print(n = Inf)


#  Remove taxa that failed to bind to the tree


failed_taxa <- c(
  "Dalhausia_africana",
  "Dalhousia_africana",
  "Dewevrei_bulabiata",
  "Mozoneurum_angolense",
  "Pithecolobium_samen",
  "Platysepalum_chevalieri",
  "Schefflerodendron_usamb",
  "Schefflerodendron_usambarense"
)

removed_taxa <- sp.list_clean %>%
  filter(species %in% failed_taxa)

cat("Failed-to-bind taxa found and removed:",
    nrow(removed_taxa), "\n")

print(removed_taxa, n = Inf)

sp.list_clean <- sp.list_clean %>%
  filter(!species %in% failed_taxa)


#  Verify the cleaned list


cat("Final number of species:", nrow(sp.list_clean), "\n")

cat("Remaining duplicated species:",
    anyDuplicated(sp.list_clean$species), "\n")


# Run phylo.maker again

phylo_result_yangambi <- phylo.maker(
  sp.list = sp.list_clean,
  tree = GBOTB.extended,
  nodes = nodes.info.1,
  scenarios = "S3"
)

yangambi_phy <- phylo_result_yangambi$scenario.3

# Extract the phylogenetic tree

# Plot tree

plot(
  yangambi_phy,
  cex = 0.5,
  no.margin = TRUE
)

# Match clean dataset to the new tree

yangambi_df$species_tree <- gsub(" ", "_", yangambi_df$species)

df_phy <- yangambi_df %>%
  filter(species_tree %in% yangambi_phy$tip.label)


# CREATE PHYLOGENETIC COVARIANCE MATRIX 

phy_cov <- vcv(yangambi_phy, corr = TRUE)

# PHYLOGENETIC SIGNAL TEST 

library(ape)
library(phytools)

# Species-level mean trait from clean data

library(dplyr)
library(ape)
library(phytools)

# Calculate mean flowering onset for each species

trait_df <- yangambi_df %>%
  group_by(species) %>%
  summarise(
    trait_value = mean(duration_doy, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(trait_value))

# Create a named trait vector

trait <- trait_df$trait_value
names(trait) <- trait_df$species

# Identify species shared by the trait data and tree
common_sp <- intersect(names(trait), yangambi_phy$tip.label)

cat("Species in trait data:", length(trait), "\n")
cat("Species in phylogeny:", length(yangambi_phy$tip.label), "\n")
cat("Matching species:", length(common_sp), "\n")

# Prune the tree to matching species
yangambi_phy <- ape::keep.tip(yangambi_phy, common_sp)

# Align trait values to the tree's tip order
trait_aligned <- trait[yangambi_phy$tip.label]

# Verify alignment
stopifnot(
  identical(names(trait_aligned), yangambi_phy$tip.label),
  !anyNA(trait_aligned)
)

cat("Aligned trait values:", length(trait_aligned), "\n")
cat("Tree tips:", length(yangambi_phy$tip.label), "\n")

# Calculate Blomberg's K with a randomization test

K_result <- phytools::phylosig(
  tree = yangambi_phy,
  x = trait_aligned,
  method = "K",
  test = TRUE,
  nsim = 999
)

print(K_result)


#Flowering onset 0.0640098 0.0750751 
#Flowering end  0.00828936 0.868869
#Flowering duration 0.0132213    0.680681 
 

##PART 2. Climate sensitivity 

##Standardization 

yangambi_df <- yangambi_df %>%
  mutate(
    temp_z   = as.numeric(scale(mean_temp)),
    precip_z = as.numeric(scale(total_precip))
  )

#Luki onset with interaction temp_z * precip_z

# Year as factor

yangambi_df$year <- as.factor(yangambi_df$year)



yangambi_onset <- phyr::pglmm(
  onset_doy ~ temp_z * precip_z +
    (1 | species__) +
    (1 | year),
  data = yangambi_df,
  family = "gaussian",
  cov_ranef = list(species = yangambi_phy),
  REML = FALSE
)

summary(yangambi_onset)

# summary(yangambi_onset)
# Linear mixed model fit by maximum likelihood
# 
# Call:onset_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
# -24322  48660  48702 
# 
# Random effects:
#              Variance  Std.Dev
# 1|species   7.673e+02   27.700
# 1|species__ 2.917e+01    5.401
# 1|year      2.573e+07 5072.831
# residual    3.770e+03   61.401
# 
# Fixed effects:
#                     Value Std.Error  Zscore Pvalue
# (Intercept)       88.3539 1201.4657  0.0735 0.9414
# temp_z            -9.4224 1344.0907 -0.0070 0.9944
# precip_z           1.9038 1199.4292  0.0016 0.9987
# temp_z:precip_z    4.5233 1248.8148  0.0036 0.9971


yangambi_onset_bobyqa <- phyr::pglmm(
  onset_doy ~ temp_z * precip_z +
    (1 | species__) +
    (1 | year),
  data = yangambi_df,
  family = "gaussian",
  cov_ranef = list(species = yangambi_phy),
  REML = FALSE,
  optimizer = "bobyqa",
  maxit = 5000,
  reltol = 1e-10,
  verbose = FALSE
)


summary(yangambi_onset_bobyqa)

# summary(yangambi_onset_bobyqa)
# Linear mixed model fit by maximum likelihood
# 
# Call:onset_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
# -24325  48665  48707 
# 
# Random effects:
#              Variance  Std.Dev
# 1|species   7.656e+02   27.670
# 1|species__ 2.782e+01    5.275
# 1|year      3.380e+07 5813.586
# residual    3.772e+03   61.418
# 
# Fixed effects:
#                     Value Std.Error  Zscore Pvalue
# (Intercept)       88.3397 1376.9059  0.0642 0.9488
# temp_z            -9.4185 1540.3601 -0.0061 0.9951
# precip_z           1.9024 1374.5745  0.0014 0.9989
# temp_z:precip_z    4.5208 1431.1716  0.0032 0.9975


##Without interaction

yangambi_onset_additive <- phyr::pglmm(
  onset_doy ~ temp_z + precip_z +
    (1 | species__) +
    (1 | year),
  data = yangambi_df,
  family = "gaussian",
  cov_ranef = list(species = yangambi_phy),
  REML = FALSE
)

summary(yangambi_onset_additive)

# summary(yangambi_onset_additive)
# Linear mixed model fit by maximum likelihood
# 
# Call:onset_doy ~ temp_z + precip_z
# 
# logLik    AIC    BIC 
# -24322  48658  48695 
# 
# Random effects:
#              Variance  Std.Dev
# 1|species   7.673e+02   27.701
# 1|species__ 2.917e+01    5.401
# 1|year      2.574e+07 5072.995
# residual    3.770e+03   61.401
# 
# Fixed effects:
#                 Value Std.Error  Zscore Pvalue
# (Intercept)   88.4769 1201.0263  0.0737 0.9413
# temp_z        -6.5605 1087.3714 -0.0060 0.9952
# precip_z       1.2977 1187.7366  0.0011 0.9991


#YANGAMBI_END PHENOLOGY

yangambi_end <- phyr::pglmm(
  end_doy ~ temp_z * precip_z +
    (1 | species__) +
    (1 | year),
  data = yangambi_df,
  family = "gaussian",
  cov_ranef = list(species = yangambi_phy)
)

summary(yangambi_end)

# summary(yangambi_end)
# Linear mixed model fit by maximum likelihood
# 
# Call:end_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
# -25370  50757  50799 
# 
# Random effects:
#              Variance Std.Dev
# 1|species   2.457e+03   49.56
# 1|species__ 1.532e+02   12.38
# 1|year      6.951e+07 8337.37
# residual    5.686e+03   75.41
# 
# Fixed effects:
#                      Value  Std.Error  Zscore Pvalue
# (Intercept)      261.93063 1974.66054  0.1326 0.8945
# temp_z            -0.25600 2209.05945 -0.0001 0.9999
# precip_z          -0.22047 1971.30331 -0.0001 0.9999
# temp_z:precip_z    1.18014 2052.47019  0.0006 0.9995
#End additive

yangambi_end_additive <- phyr::pglmm(
  end_doy ~ temp_z + precip_z +
    (1 | species__) +
    (1 | year),
  data = yangambi_df,
  family = "gaussian",
  cov_ranef = list(species = yangambi_phy)
)

summary(yangambi_end_additive)

# summary(yangambi_end_additive)
# Linear mixed model fit by restricted maximum likelihood
# 
# Call:end_doy ~ temp_z + precip_z
# 
# logLik    AIC    BIC 
# -25332  50679  50715 
# 
# Random effects:
#              Variance Std.Dev
# 1|species   2.447e+03   49.46
# 1|species__ 1.563e+02   12.50
# 1|year      6.963e+07 8344.26
# residual    5.691e+03   75.44
# 
# Fixed effects:
#                  Value  Std.Error  Zscore Pvalue
# (Intercept)  262.04590 1975.50425  0.1326 0.8945
# temp_z         0.49112 1788.54976  0.0003 0.9998
# precip_z      -0.37760 1953.63442 -0.0002 0.9998

#Adding bobyqa to stabilize the model

yangambi_end_bobyqa <- phyr::pglmm(
  end_doy ~ temp_z * precip_z +
    (1 | species__) +
    (1 | year),
  data = yangambi_data,
  family = "gaussian",
  cov_ranef = list(species = yangambi_phy),
  REML = FALSE,
  optimizer = "bobyqa",
  maxit = 5000,
  reltol = 1e-10,
  verbose = FALSE
)

summary(yangambi_end_bobyqa)

#The best model
# summary(yangambi_end_bobyqa)

# Linear mixed model fit by maximum likelihood
# 
# Call:end_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
# -25265  50545  50587 
# 
# Random effects:
#             Variance Std.Dev
# 1|species     2406.7   49.06
# 1|species__    154.4   12.43
# 1|year         295.0   17.18
# residual      5725.5   75.67
# 
# Fixed effects:
#                      Value  Std.Error  Zscore Pvalue    
# (Intercept)     262.018523  10.604120 24.7091 <2e-16 ***
# temp_z            0.341505   4.791189  0.0713 0.9432    
# precip_z         -0.046651   4.249068 -0.0110 0.9912    
# temp_z:precip_z   0.451040   4.465981  0.1010 0.9196    


##Duration

library(phyr)

yangambi_duration <- pglmm(
  duration_doy ~ temp_z * precip_z +
    (1 |species__) +
    (1 | year),
  data = yangambi_df,
  family = "gaussian",
  cov_ranef = list(species = yangambi_phy)
)

summary(yangambi_duration)

# summary(yangambi_duration)
# Linear mixed model fit by restricted maximum likelihood
# 
# Call:duration_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
# -25779  51574  51615 
# 
# Random effects:
#              Variance  Std.Dev
# 1|species   3.104e+03    55.71
# 1|species__ 2.407e+02    15.51
# 1|year      1.083e+08 10404.44
# residual    6.995e+03    83.63
# 
# Fixed effects:
#                     Value Std.Error  Zscore Pvalue
# (Intercept)      172.3099 2464.2341  0.0699 0.9443
# temp_z             9.5962 2756.7478  0.0035 0.9972
# precip_z          -2.3405 2460.0452 -0.0010 0.9992
# temp_z:precip_z   -3.5844 2561.3357 -0.0014 0.9989
##Additive

yangambi_duration_add <- pglmm (
  duration_doy ~ temp_z + precip_z +
    (1 | species__) +
    (1 | year),
  data = yangambi_df,
  family = "gaussian",
  cov_ranef = list(species = yangambi_phy)
)


summary(yangambi_duration_add)

# summary(yangambi_duration_add)
# Linear mixed model fit by restricted maximum likelihood
# 
# Call:duration_doy ~ temp_z + precip_z
# 
# logLik    AIC    BIC 
# -25792  51597  51634 
# 
# Random effects:
#              Variance  Std.Dev
# 1|species   3.105e+03    55.73
# 1|species__ 2.408e+02    15.52
# 1|year      1.083e+08 10406.76
# residual    6.992e+03    83.62
# 
# Fixed effects:
#                 Value Std.Error  Zscore Pvalue
# (Intercept)  172.2102 2463.8007  0.0699 0.9443
# temp_z         7.3286 2230.6367  0.0033 0.9974
# precip_z      -1.8605 2436.5264 -0.0008 0.9994


yangambi_duration_bobyqa <- phyr::pglmm(
  duration_doy ~ temp_z * precip_z +
    (1 | species__) +
    (1 | year),
  data = yangambi_data,
  family = "gaussian",
  cov_ranef = list(species = yangambi_phy),
  REML = FALSE,
  optimizer = "bobyqa",
  maxit = 5000,
  reltol = 1e-10,
  verbose = FALSE
)

summary(yangambi_duration_bobyqa)


# summary(yangambi_duration_bobyqa)
# Linear mixed model fit by maximum likelihood
# 
# Call:duration_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
# -25798  51613  51655 
# 
# Random effects:
#             Variance Std.Dev
# 1|species       3111   55.77
# 1|species__      240   15.49
# 1|year       3825848 1955.98
# residual        6986   83.58
# 
# Fixed effects:
#                    Value Std.Error  Zscore Pvalue
# (Intercept)     172.2833  463.4125  0.3718 0.7101
# temp_z            9.5963  518.2558  0.0185 0.9852
# precip_z         -2.3417  462.4768 -0.0051 0.9960
# temp_z:precip_z  -3.5836  481.5194 -0.0074 0.9941



yangambi_duration_bob_add <- phyr::pglmm(
  duration_doy ~ temp_z + precip_z +
    (1 | species__) +
    (1 | year),
  data = yangambi_df,
  family = "gaussian",
  cov_ranef = list(species = yangambi_phy),
  REML = FALSE,
  optimizer = "bobyqa",
  maxit = 5000,
  reltol = 1e-10,
  verbose = FALSE
)

summary(yangambi_duration_bob_add)

# summary(yangambi_duration_bob_add)
# Linear mixed model fit by maximum likelihood
# 
# Call:duration_doy ~ temp_z + precip_z
# 
# logLik    AIC    BIC 
# -25800  51615  51652 
# 
# Random effects:
#              Variance Std.Dev
# 1|species      3131.3   55.96
# 1|species__     234.1   15.30
# 1|year      4701905.6 2168.39
# residual       6984.8   83.57
# 
# Fixed effects:
#                Value Std.Error  Zscore Pvalue
# (Intercept) 172.0679  513.4976  0.3351 0.7376
# temp_z        7.3292  464.7853  0.0158 0.9874
# precip_z     -1.8634  507.6849 -0.0037 0.9971
