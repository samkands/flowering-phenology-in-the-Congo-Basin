

luki_raw <- read.csv("E:\\PAPER 2\\RESULTS FLOWERING\\Phenological data\\luki_raw.csv")
yangambi_df <- read.csv ("E:\\PAPER 2\\RESULTS FLOWERING\\Phenological data\\yangambi_df.csv")

# LOAD MODELS

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
library(phyr)


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
    trait_value = mean(duration_doy, na.rm = TRUE),
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


##TEMPORAL TRENDS

#A. LUKI

install.packages("phyr")
library(phyr)
# A. LUKI

# 1. Onset 

luki_onset <- phyr::pglmm(
  onset_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__),
  data = luki_raw,
  family = "gaussian",
  cov_ranef = list(species = phy_cov),
  REML = FALSE
)

summary(luki_onset)


# summary(luki_onset_no_year2)
# Linear mixed model fit by maximum likelihood
# 
# Call:onset_doy ~ year_c
# 
# logLik    AIC    BIC 
#  -6343  12699  12726 
# 
# Random effects:
#                       Variance  Std.Dev
# 1|species            5.748e+04 239.7511
# 1|species__          4.681e+05 684.1445
# 0 + year_c|species   6.158e+01   7.8475
# 0 + year_c|species__ 1.716e-01   0.4142
# residual             5.334e+03  73.0357
# 
# Fixed effects:
#                Value Std.Error  Zscore    Pvalue    
# (Intercept) 132.2142  366.3966  0.3609    0.7182    
# year_c       -7.4779    1.1032 -6.7783 1.216e-11 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# 2. Termination 

luki_end <- phyr::pglmm(
  end_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__),
  data = luki_raw,
  family = "gaussian",
  cov_ranef = list(species = phy_cov),
  REML = FALSE
)

summary(luki_end)


# summary(luki_end_no_year2)
# Linear mixed model fit by maximum likelihood
# 
# Call:end_doy ~ year_c
# 
# logLik    AIC    BIC 
#  -6295  12604  12631 
# 
# Random effects:
#                       Variance   Std.Dev
# 1|species            7.094e+04  266.3482
# 1|species__          2.876e+06 1696.0236
# 0 + year_c|species   1.715e+01    4.1415
# 0 + year_c|species__ 3.013e-02    0.1736
# residual             3.952e+03   62.8682
# 
# Fixed effects:
#                 Value Std.Error  Zscore    Pvalue    
# (Intercept) 270.10273 904.39954  0.2987    0.7652    
# year_c       -3.12749   0.76655 -4.0800 4.504e-05 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# 3. Duration 

luki_duration <- phyr::pglmm(
  duration_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__),
  data = luki_raw,
  family = "gaussian",
  cov_ranef = list(species = phy_cov),
  REML = FALSE
)

summary(luki_duration)


# summary(luki_duration_no_year2)
# Linear mixed model fit by maximum likelihood
# 
# Call:duration_doy ~ year_c
# 
# logLik    AIC    BIC 
#  -6454  12922  12948 
# 
# Random effects:
#                       Variance  Std.Dev
# 1|species            5.183e+04 227.6705
# 1|species__          3.268e+05 571.6562
# 0 + year_c|species   1.086e+02  10.4235
# 0 + year_c|species__ 5.912e-01   0.7689
# residual             7.008e+03  83.7155
# 
# Fixed effects:
#                Value Std.Error Zscore   Pvalue   
# (Intercept) 140.0052  306.6522 0.4566 0.647987   
# year_c        3.9040    1.4162 2.7566 0.005841 **
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# B. YANGAMBI


library(dplyr)
library(ape)
library(phytools)
library(V.PhyloMaker)


# Create cleaned species list
sp.list <- yangambi_df %>%
  dplyr::select(species, genus, family) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    family = dplyr::case_when(
      genus %in% fabaceae_genera ~ "Fabaceae",
      genus == "Bersama"         ~ "Francoaceae",
      genus == "Oncoba"          ~ "Salicaceae",
      genus == "Rhaphiostylis"   ~ "Metteniusaceae",
      genus == "Thomandersia"    ~ "Thomandersiaceae",
      tolower(family) == "leguminosae" ~ "Fabaceae",
      TRUE ~ as.character(family)
    )
  ) %>%
  as.data.frame()

# Replace spaces with underscores in species names
sp.list$species <- gsub(" ", "_", sp.list$species)

rownames(sp.list) <- NULL

# Check cleaned taxonomy
sp.list %>%
  distinct(genus, family) %>%
  arrange(family, genus)

# Count species per family
sp.list %>%
  count(family, sort = TRUE)




# Load V.PhyloMaker reference data
data(GBOTB.extended)
data(nodes.info.1)

nodes.info.1 <- nodes.info.1 %>%
  as.data.frame()

# Generate phylogeny
phylo_result <- V.PhyloMaker::phylo.maker(
  sp.list = sp.list,
  tree = GBOTB.extended,
  nodes = nodes.info.1,
  scenarios = "S3"
)

# Extract phylogenetic tree
yang_phy <- phylo_result$scenario.3

# Check number of tips
cat("Number of species in phylogeny:",
    length(yang_phy$tip.label), "\n")

# Save phylogenetic tree
ape::write.tree(
  yang_phy,
  file = "yangambi_phylogenetic_tree.newick"
)

# Save cleaned species list
write.csv(
  sp.list,
  file = "yangambi_cleaned_species_list.csv",
  row.names = FALSE
)

# Plot tree
plot(
  yang_phy,
  cex = 0.3,
  no.margin = TRUE
)




# Create standardized species names
yangambi_df$species_tree <- gsub(
  " ",
  "_",
  yangambi_df$species
)

# Identify species represented in the phylogeny
yangambi_phy <- yangambi_df %>%
  filter(species_tree %in% yang_phy$tip.label)

# Check matching
cat("Species in Yangambi dataset:",
    dplyr::n_distinct(yangambi_df$species_tree), "\n")

cat("Species in phylogeny:",
    length(yang_phy$tip.label), "\n")

cat("Species matched:",
    dplyr::n_distinct(yangambi_phy$species_tree), "\n")

cat("Species not matched:",
    length(setdiff(
      unique(yangambi_df$species_tree),
      yang_phy$tip.label
    )), "\n")



# Create phylogenetic covariance matrix
yang_phy_cov <- ape::vcv(
  yang_phy,
  corr = TRUE
)

# Check dimensions
dim(yang_phy_cov)

# Check species names
head(rownames(yang_phy_cov))

# Check matching between tree and covariance matrix
stopifnot(
  identical(
    rownames(yang_phy_cov),
    yang_phy$tip.label
  )
)


# 1. Calculate mean flowering termination for each species
trait_df <- yangambi_df %>%
  group_by(species_tree = gsub(" ", "_", species)) %>%
  summarise(
    trait_value = mean(duration_doy, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(trait_value))

# 2. Create a named trait vector
trait <- trait_df$trait_value
names(trait) <- trait_df$species_tree

# 3. Identify species shared by trait data and phylogeny
common_sp <- intersect(
  names(trait),
  yang_phy$tip.label
)

cat("Species in trait data:", length(trait), "\n")
cat("Species in phylogeny:", length(yang_phy$tip.label), "\n")
cat("Matching species:", length(common_sp), "\n")

# 4. Prune tree to matching species
yang_phy_sub <- ape::keep.tip(
  yang_phy,
  common_sp
)

# 5. Align trait values to tree tip order
trait_aligned <- trait[yang_phy_sub$tip.label]

# 6. Verify alignment
stopifnot(
  identical(
    names(trait_aligned),
    yang_phy_sub$tip.label
  ),
  !anyNA(trait_aligned)
)

cat("Aligned trait values:", length(trait_aligned), "\n")
cat("Tree tips:", length(yang_phy_sub$tip.label), "\n")

# 7. Calculate Blomberg's K
K_result <- phytools::phylosig(
  tree = yang_phy_sub,
  x = trait_aligned,
  method = "K",
  test = TRUE,
  nsim = 999
)

# Display results
print(K_result)


#Flowering onset 0.0640098 0.0750751 
#Flowering end  0.00828936 0.868869
#Flowering duration 0.0132213    0.680681 




# 1. Onset 

yangambi_onset_year <- phyr::pglmm(
  onset_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__),
  data = yangambi_df,
  family = "gaussian",
  cov_ranef = list(species = yang_phy_cov),
  REML = FALSE
)

summary(yangambi_onset_year)

# summary(yangambi_onset_no_year2)
# Linear mixed model fit by maximum likelihood
# 
# Call:onset_doy ~ year_c
# 
# logLik    AIC    BIC 
# -25844  51702  51739 
# 
# Random effects:
#                       Variance  Std.Dev
# 1|species            7.091e+05 842.0674
# 1|species__          1.470e+04 121.2542
# 0 + year_c|species   1.700e+01   4.1230
# 0 + year_c|species__ 4.968e-01   0.7049
# residual             3.176e+03  56.3539
# 
# Fixed effects:
#                 Value Std.Error  Zscore    Pvalue    
# (Intercept) 123.97456 102.40686  1.2106     0.226    
# year_c       -3.75364   0.68396 -5.4881 4.063e-08 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# 2. Termination 

yangambi_end_year <- phyr::pglmm(
  end_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__),
  data = yangambi_df,
  family = "gaussian",
  cov_ranef = list(species = yang_phy_cov),
  REML = FALSE
)

summary(yangambi_end_year)

# summary(yangambi_end_no_year2)
# Linear mixed model fit by maximum likelihood
# 
# Call:end_doy ~ year_c
# 
# logLik    AIC    BIC 
# -27038  54090  54126 
# 
# Random effects:
#                       Variance   Std.Dev
# 1|species            3.517e+06 1.875e+03
# 1|species__          2.738e+04 1.655e+02
# 0 + year_c|species   2.270e+01 4.765e+00
# 0 + year_c|species__ 4.920e-03 7.014e-02
# residual             4.724e+03 6.873e+01
# 
# Fixed effects:
#                 Value Std.Error  Zscore    Pvalue    
# (Intercept) 276.52117 158.17291  1.7482   0.08043 .  
# year_c       -2.65720   0.40349 -6.5855 4.532e-11 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#3. Duration 

yangambi_duration_year <- phyr::pglmm(
  duration_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__),
  data = yangambi_df,
  family = "gaussian",
  cov_ranef = list(species = yang_phy_cov),
  REML = FALSE
)

summary(yangambi_duration_year)


# Model Structure: duration_doy ~ year_c (No 1|year random intercept)
# Fit Criteria Metrics:
# AIC: 52575.75
# BIC: 52612.34

# Fixed Effects Table:
#                 Value   Std.Error    Zscore      Pvalue
# (Intercept) 161.443741   64.738849  2.493769  0.01263947 *
# year_c        1.274795    1.260095  1.011666  0.31169784 (Not Significant)


# summary(yangambi_duration_no_year2)
# Linear mixed model fit by maximum likelihood
# 
# Call:duration_doy ~ year_c
# 
# logLik    AIC    BIC 
# -25732  51478  51514 
# 
# Random effects:
#                       Variance Std.Dev
# 1|species            2379.7793 48.7830
# 1|species__           133.4341 11.5514
# 0 + year_c|species     11.2116  3.3484
# 0 + year_c|species__    0.9186  0.9584
# residual             6953.9349 83.3903
# 
# Fixed effects:
#                 Value Std.Error  Zscore Pvalue    
# (Intercept) 166.96873  10.01970 16.6640 <2e-16 ***
# year_c        0.97037   0.86434  1.1227 0.2616    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


save.image("E:/PAPER 2/flowering_models_year_comparison.RData")
