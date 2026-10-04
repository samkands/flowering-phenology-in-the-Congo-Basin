
##FLOWERING PHENOLOGY IN THE CONGO BASIN

install.packages("ggpubr")
library(ggpubr)
install.packages("cairoDevice")
library(cairoDevice)

library(dplyr)
library(glmmTMB)
library(ggplot2)
library(patchwork)
if(!require(VennDiagram)) install.packages("VennDiagram")
library(VennDiagram)
library(V.PhyloMaker)
library(ape)
library(patchwork)
library(phytools)

setwd("E:/PAPER 2/PGLMM")

flowering <- readRDS("luki_data.csv")

luki_data <- read.csv ("luki_data.csv")
yangambi_data <- read.csv("yangambi_data.csv")

table(flowering$year)


#Separet Luki and Yangambi

# 1. Create a dataset containing ONLY Luki

luki_data <- subset(flowering, site == "Luki")

# 2. Create a dataset containing ONLY Yangambi

yangambi_data <- subset(flowering, site == "Yangambi")

write.csv(luki_data, file = "luki_data.csv", row.names = FALSE)
write.csv(yangambi_data, file = "yangambi_data.csv", row.names = FALSE)


# reconstruct the phylogeny using phylomaker

sp.list <- luki_data %>%
  select(species, genus, family) %>%
  distinct()

# Replace spaces with underscores in species names
sp.list$species <- gsub(" ", "_", sp.list$species)

# 2. LOAD V.PHYLOMAKER DATA & UPDATE TAXONOMY

data(GBOTB.extended)
data(nodes.info.1)

# Update family names in clean data to match the megatree
sp.list$family[sp.list$genus == "Aptandra"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Cordia"] <- "Boraginaceae"
sp.list$family[sp.list$genus == "Diogoa"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Heisteria"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Octoknema"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Okoubaka"] <- "Santalaceae"
sp.list$family[sp.list$genus == "Ongokea"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Strombosia"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Strombosiopsis"] <- "Olacaceae"

# Replace NA with the correct genus name in the genus columns for 2 species
sp.list$genus[sp.list$species == "Chrysophyllum_claussensii"] <- "Chrysophyllum"
sp.list$genus[sp.list$species == "Randia_acuminata"] <- "Randia"


# 3. BUILD PHYLOGENY

phylo_result <- phylo.maker(
  sp.list = sp.list,
  tree = GBOTB.extended,
  nodes = nodes.info.1,
  scenarios = "S3"
)

# Extract the phylogenetic tree

phy <- phylo_result$scenario.3

# Plot tree
plot(
  phy,
  cex = 0.5,
  no.margin = TRUE
)

# Match clean dataset to the new tree

luki_data$species_tree <- gsub(" ", "_", luki_data$species)

df_phy <- luki_data %>%
  filter(species_tree %in% phy$tip.label)


# 4. CREATE PHYLOGENETIC COVARIANCE MATRIX 

phy_cov <- vcv(phy, corr = TRUE)



# 5. PHYLOGENETIC SIGNAL TEST 

library(ape)
library(phytools)

# Species-level mean trait from clean data

trait_df <- aggregate(duration_doy ~ species,
                      data = luki_data,
                      FUN = mean,
                      na.rm = TRUE)

# Create named vector

trait <- trait_df$duration_doy
names(trait) <- trait_df$species

# Keep only species present in both trait data and phylogeny

common_sp <- intersect(names(trait), phy$tip.label)

trait <- trait[common_sp]
phy_sub <- drop.tip(phy, setdiff(phy$tip.label, common_sp))

# Ensure order matches

trait <- trait[phy_sub$tip.label]

# Blomberg's K with randomization test

K_result <- phylosig(
  tree = phy_sub,
  x = trait,
  method = "K",
  test = TRUE,
  nsim = 999
)

K_result 

#Phylogenetic signal onset: K : 0.140331    P-value (based on 999 randomizations) :0.946947  
#Phylogenetic signal K end_doy : 0.139425  P-value (based on 999 randomizations) : 0.950951  
#Phylogenetic signal K duration: 0.13987  P-value (based on 999 randomizations) : 0.958959  

##PART 2. Climate sensitivity 

##Standardization 

luki_data <- luki_data %>%
  mutate(
    temp_z   = as.numeric(scale(mean_temp)),
    precip_z = as.numeric(scale(total_precip))
  )


install.packages("phyr")
library(phyr)

#Luki onset with interaction temp_z * precip_z

phylo_onset <- pglmm(
  onset_doy ~ temp_z * precip_z +
    (1 | species) +
    (1 | year),
  data = luki_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)


summary(phylo_onset)

# summary(phylo_onset)
# Linear mixed model fit by restricted maximum likelihood
# 
# Call:onset_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
#  -1269   2553   2582 
# 
# Random effects:
#           Variance Std.Dev
# 1|species   1.8737  1.3688
# 1|year      0.1935  0.4399
# residual    0.2014  0.4488
# 
# Fixed effects:
#                     Value Std.Error   Zscore Pvalue    
# (Intercept)     40.862123  0.267157 152.9517 <2e-16 ***
# temp_z          -0.065344  0.267116  -0.2446 0.8067    
# precip_z         0.185880  0.259191   0.7172 0.4733    
# temp_z:precip_z  0.080414  0.257538   0.3122 0.7549    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##Onset without interaction

phylo_onset_additive <- pglmm(
  onset_doy ~ temp_z + precip_z +
    (1 | species) +
    (1 | year),
  data = luki_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)

summary(phylo_onset_additive)

# Linear mixed model fit by restricted maximum likelihood
# 
# Call:onset_doy ~ temp_z + precip_z
# 
# logLik    AIC    BIC 
#  -1272   2556   2581 
# 
# Random effects:
#           Variance Std.Dev
# 1|species   1.8728  1.3685
# 1|year      0.1716  0.4142
# residual    0.2014  0.4488
# 
# Fixed effects:
#                 Value Std.Error   Zscore Pvalue    
# (Intercept) 40.799473  0.170633 239.1062 <2e-16 ***
# temp_z      -0.013018  0.195944  -0.0664 0.9470    
# precip_z     0.233224  0.198030   1.1777 0.2389    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Luki end with interaction temp_z * precip_z
phylo_end <- pglmm(
  end_doy ~ temp_z * precip_z +
    (1 | species) +
    (1 | year),
  data = luki_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)


summary(phylo_end)


# summary(phylo_end)
# Linear mixed model fit by restricted maximum likelihood
# 
# Call:end_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
#  -1292   2598   2627 
# 
# Random effects:
#           Variance Std.Dev
# 1|species   1.9848  1.4088
# 1|year      0.2104  0.4587
# residual    0.2069  0.4548
# 
# Fixed effects:
#                     Value Std.Error    Zscore Pvalue    
# (Intercept)     333.82634   0.27787 1201.3786 <2e-16 ***
# temp_z            0.12299   0.27849    0.4416 0.6588    
# precip_z          0.32008   0.27023    1.1845 0.2362    
# temp_z:precip_z   0.12149   0.26851    0.4524 0.6509    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




# LUKI END MODEL (No Interaction Term)

phylo_end_additive <- pglmm(
  end_doy ~ temp_z + precip_z +
    (1 | species) +
    (1 | year),
  data = luki_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)

summary(phylo_end_additive)

# Linear mixed model fit by restricted maximum likelihood
# 
# Call:end_doy ~ temp_z + precip_z
# 
# logLik    AIC    BIC 
#  -1295   2601   2626 
# 
# Random effects:
#           Variance Std.Dev
# 1|species   1.9846  1.4088
# 1|year      0.1896  0.4355
# residual    0.2069  0.4548
# 
# Fixed effects:
#                 Value Std.Error    Zscore  Pvalue    
# (Intercept) 333.73170   0.17765 1878.5806 < 2e-16 ***
# temp_z        0.20205   0.20593    0.9811 0.32653    
# precip_z      0.39161   0.20813    1.8816 0.05989 .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Luki duration with interaction temp_z * precip_z

phylo_duration <- pglmm(
  duration_doy ~ temp_z * precip_z +
    (1 | species) +
    (1 | year),
  data = luki_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)


summary(phylo_duration)

# Linear mixed model fit by restricted maximum likelihood
# 
# Call:duration_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
#  -2310   4634   4663 
# 
# Random effects:
#           Variance Std.Dev
# 1|species  7.69205  2.7735
# 1|year     0.01837  0.1355
# residual   0.81257  0.9014
# 
# Fixed effects:
#                      Value  Std.Error    Zscore  Pvalue    
# (Intercept)     292.964604   0.249337 1174.9766 < 2e-16 ***
# temp_z            0.188147   0.094260    1.9960 0.04593 *  
# precip_z          0.134314   0.091453    1.4687 0.14192    
# temp_z:precip_z   0.040788   0.090888    0.4488 0.65360    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##Luki DURATION MODEL (No Interaction Term)

phylo_duration_additive <- pglmm(
  duration_doy ~ temp_z + precip_z +
    (1 | species) +
    (1 | year),
  data = luki_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)


summary(phylo_onset)
summary(phylo_end)
summary(phylo_duration_additive)

# Linear mixed model fit by restricted maximum likelihood
# 
# Call:duration_doy ~ temp_z + precip_z
# 
# logLik    AIC    BIC 
#  -2312   4635   4660 
# 
# Random effects:
#           Variance Std.Dev
# 1|species  7.68360  2.7719
# 1|year     0.01597  0.1264
# residual   0.81266  0.9015
# 
# Fixed effects:
#                  Value  Std.Error    Zscore    Pvalue    
# (Intercept) 292.932921   0.238501 1228.2270 < 2.2e-16 ***
# temp_z        0.214429   0.070038    3.0616  0.002201 ** 
# precip_z      0.158216   0.070560    2.2423  0.024942 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##-------------------------------------------------------------
##YANGAMBI

yangambi_data <- subset(flowering, site == "Yangambi")

# reconstruct the phylogeny using phylomaker

sp.list <- yangambi_data %>%
  select(species, genus, family) %>%
  distinct()

# Replace spaces with underscores in species names
sp.list$species <- gsub(" ", "_", sp.list$species)

# 2. LOAD V.PHYLOMAKER DATA & UPDATE TAXONOMY

data(GBOTB.extended)
data(nodes.info.1)

# Update family names in clean data to match the megatree
sp.list$family[sp.list$genus == "Aptandra"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Cordia"] <- "Boraginaceae"
sp.list$family[sp.list$genus == "Diogoa"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Heisteria"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Octoknema"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Okoubaka"] <- "Santalaceae"
sp.list$family[sp.list$genus == "Ongokea"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Strombosia"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Strombosiopsis"] <- "Olacaceae"

# Replace NA with the correct genus name in the genus columns for 2 species
sp.list$genus[sp.list$species == "Chrysophyllum_claussensii"] <- "Chrysophyllum"
sp.list$genus[sp.list$species == "Randia_acuminata"] <- "Randia"


# BUILD PHYLOGENY

phylo_result <- phylo.maker(
  sp.list = sp.list,
  tree = GBOTB.extended,
  nodes = nodes.info.1,
  scenarios = "S3"
)

# Extract the phylogenetic tree

phy_yang <- phylo_result$scenario.3

# Plot tree
plot(
  phy_yang,
  cex = 0.5,
  no.margin = TRUE
)

# Match clean dataset to the new tree

yangambi_data$species_tree <- gsub(" ", "_", yangambi_data$species)

df_phy <- yangambi_data %>%
  filter(species_tree %in% phy$tip.label)


# CREATE PHYLOGENETIC COVARIANCE MATRIX 

phy_cov <- vcv(phy_yang, corr = TRUE)



# PHYLOGENETIC SIGNAL TEST 

library(ape)
library(phytools)

# Species-level mean trait from clean data

trait_df <- aggregate(end_doy ~ species,
                      data = yangambi_data,
                      FUN = mean,
                      na.rm = TRUE)

# Create named vector

trait <- trait_df$end_doy
names(trait) <- trait_df$species

# Keep only species present in both trait data and phylogeny

common_sp <- intersect(names(trait), phy$tip.label)

trait <- trait[common_sp]
phy_sub <- drop.tip(phy_yang, setdiff(phy_yang$tip.label, common_sp))

# Ensure order matches

trait <- trait[phy_sub$tip.label]

# Blomberg's K with randomization test

K_result <- phylosig(
  tree = phy_sub,
  x = trait,
  method = "K",
  test = TRUE,
  nsim = 999
)

K_result 

#Phylogenetic signal onset: K : 0.128618      P-value (based on 999 randomizations) : 0.003003  
#Phylogenetic signal K end_doy : 0.12619   P-value (based on 999 randomizations) : 0.004004   
#Phylogenetic signal K duration: 0.127505  P-value (based on 999 randomizations) : 0.001001   

##Climate sensitivity 

##Standardization 

yangambi_data <- yangambi_data %>%
  mutate(
    temp_z   = as.numeric(scale(mean_temp)),
    precip_z = as.numeric(scale(total_precip))
  )


##Yangambi onset with interaction (temp_z * precip_z)
yangambi_onset <- pglmm(
  onset_doy ~ temp_z * precip_z +
    (1 | species) +
    (1 | year),
  data = yangambi_data,
  family = "gaussian",
  cov_ranef = list(species = phy_yang)
)

summary(yangambi_onset)

# Linear mixed model fit by restricted maximum likelihood
# 
# Call:onset_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
#  -9197  18409  18448 
# 
# Random effects:
#           Variance Std.Dev
# 1|species   1.0378  1.0187
# 1|year      0.2621  0.5120
# residual    0.9621  0.9808
# 
# Fixed effects:
#                     Value Std.Error   Zscore Pvalue    
# (Intercept)     42.296348  0.134183 315.2140 <2e-16 ***
# temp_z           0.098275  0.147821   0.6648 0.5062    
# precip_z         0.200754  0.141348   1.4203 0.1555    
# temp_z:precip_z -0.058684  0.091333  -0.6425 0.5205    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##Yangambi onset without interaction

yangambi_onset_mode <- pglmm(
  onset_doy ~ temp_z + precip_z +
    (1 | species) +
    (1 | year),
  data = yangambi_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)

# Linear mixed model fit by restricted maximum likelihood
# 
# Call:onset_doy ~ temp_z + precip_z
# 
# logLik    AIC    BIC 
#  -9201  18413  18447 
# 
# Random effects:
#           Variance Std.Dev
# 1|species   1.0338  1.0168
# 1|year      0.2497  0.4997
# residual    0.9624  0.9810
# 
# Fixed effects:
#                 Value Std.Error   Zscore Pvalue    
# (Intercept) 42.164216  0.193683 217.6968 <2e-16 ***
# temp_z       0.067214  0.153668   0.4374 0.6618    
# precip_z     0.394872  0.289863   1.3623 0.1731    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##Yangambi end with interaction temp_z * precip_z

yangambi_end <- pglmm(
  end_doy ~ temp_z * precip_z +
    (1 | species) +
    (1 | year),
  data = yangambi_data,
  family = "gaussian",
  cov_ranef = list(species = phy_yang)
)


summary(yangambi_end)

# summary(yangambi_end)
# Linear mixed model fit by restricted maximum likelihood
# 
# Call:end_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
# -10008  20030  20069 
# 
# Random effects:
#           Variance Std.Dev
# 1|species    1.347  1.1605
# 1|year       0.485  0.6964
# residual     1.257  1.1213
# 
# Fixed effects:
#                       Value   Std.Error    Zscore  Pvalue    
# (Intercept)     329.2719858   0.1793805 1835.6059 < 2e-16 ***
# temp_z           -0.3718757   0.2002623   -1.8569 0.06332 .  
# precip_z         -0.2000233   0.1918123   -1.0428 0.29704    
# temp_z:precip_z  -0.0021401   0.1239015   -0.0173 0.98622    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


## Yangambi end without interaction

yangambi_end_mode <- pglmm(
  end_doy ~ temp_z + precip_z +
    (1 | species) +
    (1 | year),
  data = yangambi_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)

# yangambi_end_mode 
# Linear mixed model fit by restricted maximum likelihood
# 
# Call:end_doy ~ temp_z + precip_z
# 
# logLik    AIC    BIC 
# -10011  20035  20068 
# 
# Random effects:
#           Variance Std.Dev
# 1|species   1.3484   1.161
# 1|year      0.4623   0.680
# residual    1.2570   1.121
# 
# Fixed effects:
#                 Value Std.Error    Zscore  Pvalue    
# (Intercept) 329.56758   0.26120 1261.7462 < 2e-16 ***
# temp_z       -0.43850   0.20827   -2.1055 0.03525 *  
# precip_z     -0.42672   0.39352   -1.0844 0.27820    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Yangambi duration with interaction (temp_z*precip_z)
yangambi_duration <- pglmm(
  duration_doy ~ temp_z * precip_z +
    (1 | species) +
    (1 | year),
  data = yangambi_data,
  family = "gaussian",
  cov_ranef = list(species = phy_yang)
)


summary(yangambi_duration)

# summary(yangambi_duration)
# Linear mixed model fit by restricted maximum likelihood
# 
# Call:duration_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
# -13810  27634  27673 
# 
# Random effects:
#           Variance Std.Dev
# 1|species    4.757   2.181
# 1|year       1.370   1.170
# residual     4.418   2.102
# 
# Fixed effects:
#                      Value  Std.Error   Zscore Pvalue    
# (Intercept)     286.976065   0.304576 942.2151 <2e-16 ***
# temp_z           -0.469084   0.337358  -1.3905 0.1644    
# precip_z         -0.400424   0.322800  -1.2405 0.2148    
# temp_z:precip_z   0.056228   0.208553   0.2696 0.7875    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##yangambi duration without interaction

yangambi_duration_mode <- pglmm(
  duration_doy ~ temp_z + precip_z +
    (1 | species) +
    (1 | year),
  data = yangambi_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)



# yangambi_duration_mode 
# Linear mixed model fit by restricted maximum likelihood
# 
# Call:duration_doy ~ temp_z + precip_z
# 
# logLik    AIC    BIC 
# -13814  27640  27674 
# 
# Random effects:
#           Variance Std.Dev
# 1|species    4.705   2.169
# 1|year       1.329   1.153
# residual     4.423   2.103
# 
# Fixed effects:
#                 Value Std.Error   Zscore Pvalue    
# (Intercept) 287.40423   0.44493 645.9574 <2e-16 ***
# temp_z       -0.50524   0.35389  -1.4277 0.1534    
# precip_z     -0.82120   0.66804  -1.2293 0.2190    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


### PLOTS FOR CLIMATE SENDITIVITY
### CONSOLIDATED REPLICATED PHENOLOGY FOREST PLOTS

library(ggplot2)
library(dplyr)
library(grid)

# 1. Workspace Configuration
setwd("E:/PAPER 2/PGLMM")
output_dir <- "E:/PAPER 2/PGLMM/"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


# COMPREHENSIVE REPOSITORY DATA FRAME (Directly from your verified summaries)

master_effects <- data.frame(
  Site     = c(rep("Luki", 9), rep("Yangambi", 9)),
  Metric   = c(rep("Onset", 3), rep("End", 3), rep("Duration", 3), rep("Onset", 3), rep("End", 3), rep("Duration", 3)),
  Term     = rep(c("Temperature", "Precipitation", "Temp × Precip"), 6),
  Estimate = c(
    # Luki: Onset, End, Duration
    -0.065344, 0.185880, 0.080414,   0.122990, 0.320080, 0.121490,   0.188147, 0.134314, 0.040788,
    # Yangambi: Onset, End, Duration
    0.098275, 0.200754, -0.058684,  -0.371876, -0.200023, -0.002140,  -0.469084, -0.400424, 0.056228
  ),
  StdError = c(
    # Luki: Onset, End, Duration
    0.267116, 0.259191, 0.257538,   0.278490, 0.270230, 0.268510,   0.094260, 0.091453, 0.090888,
    # Yangambi: Onset, End, Duration
    0.147821, 0.141348, 0.091333,   0.200262, 0.191812, 0.123901,   0.337358, 0.322800, 0.208553
  ),
  PValue   = c(
    # Luki: Onset, End, Duration
    0.806700, 0.473300, 0.754900,   0.658800, 0.236200, 0.650900,   0.045930, 0.141920, 0.653600,
    # Yangambi: Onset, End, Duration
    0.506200, 0.155500, 0.520500,   0.063320, 0.297040, 0.986220,   0.164400, 0.214800, 0.787500
  )
)


# REPLICATION PLOTTING ENGINE

build_forest_plot <- function(target_metric) {
  
  # Calculate exact row sizes from your raw source objects
  true_luki_n     = nrow(luki_data)
  true_yangambi_n = nrow(yangambi_data)
  
  plot_data <- master_effects %>%
    filter(Metric == target_metric) %>%
    mutate(
      Term = factor(Term, levels = c("Temp × Precip", "Precipitation", "Temperature")),
      Sig_Label = case_when(
        PValue < 0.05                 ~ "*",
        PValue >= 0.05 & PValue < 0.1 ~ "†",
        TRUE                          ~ ""
      ),
      Label_Pos = Estimate + (1.96 * StdError) + 0.04
    )
  
  legend_labels <- c(
    "Luki"     = paste0("Luki (n = ", true_luki_n, ")"),
    "Yangambi" = paste0("Yangambi (n = ", true_yangambi_n, ")")
  )
  
  # HARDCODED CLEAN NUMERIC STRINGS: Bypasses float evaluation bugs completely
  if(target_metric == "Duration") {
    x_limits <- c(-1.6, 1.6)
    x_breaks <- c(-1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2)
    x_labels <- c("-1.2", "-0.8", "-0.4", "0.0", "0.4", "0.8", "1.2")
  } else {
    x_limits <- c(-1.0, 1.0)
    x_breaks <- c(-0.8, -0.4, 0.0, 0.4, 0.8)
    x_labels <- c("-0.8", "-0.4", "0.0", "0.4", "0.8")
  }
  
  forest_p <- ggplot(plot_data, aes(x = Estimate, y = Term, color = Site, fill = Site)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    
    geom_errorbarh(aes(xmin = Estimate - 1.96 * StdError, xmax = Estimate + 1.96 * StdError), 
                   height = 0.12, linewidth = 0.8, position = position_dodge(width = 0.3)) +
    
    geom_point(size = 3.5, shape = 21, color = "black", position = position_dodge(width = 0.3)) +
    
    geom_text(aes(x = Label_Pos, label = Sig_Label),
              position = position_dodge(width = 0.3), 
              size = 5.5, fontface = "bold", color = "darkred", show.legend = FALSE, hjust = 0, vjust = 0.3) +
    
    scale_color_manual(values = c("Luki" = "#3A93D1", "Yangambi" = "#F4D03F"), labels = legend_labels) +
    scale_fill_manual(values = c("Luki" = "#3A93D1", "Yangambi" = "#F4D03F"), labels = legend_labels) +
    
    # Overwrites automated limits using clean explicit string overrides
    scale_x_continuous(limits = x_limits, breaks = x_breaks, labels = x_labels) +
    labs(
      title = paste0(target_metric, " Dynamics (Yangambi vs. Luki)"), 
      x = expression(bold("Standardised Effect Size ("*beta*" ± 95% CI)")), 
      y = NULL
    ) +
    
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.15, margin = margin(b = 15)),
      axis.text.x = element_text(color = "black", size = 10),
      axis.text.y = element_text(color = "black", size = 10),
      axis.title.x = element_text(face = "bold", size = 11, margin = margin(t = 10)),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 10),
      plot.margin = margin(t = 15, r = 50, b = 15, l = 15, unit = "pt")
    )
  
  pdf_filename <- paste0(output_dir, "figure_", tolower(target_metric), "_forest_plot.pdf")
  pdf(file = pdf_filename, width = 7.5, height = 5, useDingbats = FALSE)
  print(forest_p)
  dev.off()
  
  message("Saved tailored vector forest plot to: ", pdf_filename)
}


# EXECUTE CODE GENERATIONS

build_forest_plot("Onset")
build_forest_plot("End")
build_forest_plot("Duration")


##--------------------------------------------------------------------
##PART 2. Temporal trends


# Flowering onset – temporal trend

flower_luki_onset_year <- pglmm(
  onset_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = luki_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)

summary(flower_luki_onset_year)
# > flower_luki_onset_year
# Linear mixed model fit by restricted maximum likelihood
# 
# Call:onset_doy ~ year_c
# 
# logLik    AIC    BIC 
#  -1019   2054   2087 
# 
# Random effects:
#                       Variance  Std.Dev
# 1|species            4.319e+00 2.078178
# 1|species__          4.865e-02 0.220576
# 0 + year_c|species   8.438e-03 0.091856
# 0 + year_c|species__ 2.428e-06 0.001558
# 1|year               2.115e-01 0.459865
# residual             9.986e-02 0.316004
# 
# Fixed effects:
#                 Value Std.Error  Zscore Pvalue    
# (Intercept) 41.170607  0.709215 58.0510 <2e-16 ***
# year_c      -0.026296  0.044620 -0.5893 0.5556    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

flower_luki_end_year <- pglmm(
  end_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = luki_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)

summary(flower_luki_end_year)
# > flower_luki_end_year
# Linear mixed model fit by restricted maximum likelihood
# 
# Call:end_doy ~ year_c
# 
# logLik    AIC    BIC 
#  -1019   2053   2087 
# 
# Random effects:
#                       Variance   Std.Dev
# 1|species            4.664e+00 2.1597103
# 1|species__          3.717e-02 0.1927935
# 0 + year_c|species   9.113e-03 0.0954634
# 0 + year_c|species__ 1.478e-07 0.0003844
# 1|year               2.490e-01 0.4989701
# residual             9.823e-02 0.3134091
# 
# Fixed effects:
#                  Value  Std.Error   Zscore Pvalue    
# (Intercept) 333.342938   0.762425 437.2139 <2e-16 ***
# year_c        0.027096   0.048328   0.5607  0.575    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

flower_luki_duration_year <- pglmm(
  duration_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = luki_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)
summary(flower_luki_duration_year)

# > flower_luki_duration_year
# Linear mixed model fit by restricted maximum likelihood
# 
# Call:duration_doy ~ year_c
# 
# logLik    AIC    BIC 
#  -2036   4087   4121 
# 
# Random effects:
#                       Variance  Std.Dev
# 1|species            1.792e+01 4.233731
# 1|species__          1.393e-01 0.373271
# 0 + year_c|species   3.514e-02 0.187446
# 0 + year_c|species__ 7.548e-06 0.002747
# 1|year               3.514e-03 0.059282
# residual             3.917e-01 0.625854
# 
# Fixed effects:
#                  Value  Std.Error   Zscore    Pvalue    
# (Intercept) 292.172824   0.443006 659.5236 < 2.2e-16 ***
# year_c        0.052948   0.017723   2.9875  0.002812 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##------------------------------------------------

# YANGAMBI: TEMPORAL TREND MODELS

# 1. Flowering onset – temporal trend

flower_yangambi_onset_year <- pglmm(
  onset_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = yangambi_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)

summary(flower_yangambi_onset_year)
# > flower_yangambi_onset_year
# Linear mixed model fit by restricted maximum likelihood
# 
# Call:onset_doy ~ year_c
# 
# logLik    AIC    BIC 
#  -8926  17867  17912 
# 
# Random effects:
#                       Variance Std.Dev
# 1|species            0.9828777 0.99140
# 1|species__          0.0226064 0.15035
# 0 + year_c|species   0.0043860 0.06623
# 0 + year_c|species__ 0.0001675 0.01294
# 1|year               0.3696757 0.60801
# residual             0.7875269 0.88743
# 
# Fixed effects:
#                 Value Std.Error   Zscore Pvalue    
# (Intercept) 42.147336  0.284777 148.0009 <2e-16 ***
# year_c       0.036008  0.024468   1.4716 0.1411    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# 2. Flowering end – temporal trend


flower_yangambi_end_year <- pglmm(
  end_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = yangambi_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)

summary(flower_yangambi_end_year)

# > flower_yangambi_end_year
# Linear mixed model fit by restricted maximum likelihood
# 
# Call:end_doy ~ year_c
# 
# logLik    AIC    BIC 
#  -9737  19490  19535 
# 
# Random effects:
#                       Variance Std.Dev
# 1|species            1.2855215 1.13381
# 1|species__          0.0281063 0.16765
# 0 + year_c|species   0.0057203 0.07563
# 0 + year_c|species__ 0.0002288 0.01513
# 1|year               0.6774482 0.82307
# residual             1.0285490 1.01417
# 
# Fixed effects:
#                  Value  Std.Error   Zscore Pvalue    
# (Intercept) 329.438178   0.373889 881.1131 <2e-16 ***
# year_c       -0.035244   0.032247  -1.0929 0.2744    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# 3. Flowering duration – temporal trend

flower_yangambi_duration_year <- pglmm(
  duration_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = yangambi_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)


summary(flower_yangambi_duration_year)
# > flower_yangambi_duration_year
# Linear mixed model fit by restricted maximum likelihood
# 
# Call:duration_doy ~ year_c
# 
# logLik    AIC    BIC 
# -13564  27145  27189 
# 
# Random effects:
#                       Variance Std.Dev
# 1|species            4.469e+00 2.11403
# 1|species__          1.217e-01 0.34880
# 0 + year_c|species   2.036e-02 0.14269
# 0 + year_c|species__ 7.819e-04 0.02796
# 1|year               6.927e+01 8.32289
# residual             3.600e+00 1.89748
# 
# Fixed effects:
#                  Value  Std.Error  Zscore Pvalue    
# (Intercept) 287.302105   3.517729 81.6726 <2e-16 ***
# year_c       -0.075728   0.300860 -0.2517 0.8013    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Figure 2. Temporal trends in flowering phenology across both sites. 

library(phyr)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)

output_dir <- "E:/PAPER 2/PGLMM"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

luki_data <- flowering %>%
  filter(site == "Luki")

yangambi_data <- flowering %>%
  filter(site == "Yangambi")

flower_luki_onset_year <- pglmm(
  onset_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = luki_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)

flower_luki_end_year <- pglmm(
  end_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = luki_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)

flower_luki_duration_year <- pglmm(
  duration_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = luki_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)

flower_yangambi_onset_year <- pglmm(
  onset_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = yangambi_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)

flower_yangambi_end_year <- pglmm(
  end_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = yangambi_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)

flower_yangambi_duration_year <- pglmm(
  duration_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = yangambi_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)

summary(flower_luki_onset_year)
summary(flower_luki_end_year)
summary(flower_luki_duration_year)

summary(flower_yangambi_onset_year)
summary(flower_yangambi_end_year)
summary(flower_yangambi_duration_year)

combined_points <- bind_rows(
  luki_data %>%
    mutate(site = "Luki"),
  yangambi_data %>%
    mutate(site = "Yangambi")
) %>%
  mutate(
    year = as.numeric(as.character(year)),
    year_c = as.numeric(year_c)
  )

year_lookup <- combined_points %>%
  select(year, year_c) %>%
  distinct() %>%
  arrange(year)

model_parameters <- data.frame(
  site = c(
    "Luki",
    "Yangambi"
  ),
  onset_intercept = c(
    41.170607,
    42.147336
  ),
  onset_beta = c(
    -0.026296,
    0.036008
  ),
  onset_p = c(
    0.5556,
    0.1411
  ),
  end_intercept = c(
    333.342938,
    329.438178
  ),
  end_beta = c(
    0.027096,
    -0.035244
  ),
  end_p = c(
    0.5750,
    0.2744
  ),
  duration_intercept = c(
    292.172824,
    287.302105
  ),
  duration_beta = c(
    0.052948,
    -0.075728
  ),
  duration_p = c(
    0.002812,
    0.8013
  )
)

pred_grid <- expand.grid(
  year = year_lookup$year,
  site = c("Yangambi", "Luki"),
  stringsAsFactors = FALSE
) %>%
  left_join(
    year_lookup,
    by = "year"
  ) %>%
  left_join(
    model_parameters,
    by = "site"
  ) %>%
  mutate(
    onset_fit =
      onset_intercept +
      onset_beta * year_c,
    end_fit =
      end_intercept +
      end_beta * year_c,
    duration_fit =
      duration_intercept +
      duration_beta * year_c
  )

site_colors <- c(
  "Luki" = "#D55E00",
  "Yangambi" = "#0072B2"
)

theme_publication <- function() {
  
  theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 8)
      ),
      axis.title = element_text(
        size = 10.5,
        face = "bold"
      ),
      axis.text = element_text(
        size = 9,
        color = "black"
      ),
      panel.grid.major.x = element_line(
        color = "grey90",
        linewidth = 0.35
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(
        size = 9
      ),
      legend.position = "bottom",
      plot.margin = margin(
        10,
        10,
        10,
        10
      )
    )
}

p_A <- ggplot() +
  
  geom_point(
    data = combined_points,
    aes(
      x = year,
      y = onset_doy,
      color = site
    ),
    alpha = 0.18,
    size = 1.5
  ) +
  
  geom_line(
    data = pred_grid,
    aes(
      x = year,
      y = onset_fit,
      color = site
    ),
    linewidth = 1.1
  ) +
  
  scale_color_manual(
    values = site_colors
  ) +
  
  scale_x_continuous(
    breaks = c(
      1940,
      1945,
      1950,
      1955
    ),
    expand = expansion(
      mult = c(
        0.01,
        0.02
      )
    )
  ) +
  
  labs(
    title = "A. Flowering onset",
    x = "Calendar year",
    y = "Flowering onset (day of year)"
  ) +
  
  annotate(
    "text",
    x = 1937.5,
    y = 46,
    label = "Yangambi: p = 0.1411\nLuki: p = 0.5556",
    hjust = 0,
    vjust = 1,
    size = 3.1,
    color = "grey25"
  ) +
  
  theme_publication()

p_B <- ggplot() +
  
  geom_point(
    data = combined_points,
    aes(
      x = year,
      y = end_doy,
      color = site
    ),
    alpha = 0.18,
    size = 1.5
  ) +
  
  geom_line(
    data = pred_grid,
    aes(
      x = year,
      y = end_fit,
      color = site
    ),
    linewidth = 1.1
  ) +
  
  scale_color_manual(
    values = site_colors
  ) +
  
  scale_x_continuous(
    breaks = c(
      1940,
      1945,
      1950,
      1955
    ),
    expand = expansion(
      mult = c(
        0.01,
        0.02
      )
    )
  ) +
  
  labs(
    title = "B. Flowering termination",
    x = "Calendar year",
    y = "Flowering termination (day of year)"
  ) +
  
  annotate(
    "text",
    x = 1947,
    y = 337,
    label = "Yangambi: p = 0.2744\nLuki: p = 0.5750",
    hjust = 0,
    vjust = 1,
    size = 3.1,
    color = "grey25"
  ) +
  
  theme_publication()

p_C <- ggplot() +
  
  geom_point(
    data = combined_points,
    aes(
      x = year,
      y = duration_doy,
      color = site
    ),
    alpha = 0.18,
    size = 1.5
  ) +
  
  geom_line(
    data = pred_grid,
    aes(
      x = year,
      y = duration_fit,
      color = site
    ),
    linewidth = 1.1
  ) +
  
  scale_color_manual(
    values = site_colors
  ) +
  
  scale_x_continuous(
    breaks = c(
      1940,
      1945,
      1950,
      1955
    ),
    expand = expansion(
      mult = c(
        0.01,
        0.02
      )
    )
  ) +
  
  labs(
    title = "C. Flowering duration",
    x = "Calendar year",
    y = "Flowering duration (days)"
  ) +
  
  annotate(
    "text",
    x = 1947,
    y = 300,
    label = "Yangambi: p = 0.8013\nLuki: p = 0.0028 **",
    hjust = 0,
    vjust = 1,
    size = 3.1,
    color = "grey25"
  ) +
  
  theme_publication()

p_A
p_B
p_C

trend_multi_panel <- grid.arrange(
  p_A,
  p_B,
  p_C,
  ncol = 3
)

trend_multi_panel

target_file <- file.path(
  output_dir,
  "Figure2_Phenological_Temporal_Trends.pdf"
)

pdf(
  target_file,
  width = 14,
  height = 5,
  useDingbats = FALSE
)

grid.arrange(
  p_A,
  p_B,
  p_C,
  ncol = 3
)

dev.off()

png_file <- file.path(
  output_dir,
  "Figure2_Phenological_Temporal_Trends.png"
)

ggsave(
  filename = png_file,
  plot = trend_multi_panel,
  width = 14,
  height = 5,
  units = "in",
  dpi = 600
)

print(target_file)
print(png_file)


##FIGURS 4-9

#Plots of  FLOWERING TEMPORAL AND FAMILY-SPECIFIC ANALYSIS



# BLOCK 1 — LOAD REQUIRED PACKAGES


library(phyr)
library(dplyr)
library(ggplot2)
library(openxlsx)



# BLOCK 2 — DEFINE OUTPUT FOLDER


output_folder <- "E:/PAPER 2/PGLMM"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}



# BLOCK 3 — PREPARE SITE-SPECIFIC DATA


# We analyse Luki and Yangambi separately.
#
# year_c is centred separately within each site:
#
#
# This makes zero correspond to the mean observation year at each site
# and avoids using a pooled temporal centre across the two sites.

luki_data <- flowering %>%
  filter(site == "Luki") %>%
  mutate(
    year_c = year - mean(year, na.rm = TRUE)
  )

yangambi_data <- flowering %>%
  filter(site == "Yangambi") %>%
  mutate(
    year_c = year - mean(year, na.rm = TRUE)
  )


# Check centring
cat("\nMean centred year — Luki:\n")
print(mean(luki_data$year_c, na.rm = TRUE))

cat("\nMean centred year — Yangambi:\n")
print(mean(yangambi_data$year_c, na.rm = TRUE))



# BLOCK 4 — PREPARE YANGAMBI PHYLOGENY


# Keep only species occurring at Yangambi.

yangambi_species <- unique(yangambi_data$species)

phy_yang <- drop.tip(
  phy,
  setdiff(phy$tip.label, yangambi_species)
)

cat("\nNumber of species in Yangambi data:",
    length(unique(yangambi_data$species)), "\n")

cat("Number of tips in Yangambi phylogeny:",
    length(phy_yang$tip.label), "\n")



# BLOCK 5 — FIT LUKI ONSET TEMPORAL MODEL


# Response:
# onset_doy
#
# Fixed effect:
# year_c
#
# Random effects:
# (1 | species__)              = species intercept with phylogenetic structure
# (0 + year_c | species__)     = phylogenetically structured species slopes
# (1 | year)                   = shared year-to-year variation

flower_luki_onset_year <- pglmm(
  onset_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = luki_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)



# BLOCK 6 — FIT LUKI END TEMPORAL MODEL


flower_luki_end_year <- pglmm(
  end_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = luki_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)



# BLOCK 7 — FIT LUKI DURATION TEMPORAL MODEL


flower_luki_duration_year <- pglmm(
  duration_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = luki_data,
  family = "gaussian",
  cov_ranef = list(species = phy)
)



# BLOCK 8 — FIT YANGAMBI ONSET TEMPORAL MODEL


flower_yangambi_onset_year <- pglmm(
  onset_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = yangambi_data,
  family = "gaussian",
  cov_ranef = list(species = phy_yang)
)



# BLOCK 9 — FIT YANGAMBI END TEMPORAL MODEL


flower_yangambi_end_year <- pglmm(
  end_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = yangambi_data,
  family = "gaussian",
  cov_ranef = list(species = phy_yang)
)



# BLOCK 10 — FIT YANGAMBI DURATION TEMPORAL MODEL


flower_yangambi_duration_year <- pglmm(
  duration_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = yangambi_data,
  family = "gaussian",
  cov_ranef = list(species = phy_yang)
)



# BLOCK 11 — CHECK MODEL SUMMARIES


summary(flower_luki_onset_year)
summary(flower_luki_end_year)
summary(flower_luki_duration_year)

summary(flower_yangambi_onset_year)
summary(flower_yangambi_end_year)
summary(flower_yangambi_duration_year)



# BLOCK 12 — EXTRACT COMMUNITY-LEVEL TEMPORAL EFFECTS


# The fixed year_c coefficient is the primary community-level temporal test.
#
# We use matrix position [2,1] because B.se, B.zscore and B.pvalue
# may not have useful row names in communityPGLMM objects.

extract_pglmm_result <- function(
    model,
    site_name,
    response_name
) {
  
  data.frame(
    site = site_name,
    response = response_name,
    beta = as.numeric(model$B[2, 1]),
    SE = as.numeric(model$B.se[2, 1]),
    Z = as.numeric(model$B.zscore[2, 1]),
    P = as.numeric(model$B.pvalue[2, 1])
  )
}


pglmm_temporal_results <- bind_rows(
  
  extract_pglmm_result(
    flower_luki_onset_year,
    "Luki",
    "Onset"
  ),
  
  extract_pglmm_result(
    flower_luki_end_year,
    "Luki",
    "End"
  ),
  
  extract_pglmm_result(
    flower_luki_duration_year,
    "Luki",
    "Duration"
  ),
  
  extract_pglmm_result(
    flower_yangambi_onset_year,
    "Yangambi",
    "Onset"
  ),
  
  extract_pglmm_result(
    flower_yangambi_end_year,
    "Yangambi",
    "End"
  ),
  
  extract_pglmm_result(
    flower_yangambi_duration_year,
    "Yangambi",
    "Duration"
  )
)


print(pglmm_temporal_results)



# BLOCK 13 — FUNCTION TO EXTRACT SPECIES-SPECIFIC TEMPORAL SLOPES


# IMPORTANT:
#
# ranef(model) does not directly provide the complete species-specific
# temporal slope table required here.
#
# The fitted communityPGLMM object contains:
#
# St  = random-effect block structure
# Zt  = random-effect design matrix
# s2r = fitted random-effect variances
# iV  = inverse covariance structure
# H   = residual component
# B   = fixed effects
#
# The total species slope is:
#
# fixed year slope
# +
# identity species slope
# +
# phylogenetic species slope


extract_species_slopes <- function(
    model,
    site_name,
    response_name
) {
  
  
  # STEP 1 — Identify random-effect block sizes
  
  
  block_sizes <- apply(
    model$St != 0,
    1,
    sum
  )
  
  
  
  # STEP 2 — Find beginning and ending positions of each random-effect block
  
  
  block_start <- c(
    1,
    cumsum(block_sizes)[-length(block_sizes)] + 1
  )
  
  block_end <- cumsum(block_sizes)
  
  
  
  # STEP 3 — Obtain unique species names
  
  
  # IMPORTANT:
  # We use levels(), not as.character().
  #
  # as.character() would return one value per observation.
  # levels() returns one value per species.
  
  species_factor <- model$random.effects[[3]]$species
  
  species_names <- levels(species_factor)
  
  n_species <- length(species_names)
  
  
  
  # STEP 4 — Identify the two species-slope variance components
  
  
  variance_names <- names(model$s2r)
  
  
  identity_index <- grep(
    "0 \\+ year_c\\|species$",
    variance_names
  )
  
  
  phylogenetic_index <- grep(
    "0 \\+ year_c\\|species__",
    variance_names
  )
  
  
  identity_variance <- as.numeric(
    model$s2r[identity_index]
  )
  
  
  phylogenetic_variance <- as.numeric(
    model$s2r[phylogenetic_index]
  )
  
  
  
  # STEP 5 — Extract the fixed temporal slope
  
  
  fixed_slope <- as.numeric(
    model$B["year_c", 1]
  )
  
  
  
  # STEP 6 — Extract identity species-slope block
  
  
  Zt_identity <- model$Zt[
    block_start[3]:block_end[3],
    ,
    drop = FALSE
  ]
  
  
  
  # STEP 7 — Extract phylogenetic species-slope block
  
  
  Zt_phylogenetic <- model$Zt[
    block_start[4]:block_end[4],
    ,
    drop = FALSE
  ]
  
  
  
  # STEP 8 — Calculate identity species slope
  
  
  identity_slope <- identity_variance *
    as.numeric(
      Zt_identity %*%
        model$iV %*%
        model$H
    )
  
  
  
  # STEP 9 — Extract fitted phylogenetic covariance matrix
  
  
  phylo_covariance <-
    model$random.effects[[4]]$covar
  
  
  
  # STEP 10 — Calculate phylogenetic species slope
  
  
  phylogenetic_slope <-
    phylogenetic_variance *
    as.numeric(
      t(chol(phylo_covariance)) %*%
        Zt_phylogenetic %*%
        model$iV %*%
        model$H
    )
  
  
  
  # STEP 11 — Calculate total species-specific temporal slope
  
  
  total_slope <-
    fixed_slope +
    identity_slope +
    phylogenetic_slope
  
  
  
  # STEP 12 — Translate slope sign into biological direction
  
  
  direction <- case_when(
    
    response_name == "Onset" &
      total_slope < 0 ~
      "Advance",
    
    response_name == "Onset" &
      total_slope > 0 ~
      "Delay",
    
    response_name == "End" &
      total_slope < 0 ~
      "Earlier termination",
    
    response_name == "End" &
      total_slope > 0 ~
      "Later termination",
    
    response_name == "Duration" &
      total_slope < 0 ~
      "Decrease",
    
    response_name == "Duration" &
      total_slope > 0 ~
      "Increase",
    
    TRUE ~
      "No change"
  )
  
  
  
  # STEP 13 — Return species-level table
  
  
  data.frame(
    
    site = site_name,
    
    response = response_name,
    
    species = species_names,
    
    fixed_slope = fixed_slope,
    
    identity_slope = identity_slope,
    
    phylogenetic_slope = phylogenetic_slope,
    
    total_slope = total_slope,
    
    direction = direction,
    
    stringsAsFactors = FALSE
  )
}



# BLOCK 14 — EXTRACT LUKI SPECIES SLOPES


luki_onset <- extract_species_slopes(
  flower_luki_onset_year,
  "Luki",
  "Onset"
)

luki_end <- extract_species_slopes(
  flower_luki_end_year,
  "Luki",
  "End"
)

luki_duration <- extract_species_slopes(
  flower_luki_duration_year,
  "Luki",
  "Duration"
)



# BLOCK 15 — EXTRACT YANGAMBI SPECIES SLOPES


yangambi_onset <- extract_species_slopes(
  flower_yangambi_onset_year,
  "Yangambi",
  "Onset"
)

yangambi_end <- extract_species_slopes(
  flower_yangambi_end_year,
  "Yangambi",
  "End"
)

yangambi_duration <- extract_species_slopes(
  flower_yangambi_duration_year,
  "Yangambi",
  "Duration"
)



# BLOCK 16 — CHECK SPECIES COUNTS


cat("\nLuki onset species:",
    nrow(luki_onset), "\n")

cat("Luki end species:",
    nrow(luki_end), "\n")

cat("Luki duration species:",
    nrow(luki_duration), "\n")

cat("\nYangambi onset species:",
    nrow(yangambi_onset), "\n")

cat("Yangambi end species:",
    nrow(yangambi_end), "\n")

cat("Yangambi duration species:",
    nrow(yangambi_duration), "\n")



# BLOCK 17 — COMBINE ALL SPECIES RESULTS


all_species_results <- bind_rows(
  
  luki_onset,
  luki_end,
  luki_duration,
  
  yangambi_onset,
  yangambi_end,
  yangambi_duration
)



# BLOCK 18 — SPECIES-LEVEL TEMPORAL SUMMARY


species_temporal_summary <- all_species_results %>%
  
  group_by(
    site,
    response
  ) %>%
  
  summarise(
    
    n_species =
      n_distinct(species),
    
    mean_slope =
      mean(
        total_slope,
        na.rm = TRUE
      ),
    
    median_slope =
      median(
        total_slope,
        na.rm = TRUE
      ),
    
    Q1 =
      as.numeric(
        quantile(
          total_slope,
          0.25,
          na.rm = TRUE
        )
      ),
    
    Q3 =
      as.numeric(
        quantile(
          total_slope,
          0.75,
          na.rm = TRUE
        )
      ),
    
    IQR =
      IQR(
        total_slope,
        na.rm = TRUE
      ),
    
    SD =
      sd(
        total_slope,
        na.rm = TRUE
      ),
    
    negative_n =
      sum(
        total_slope < 0,
        na.rm = TRUE
      ),
    
    positive_n =
      sum(
        total_slope > 0,
        na.rm = TRUE
      ),
    
    negative_percent =
      100 * negative_n / n_species,
    
    positive_percent =
      100 * positive_n / n_species,
    
    .groups = "drop"
  )


print(species_temporal_summary)



# BLOCK 19 — SPECIES DIRECTION COUNTS


species_direction_counts <- all_species_results %>%
  
  group_by(
    site,
    response,
    direction
  ) %>%
  
  summarise(
    n_species = n_distinct(species),
    .groups = "drop"
  ) %>%
  
  group_by(
    site,
    response
  ) %>%
  
  mutate(
    percent =
      100 * n_species /
      sum(n_species)
  ) %>%
  
  ungroup()


print(species_direction_counts)



# BLOCK 20 — WILCOXON SIGNED-RANK TESTS


# This is a secondary analysis.
#
# It tests whether the distribution of estimated species slopes is shifted
# from zero.
#
# It does NOT test whether each individual species has a significant slope.

wilcoxon_species_slopes <- all_species_results %>%
  
  group_by(
    site,
    response
  ) %>%
  
  summarise(
    
    n_species = n(),
    
    W =
      unname(
        wilcox.test(
          total_slope,
          mu = 0,
          exact = FALSE
        )$statistic
      ),
    
    P =
      wilcox.test(
        total_slope,
        mu = 0,
        exact = FALSE
      )$p.value,
    
    median_slope =
      median(
        total_slope,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    
    P_formatted =
      ifelse(
        P < 0.001,
        "<0.001",
        sprintf(
          "%.4f",
          P
        )
      )
  )


print(wilcoxon_species_slopes)



# BLOCK 21 — CREATE SPECIES × FAMILY LOOKUP TABLES


luki_species_family <- luki_data %>%
  
  select(
    species,
    family
  ) %>%
  
  distinct()


yangambi_species_family <- yangambi_data %>%
  
  select(
    species,
    family
  ) %>%
  
  distinct()



# BLOCK 22 — CHECK THAT SPECIES HAVE ONE FAMILY


luki_family_check <- luki_species_family %>%
  
  group_by(species) %>%
  
  summarise(
    n_family =
      n_distinct(family),
    .groups = "drop"
  ) %>%
  
  filter(
    n_family > 1
  )


yangambi_family_check <- yangambi_species_family %>%
  
  group_by(species) %>%
  
  summarise(
    n_family =
      n_distinct(family),
    .groups = "drop"
  ) %>%
  
  filter(
    n_family > 1
  )


print(luki_family_check)
print(yangambi_family_check)



# BLOCK 23 — JOIN FAMILY INFORMATION TO SPECIES SLOPES


luki_onset_family <- luki_onset %>%
  
  left_join(
    luki_species_family,
    by = "species"
  )


luki_end_family <- luki_end %>%
  
  left_join(
    luki_species_family,
    by = "species"
  )


luki_duration_family <- luki_duration %>%
  
  left_join(
    luki_species_family,
    by = "species"
  )


yangambi_onset_family <- yangambi_onset %>%
  
  left_join(
    yangambi_species_family,
    by = "species"
  )


yangambi_end_family <- yangambi_end %>%
  
  left_join(
    yangambi_species_family,
    by = "species"
  )


yangambi_duration_family <- yangambi_duration %>%
  
  left_join(
    yangambi_species_family,
    by = "species"
  )



# BLOCK 24 — FAMILY SUMMARY FUNCTION


calculate_family_distribution <- function(
    species_data,
    site_name,
    response_name
) {
  
  required_columns <- c(
    "species",
    "family",
    "total_slope"
  )
  
  
  missing_columns <- setdiff(
    required_columns,
    names(species_data)
  )
  
  
  if (length(missing_columns) > 0) {
    
    stop(
      paste(
        "Missing columns:",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )
  }
  
  
  result <- species_data %>%
    
    filter(
      !is.na(family),
      family != "",
      !is.na(total_slope)
    ) %>%
    
    group_by(family) %>%
    
    summarise(
      
      n_species =
        n_distinct(species),
      
      mean_slope =
        mean(
          total_slope,
          na.rm = TRUE
        ),
      
      median_slope =
        median(
          total_slope,
          na.rm = TRUE
        ),
      
      Q1 =
        as.numeric(
          quantile(
            total_slope,
            0.25,
            na.rm = TRUE
          )
        ),
      
      Q3 =
        as.numeric(
          quantile(
            total_slope,
            0.75,
            na.rm = TRUE
          )
        ),
      
      IQR =
        IQR(
          total_slope,
          na.rm = TRUE
        ),
      
      SD =
        sd(
          total_slope,
          na.rm = TRUE
        ),
      
      min_slope =
        min(
          total_slope,
          na.rm = TRUE
        ),
      
      max_slope =
        max(
          total_slope,
          na.rm = TRUE
        ),
      
      negative_n =
        sum(
          total_slope < 0,
          na.rm = TRUE
        ),
      
      positive_n =
        sum(
          total_slope > 0,
          na.rm = TRUE
        ),
      
      zero_n =
        sum(
          total_slope == 0,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    ) %>%
    
    mutate(
      
      site =
        site_name,
      
      response =
        response_name,
      
      negative_percent =
        100 *
        negative_n /
        n_species,
      
      positive_percent =
        100 *
        positive_n /
        n_species,
      
      dominant_direction =
        case_when(
          
          negative_n >
            positive_n ~
            "Negative",
          
          positive_n >
            negative_n ~
            "Positive",
          
          TRUE ~
            "Equal"
        )
    ) %>%
    
    select(
      
      site,
      response,
      family,
      n_species,
      
      mean_slope,
      median_slope,
      
      Q1,
      Q3,
      IQR,
      SD,
      
      min_slope,
      max_slope,
      
      negative_n,
      positive_n,
      zero_n,
      
      negative_percent,
      positive_percent,
      
      dominant_direction
    )
  
  
  return(result)
}



# BLOCK 25 — CALCULATE ALL FAMILY SUMMARIES


family_luki_onset <-
  calculate_family_distribution(
    luki_onset_family,
    "Luki",
    "Onset"
  )


family_luki_end <-
  calculate_family_distribution(
    luki_end_family,
    "Luki",
    "End"
  )


family_luki_duration <-
  calculate_family_distribution(
    luki_duration_family,
    "Luki",
    "Duration"
  )


family_yangambi_onset <-
  calculate_family_distribution(
    yangambi_onset_family,
    "Yangambi",
    "Onset"
  )


family_yangambi_end <-
  calculate_family_distribution(
    yangambi_end_family,
    "Yangambi",
    "End"
  )


family_yangambi_duration <-
  calculate_family_distribution(
    yangambi_duration_family,
    "Yangambi",
    "Duration"
  )



# BLOCK 26 — COMBINE ALL FAMILY RESULTS


family_temporal_statistics <- bind_rows(
  
  family_luki_onset,
  family_luki_end,
  family_luki_duration,
  
  family_yangambi_onset,
  family_yangambi_end,
  family_yangambi_duration
)



# BLOCK 27 — CHECK FAMILY AND SPECIES COUNTS


family_counts <- family_temporal_statistics %>%
  
  group_by(
    site,
    response
  ) %>%
  
  summarise(
    
    n_families =
      n_distinct(family),
    
    n_species =
      sum(n_species),
    
    .groups = "drop"
  )


print(family_counts)



# Luki      = 140 species across 38 families
# Yangambi  = 595 species across 81 families



# BLOCK 28 — FINAL FAMILY FOREST-PLOT FUNCTION


plot_family_forest <- function(
    data,
    site_name,
    response_name
) {
  
  
  # Select relevant site and response
  
  
  plot_data <- data %>%
    
    filter(
      site == site_name,
      response == response_name
    ) %>%
    
    filter(
      !is.na(family),
      family != "",
      !is.na(mean_slope),
      !is.na(Q1),
      !is.na(Q3)
    )
  
  
  if (nrow(plot_data) == 0) {
    
    stop(
      paste(
        "No data available for",
        site_name,
        response_name
      )
    )
  }
  
  
  
  # Number of families and species
  
  
  n_families <-
    n_distinct(plot_data$family)
  
  n_species <-
    sum(
      plot_data$n_species,
      na.rm = TRUE
    )
  
  
  
  # Response-specific biological terminology
  
  
  if (response_name == "Onset") {
    
    negative_label <-
      "Advanced trend"
    
    positive_label <-
      "Delayed trend"
    
  } else if (response_name == "End") {
    
    negative_label <-
      "Earlier termination"
    
    positive_label <-
      "Later termination"
    
  } else if (response_name == "Duration") {
    
    negative_label <-
      "Shorter duration"
    
    positive_label <-
      "Longer duration"
    
  } else {
    
    stop(
      "response_name must be 'Onset', 'End', or 'Duration'"
    )
  }
  
  
  
  # Classify families according to mean slope
  
  
  plot_data <- plot_data %>%
    
    mutate(
      
      phenological_shift =
        case_when(
          
          mean_slope < 0 ~
            negative_label,
          
          mean_slope > 0 ~
            positive_label,
          
          TRUE ~
            NA_character_
        )
    ) %>%
    
    arrange(mean_slope)
  
  
  
  # negative / advanced
  #        ↓
  # positive / delayed
  
  plot_data$Plant_Family <- factor(
    plot_data$family,
    levels = plot_data$family
  )
  
  
  
  # X-axis limits
  
  
  all_x <- c(
    plot_data$Q1,
    plot_data$Q3,
    plot_data$mean_slope
  )
  
  
  x_range <-
    range(
      all_x,
      na.rm = TRUE
    )
  
  
  x_padding <-
    diff(x_range) * 0.15
  
  
  if (
    !is.finite(x_padding) ||
    x_padding == 0
  ) {
    
    x_padding <- 0.05
  }
  
  
  x_min <-
    x_range[1] -
    x_padding
  
  
  x_max <-
    x_range[2] +
    x_padding
  
  
  
  # Construct plot
  
  
  p <- ggplot(
    
    plot_data,
    
    aes(
      x = mean_slope,
      y = Plant_Family,
      color = phenological_shift
    )
    
  ) +
    
    
    # Zero reference line
    
    geom_vline(
      xintercept = 0,
      color = "grey45",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    
    
    # Interquartile range
    
    geom_segment(
      
      aes(
        x = Q1,
        xend = Q3,
        y = Plant_Family,
        yend = Plant_Family
      ),
      
      linewidth = 0.8
    ) +
    
    
    # Family mean
    
    geom_point(
      size = 2.8
    ) +
    
    
    # Phenological shift legend
    
    scale_color_manual(
      
      name =
        "Phenological shifts",
      
      values = c(
        
        setNames(
          "#2C7FB8",
          negative_label
        ),
        
        setNames(
          "#E66101",
          positive_label
        )
      )
    ) +
    
    
    # X-axis
    
    scale_x_continuous(
      
      limits =
        c(
          x_min,
          x_max
        ),
      
      breaks =
        pretty(
          c(
            x_min,
            x_max
          ),
          n = 5
        ),
      
      expand =
        c(
          0,
          0
        )
    ) +
    
    
    # Labels
    
    labs(
      
      x =
        expression(
          "Mean family-specific temporal slope (days year"^{-1}*")"
        ),
      
      y =
        "Plant family",
      
      title =
        paste(
          site_name,
          "— Flowering",
          tolower(response_name),
          "temporal trends across families"
        )
    ) +
    
    
    # Space above plot for species/family count
    
    coord_cartesian(
      
      ylim =
        c(
          0.5,
          n_families + 2
        ),
      
      clip =
        "off"
    ) +
    
    
    # Species/family information
    
    annotate(
      
      "text",
      
      x =
        x_min,
      
      y =
        n_families + 1.25,
      
      label =
        paste0(
          n_species,
          " species across ",
          n_families,
          " families"
        ),
      
      hjust =
        0,
      
      vjust =
        0.5,
      
      size =
        3.4,
      
      fontface =
        "bold",
      
      color =
        "grey25"
    ) +
    
    
    # Theme
    
    theme_bw(
      base_size = 11
    ) +
    
    
    theme(
      
      plot.title =
        element_text(
          face = "bold",
          size = 12,
          hjust = 0.5,
          margin =
            margin(
              b = 14
            )
        ),
      
      axis.title.x =
        element_text(
          size = 11,
          margin =
            margin(
              t = 10
            )
        ),
      
      axis.title.y =
        element_text(
          size = 11,
          margin =
            margin(
              r = 10
            )
        ),
      
      axis.text.x =
        element_text(
          color = "black",
          size = 10
        ),
      
      axis.text.y =
        element_text(
          color = "black",
          size =
            ifelse(
              n_families > 70,
              7,
              8.5
            ),
          face = "italic"
        ),
      
      panel.grid.major.y =
        element_blank(),
      
      panel.grid.minor =
        element_blank(),
      
      panel.grid.major.x =
        element_blank(),
      
      legend.title =
        element_text(
          face = "bold",
          size = 9.5
        ),
      
      legend.text =
        element_text(
          size = 9
        ),
      
      legend.position =
        "right",
      
      plot.margin =
        margin(
          15,
          15,
          15,
          15
        )
    )
  
  
  return(p)
}



# BLOCK 29 — CREATE ALL SIX FAMILY FIGURES


family_luki_onset_plot <-
  plot_family_forest(
    family_luki_onset,
    "Luki",
    "Onset"
  )


family_luki_end_plot <-
  plot_family_forest(
    family_luki_end,
    "Luki",
    "End"
  )


family_luki_duration_plot <-
  plot_family_forest(
    family_luki_duration,
    "Luki",
    "Duration"
  )


family_yangambi_onset_plot <-
  plot_family_forest(
    family_yangambi_onset,
    "Yangambi",
    "Onset"
  )


family_yangambi_end_plot <-
  plot_family_forest(
    family_yangambi_end,
    "Yangambi",
    "End"
  )


family_yangambi_duration_plot <-
  plot_family_forest(
    family_yangambi_duration,
    "Yangambi",
    "Duration"
  )



# BLOCK 30 — VIEW THE FIGURES


print(family_luki_onset_plot)
print(family_luki_end_plot)
print(family_luki_duration_plot)

print(family_yangambi_onset_plot)
print(family_yangambi_end_plot)
print(family_yangambi_duration_plot)



# BLOCK 31 — FINAL COMMUNITY + SPECIES SUMMARY TABLE


temporal_results_final <- pglmm_temporal_results %>%
  
  left_join(
    species_temporal_summary,
    by = c(
      "site",
      "response"
    )
  ) %>%
  
  select(
    
    site,
    response,
    n_species,
    
    beta,
    SE,
    Z,
    P,
    
    mean_slope,
    median_slope,
    
    Q1,
    Q3,
    IQR,
    SD,
    
    negative_n,
    positive_n,
    
    negative_percent,
    positive_percent
  )



# BLOCK 32 — ADD BIOLOGICAL DIRECTION AND SIGNIFICANCE


temporal_table <- temporal_results_final %>%
  
  mutate(
    
    Direction =
      case_when(
        
        response == "Onset" &
          beta < 0 ~
          "Earlier flowering",
        
        response == "Onset" &
          beta > 0 ~
          "Later flowering",
        
        response == "End" &
          beta < 0 ~
          "Earlier termination",
        
        response == "End" &
          beta > 0 ~
          "Later termination",
        
        response == "Duration" &
          beta < 0 ~
          "Shorter duration",
        
        response == "Duration" &
          beta > 0 ~
          "Longer duration"
      ),
    
    Significance =
      case_when(
        
        P < 0.001 ~ "***",
        
        P < 0.01 ~ "**",
        
        P < 0.05 ~ "*",
        
        TRUE ~ "ns"
      )
  )


print(temporal_table)



# BLOCK 33 — FORMATTED MANUSCRIPT TABLE


temporal_table_formatted <- temporal_table %>%
  
  mutate(
    
    beta =
      round(
        beta,
        4
      ),
    
    SE =
      round(
        SE,
        4
      ),
    
    Z =
      round(
        Z,
        3
      ),
    
    P =
      ifelse(
        P < 0.001,
        "<0.001",
        sprintf(
          "%.4f",
          P
        )
      ),
    
    negative_percent =
      round(
        negative_percent,
        1
      ),
    
    positive_percent =
      round(
        positive_percent,
        1
      )
  )


print(temporal_table_formatted)



# BLOCK 34 — EXPORT EXCEL WORKBOOK


wb <- createWorkbook()


write_sheet <- function(
    wb,
    sheet_name,
    data
) {
  
  addWorksheet(
    wb,
    sheet_name
  )
  
  
  header_style <- createStyle(
    textDecoration = "bold",
    halign = "center",
    border = "Bottom"
  )
  
  
  writeData(
    wb,
    sheet = sheet_name,
    x = data,
    headerStyle = header_style
  )
  
  
  freezePane(
    wb,
    sheet = sheet_name,
    firstRow = TRUE
  )
  
  
  setColWidths(
    wb,
    sheet = sheet_name,
    cols = seq_len(ncol(data)),
    widths = "auto"
  )
}


# Main temporal results

write_sheet(
  wb,
  "Temporal_Summary",
  temporal_table
)


write_sheet(
  wb,
  "Complete_Summary",
  temporal_results_final
)


# All species-specific slopes

write_sheet(
  wb,
  "All_Species_Slopes",
  all_species_results
)


# Luki family results

write_sheet(
  wb,
  "Luki_Onset",
  family_luki_onset
)


write_sheet(
  wb,
  "Luki_End",
  family_luki_end
)


write_sheet(
  wb,
  "Luki_Duration",
  family_luki_duration
)


# Yangambi family results

write_sheet(
  wb,
  "Yangambi_Onset",
  family_yangambi_onset
)


write_sheet(
  wb,
  "Yangambi_End",
  family_yangambi_end
)


write_sheet(
  wb,
  "Yangambi_Duration",
  family_yangambi_duration
)


# Combined family table

write_sheet(
  wb,
  "All_Family_Statistics",
  family_temporal_statistics
)


# Wilcoxon results

write_sheet(
  wb,
  "Wilcoxon_Species_Slopes",
  wilcoxon_species_slopes
)


# Family counts

write_sheet(
  wb,
  "Family_Counts",
  family_counts
)


# Save workbook

output_file <- file.path(
  output_folder,
  "Flowering_Temporal_and_Family_Results.xlsx"
)


saveWorkbook(
  wb,
  output_file,
  overwrite = TRUE
)


cat(
  "\nExcel results saved to:\n",
  output_file,
  "\n"
)



# BLOCK 35 — EXPORT LUKI FIGURES


ggsave(
  file.path(
    output_folder,
    "Figure_Family_Luki_Onset.pdf"
  ),
  family_luki_onset_plot,
  device = cairo_pdf,
  width = 9,
  height = 12,
  units = "in"
)


ggsave(
  file.path(
    output_folder,
    "Figure_Family_Luki_End.pdf"
  ),
  family_luki_end_plot,
  device = cairo_pdf,
  width = 9,
  height = 12,
  units = "in"
)


ggsave(
  file.path(
    output_folder,
    "Figure_Family_Luki_Duration.pdf"
  ),
  family_luki_duration_plot,
  device = cairo_pdf,
  width = 9,
  height = 12,
  units = "in"
)



# BLOCK 36 — EXPORT YANGAMBI FIGURES


ggsave(
  file.path(
    output_folder,
    "Figure_Family_Yangambi_Onset.pdf"
  ),
  family_yangambi_onset_plot,
  device = cairo_pdf,
  width = 9,
  height = 23,
  units = "in"
)


ggsave(
  file.path(
    output_folder,
    "Figure_Family_Yangambi_End.pdf"
  ),
  family_yangambi_end_plot,
  device = cairo_pdf,
  width = 9,
  height = 23,
  units = "in"
)


ggsave(
  file.path(
    output_folder,
    "Figure_Family_Yangambi_Duration.pdf"
  ),
  family_yangambi_duration_plot,
  device = cairo_pdf,
  width = 9,
  height = 23,
  units = "in"
)



# BLOCK 37 — FINAL CHECKS


cat("\n==============================================\n")
cat("FINAL ANALYSIS CHECK\n")
cat("==============================================\n")

cat(
  "\nLuki:",
  n_distinct(luki_onset$species),
  "species |",
  n_distinct(luki_onset_family$family),
  "families\n"
)

cat(
  "Yangambi:",
  n_distinct(yangambi_onset$species),
  "species |",
  n_distinct(yangambi_onset_family$family),
  "families\n"
)

cat(
  "\nResults folder:\n",
  output_folder,
  "\n"
)

cat("\nAnalysis completed.\n")