

##Yangambi_onset_additive for yangambi_onset
##yangambi_duration_additive for yangambi_duration


# LUKI DATA MODELS


luki_onset_additive <- pglmm(
  onset_doy ~ temp_z + precip_z + (1 | species) +
    (1 | year),
  data = luki_raw, 
  family = "gaussian",
  cov_ranef = list(species = phy), 
  REML = FALSE
)

summary(luki_onset_additive)

# summary(luki_onset_additive)
# Linear mixed model fit by maximum likelihood
# 
# Call:onset_doy ~ temp_z + precip_z
# 
# logLik    AIC    BIC 
#  -6059  12131  12153 
# 
# Random effects:
#           Variance Std.Dev
# 1|species   3444.3   58.69
# 1|year       981.7   31.33
# residual    6042.5   77.73
# 
# Fixed effects:
#                Value Std.Error  Zscore  Pvalue    
# (Intercept) 118.9325   11.1342 10.6817 < 2e-16 ***
# temp_z       17.9776   10.1897  1.7643 0.07768 .  
# precip_z      3.1993    9.9676  0.3210 0.74824    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

luki_onset_interaction <- pglmm(
  onset_doy ~ temp_z * precip_z + (1 | species) +
    (1 | year),
  data = luki_raw, 
  family = "gaussian",
  cov_ranef = list(species = phy), 
  REML = FALSE
)

summary(luki_onset_interaction)
# summary(luki_onset_interaction)
# Linear mixed model fit by maximum likelihood
# 
# Call:onset_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
#  -6058  12130  12157 
# 
# Random effects:
#           Variance Std.Dev
# 1|species     3422   58.50
# 1|year         738   27.17
# residual      6048   77.77
# 
# Fixed effects:
#                    Value Std.Error Zscore  Pvalue    
# (Intercept)     110.2155   11.2769 9.7736 < 2e-16 ***
# temp_z           21.4792    9.1559 2.3459 0.01898 *  
# precip_z          4.6933    8.7824 0.5344 0.59307    
# temp_z:precip_z  30.2494   17.5536 1.7233 0.08484 .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Model End phenology


luki_end_additive <- pglmm(
  end_doy ~ temp_z + precip_z + (1 | species) +
    (1 | year),
  data = luki_raw, 
  family = "gaussian",
  cov_ranef = list(species = phy), 
  REML = FALSE
)

summary(luki_end_additive)

# summary(luki_end_additive)
# Linear mixed model fit by maximum likelihood
# 
# Call:end_doy ~ temp_z + precip_z
# 
# logLik    AIC    BIC 
#  -5914  11841  11863 
# 
# Random effects:
#           Variance Std.Dev
# 1|species   3828.1   61.87
# 1|year       419.9   20.49
# residual    4387.9   66.24
# 
# Fixed effects:
#                 Value Std.Error  Zscore Pvalue    
# (Intercept) 275.94511   8.58536 32.1414 <2e-16 ***
# temp_z        0.63408   6.81134  0.0931 0.9258    
# precip_z      4.68132   6.67132  0.7017 0.4829    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

luki_end_interaction <- pglmm(
  end_doy ~ temp_z * precip_z + (1 | species) +
    (1 | year),
  data = luki_raw, 
  family = "gaussian",
  cov_ranef = list(species = phy), 
  REML = FALSE
)

summary(luki_end_interaction)

# summary(luki_end_interaction)
# Linear mixed model fit by maximum likelihood
# 
# Call:end_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
#  -5913  11840  11866 
# 
# Random effects:
#           Variance Std.Dev
# 1|species   3854.0   62.08
# 1|year       295.2   17.18
# residual    4385.9   66.23
# 
# Fixed effects:
#                    Value Std.Error  Zscore  Pvalue    
# (Intercept)     269.8673    8.5661 31.5041 < 2e-16 ***
# temp_z            3.1334    5.9847  0.5236 0.60057    
# precip_z          5.6599    5.7499  0.9843 0.32495    
# temp_z:precip_z  21.2658   11.4768  1.8529 0.06389 .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##Model Duration

luki_duration_additive <- pglmm(
  duration_doy ~ temp_z + precip_z + (1 | species) +
    (1 | year),
  data = luki_raw, 
  family = "gaussian",
  cov_ranef = list(species = phy), 
  REML = FALSE
)

summary(luki_duration_additive)

# summary(luki_duration_additive)
# Linear mixed model fit by maximum likelihood
# 
# Call:duration_doy ~ temp_z + precip_z
# 
# logLik    AIC    BIC 
#  -6221  12453  12476 
# 
# Random effects:
#           Variance Std.Dev
# 1|species   3029.8   55.04
# 1|year       146.5   12.10
# residual    8828.0   93.96
# 
# Fixed effects:
#                Value Std.Error  Zscore    Pvalue    
# (Intercept) 159.5269    6.8869 23.1638 < 2.2e-16 ***
# temp_z      -17.0588    4.9075 -3.4761 0.0005088 ***
# precip_z      1.1005    4.8468  0.2271 0.8203771    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

luki_duration_interaction <- pglmm(
  duration_doy ~ temp_z * precip_z + (1 | species) +
    (1 | year),
  data = luki_raw, 
  family = "gaussian",
  cov_ranef = list(species = phy), 
  REML = FALSE
)

summary(luki_duration_interaction)

# summary(luki_duration_interaction)
# Linear mixed model fit by maximum likelihood
# 
# Call:duration_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
#  -6220  12454  12481 
# 
# Random effects:
#           Variance Std.Dev
# 1|species   3081.0   55.51
# 1|year       120.7   10.99
# residual    8816.3   93.89
# 
# Fixed effects:
#                     Value Std.Error  Zscore    Pvalue    
# (Intercept)     162.05434   7.22212 22.4386 < 2.2e-16 ***
# temp_z          -18.08940   4.76057 -3.7998 0.0001448 ***
# precip_z          0.63041   4.60727  0.1368 0.8911664    
# temp_z:precip_z  -9.19738   9.14376 -1.0059 0.3144807    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# YANGAMBI DATA MODELS

#Onset

yangambi_onset_additive <- pglmm(
  onset_doy ~ temp_z + precip_z + (1 | species) +
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
# -24243  48497  48529 
# 
# Random effects:
#           Variance Std.Dev
# 1|species    911.8   30.20
# 1|year       849.7   29.15
# residual    3789.7   61.56
# 
# Fixed effects:
#               Value Std.Error  Zscore Pvalue    
# (Intercept) 89.8819    7.1032 12.6538 <2e-16 ***
# temp_z      -6.3817    6.3386 -1.0068 0.3140    
# precip_z     1.4552    6.9058  0.2107 0.8331    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

yangambi_onset_interaction <- pglmm(
  onset_doy ~ temp_z * precip_z + (1 | species) +
    (1 | year),
  data = yangambi_df, 
  family = "gaussian",
  cov_ranef = list(species = yangambi_phy), 
  REML = FALSE
)

summary(yangambi_onset_interaction)

# summary(yangambi_onset_interaction)
# Linear mixed model fit by maximum likelihood
# 
# Call:onset_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
# -24242  48499  48536 
# 
# Random effects:
#           Variance Std.Dev
# 1|species    894.4   29.91
# 1|year       880.5   29.67
# residual    3795.1   61.60
# 
# Fixed effects:
#                   Value Std.Error  Zscore Pvalue    
# (Intercept)     89.7227    7.2243 12.4196 <2e-16 ***
# temp_z          -9.0055    7.9566 -1.1318 0.2577    
# precip_z         1.9770    7.0894  0.2789 0.7803    
# temp_z:precip_z  4.1648    7.3983  0.5629 0.5735    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Model end 

yangambi_end_additive <- pglmm(
  end_doy ~ temp_z + precip_z + (1 | species) +
    (1 | year),
  data = yangambi_df, 
  family = "gaussian",
  cov_ranef = list(species = yangambi_phy), 
  REML = FALSE
)

summary(yangambi_end_additive)

# summary(yangambi_end_additive)
# Linear mixed model fit by maximum likelihood
# 
# Call:end_doy ~ temp_z + precip_z
# 
# logLik    AIC    BIC 
# -25277  50567  50598 
# 
# Random effects:
#           Variance Std.Dev
# 1|species   3064.3   55.36
# 1|year       290.5   17.04
# residual    5725.0   75.66
# 
# Fixed effects:
#                 Value Std.Error  Zscore Pvalue    
# (Intercept) 254.11674   4.84909 52.4051 <2e-16 ***
# temp_z        0.56227   3.87833  0.1450 0.8847    
# precip_z     -0.19682   4.19351 -0.0469 0.9626    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

yangambi_end_interaction <- pglmm(
  end_doy ~ temp_z * precip_z + (1 | species) +
    (1 | year),
  data = yangambi_df, 
  family = "gaussian",
  cov_ranef = list(species = yangambi_phy), 
  REML = FALSE
)

summary(yangambi_end_interaction)
# summary(yangambi_end_interaction)
# Linear mixed model fit by maximum likelihood
# 
# Call:end_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
# -25277  50569  50605 
# 
# Random effects:
#           Variance Std.Dev
# 1|species   3062.6   55.34
# 1|year       290.3   17.04
# residual    5725.3   75.67
# 
# Fixed effects:
#                     Value Std.Error  Zscore Pvalue    
# (Intercept)     254.09886   4.84955 52.3964 <2e-16 ***
# temp_z            0.13951   4.75712  0.0293 0.9766    
# precip_z         -0.12480   4.21850 -0.0296 0.9764    
# temp_z:precip_z   0.68032   4.43438  0.1534 0.8781    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Model duration

yangambi_duration_additive <- pglmm(
  duration_doy ~ temp_z + precip_z + (1 | species) +
    (1 | year),
  data = yangambi_df, 
  family = "gaussian",
  cov_ranef = list(species = yangambi_phy), 
  REML = FALSE
)


# summary(yangambi_duration_additive)
# Linear mixed model fit by maximum likelihood
# 
# Call:duration_doy ~ temp_z + precip_z
# 
# logLik    AIC    BIC 
# -25744  51499  51531 
# 
# Random effects:
#           Variance Std.Dev
# 1|species     4188   64.71
# 1|year         320   17.89
# residual      7030   83.84
# 
# Fixed effects:
#                Value Std.Error  Zscore Pvalue    
# (Intercept) 162.3399    5.2484 30.9311 <2e-16 ***
# temp_z        6.5289    4.0963  1.5939 0.1110    
# precip_z     -2.0809    4.4248 -0.4703 0.6382    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

yangambi_duration_interaction <- pglmm(
  duration_doy ~ temp_z * precip_z + (1 | species) +
    (1 | year),
  data = yangambi_df, 
  family = "gaussian",
  cov_ranef = list(species = yangambi_phy), 
  REML = FALSE
)


summary(yangambi_duration_interaction)

# summary(yangambi_duration_interaction)
# Linear mixed model fit by maximum likelihood
# 
# Call:duration_doy ~ temp_z * precip_z
# 
# logLik    AIC    BIC 
# -25743  51501  51538 
# 
# Random effects:
#           Variance Std.Dev
# 1|species     4177   64.63
# 1|year         313   17.69
# residual      7032   83.86
# 
# Fixed effects:
#                    Value Std.Error  Zscore  Pvalue    
# (Intercept)     162.4549    5.2108 31.1767 < 2e-16 ***
# temp_z            8.3759    4.9730  1.6843 0.09213 .  
# precip_z         -2.3868    4.4067 -0.5416 0.58807    
# temp_z:precip_z  -3.0005    4.6374 -0.6470 0.51762    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
