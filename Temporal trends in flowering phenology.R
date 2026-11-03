
library(phyr)

# 12 MODELS: LUKI AND YANGAMBI


# A. LUKI

luki_raw$year_f <- as.factor(luki_raw$year)


# 1. Onset WITH year random intercept

luki_onset_year1 <- phyr::pglmm(
  onset_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = luki_raw,
  family = "gaussian",
  cov_ranef = list(species = phy),
  REML = FALSE
)

summary(luki_onset_year1)

# summary(luki_onset_year1)
# Linear mixed model fit by maximum likelihood
# 
# Call:onset_doy ~ year_c
# 
# logLik    AIC    BIC 
#  -6109  12234  12265 
# 
# Random effects:
#                       Variance   Std.Dev
# 1|species            2.588e+03   50.8769
# 1|species__          4.132e+02   20.3282
# 0 + year_c|species   1.475e+01    3.8404
# 0 + year_c|species__ 5.399e-02    0.2324
# 1|year               6.085e+07 7800.8189
# residual             5.730e+03   75.6975
# 
# Fixed effects:
#                 Value Std.Error  Zscore Pvalue
# (Intercept)  146.7171 4400.2757  0.0333 0.9734
# year_c        -7.4721  743.7794 -0.0100 0.9920


# 2. Onset WITHOUT year random intercept

luki_onset_no_year2 <- phyr::pglmm(
  onset_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__),
  data = luki_raw,
  family = "gaussian",
  cov_ranef = list(species = phy),
  REML = FALSE
)

summary(luki_onset_no_year2)


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


# WHY MODEL 1 WAS REJECTED 

# STRUCTURAL COLLINEARITY: 
#    Including both a continuous linear trend (year_c) and a random intercept 
#    (1|year)  forces them to compete for the exact same variance.

# OPTIMISER BREAKDOWN: 
#    The model failed to distinguish between the overall trend and random annual 
#    noise. This caused the 1|year variance to explode to an artificial 
#    60.85 million witha SD of over 7,800 days.

# LOSS OF STATISTICAL POWER: 
#    The fixed-effect standard error for year_c jump from 1.10 to 743.78.

#    We rejected the lower AIC because the model parameters 
#    are physically and biologically impossible. An intercept standard deviation of 
#    ~7,800 days on a variable between 1 and 365 proof that the model's 
#    math is broken, making its AIC score completely invalid for selection.


# 3. Termination WITH year random intercept

luki_end_year1 <- phyr::pglmm(
  end_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = luki_raw,
  family = "gaussian",
  cov_ranef = list(species = phy),
  REML = FALSE
)

summary(luki_end_year1)

# summary(luki_end_year1)
# Linear mixed model fit by maximum likelihood
# 
# Call:end_doy ~ year_c
# 
# logLik    AIC    BIC 
#  -5967  11951  11981 
# 
# Random effects:
#                       Variance   Std.Dev
# 1|species            2.683e+03 5.180e+01
# 1|species__          4.205e+02 2.051e+01
# 0 + year_c|species   1.144e+01 3.382e+00
# 0 + year_c|species__ 1.015e-03 3.185e-02
# 1|year               3.915e+07 6.257e+03
# residual             4.214e+03 6.491e+01
# 
# Fixed effects:
#                 Value Std.Error  Zscore Pvalue
# (Intercept)  282.6793 3529.3895  0.0801 0.9362
# year_c        -3.0206  596.5719 -0.0051 0.9960

# DIAGNOSTIC NOTE:
# This  luki_end_year1 suffers from the exact same mathematical 
# collapse as the onset model. The 1|year variance exploded to 3.915e+07 
# Std.Dev ~6,257 days, blowing out the year_c Standard Error to 596.57.


# 4. Termination WITHOUT year random intercept

luki_end_no_year2 <- phyr::pglmm(
  end_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__),
  data = luki_raw,
  family = "gaussian",
  cov_ranef = list(species = phy),
  REML = FALSE
)

summary(luki_end_no_year2)




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




# 5. Duration WITH year random intercept

luki_duration_year1 <- phyr::pglmm(
  duration_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = luki_raw,
  family = "gaussian",
  cov_ranef = list(species = phy),
  REML = FALSE
)

summary(luki_duration_year1)


# summary(luki_duration_year1)
# Linear mixed model fit by maximum likelihood
# 
# Call:duration_doy ~ year_c
# 
# logLik    AIC    BIC 
#  -6279  12573  12604 
# 
# Random effects:
#                       Variance   Std.Dev
# 1|species            2.623e+03 5.121e+01
# 1|species__          3.184e+01 5.643e+00
# 0 + year_c|species   3.601e+01 6.000e+00
# 0 + year_c|species__ 1.371e-03 3.703e-02
# 1|year               8.791e+07 9.376e+03
# residual             8.212e+03 9.062e+01
# 
# Fixed effects:
#                 Value Std.Error Zscore Pvalue
# (Intercept)  138.2054 5288.7441 0.0261 0.9792
# year_c         4.2285  893.9606 0.0047 0.9962

# DIAGNOSTIC NOTE:
# luki_duration_year1 has an even worse optimization failure, with the 
# 1|year variance hitting an astronomical 8.791e+07 SD ~9,376 days. 
# year_c Standard Error to 893.96, 

# 6. Duration WITHOUT year random intercept

luki_duration_no_year2 <- phyr::pglmm(
  duration_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__),
  data = luki_raw,
  family = "gaussian",
  cov_ranef = list(species = phy),
  REML = FALSE
)

summary(luki_duration_no_year2)


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


# 7. Onset WITH year random intercept

yangambi_onset_year1 <- phyr::pglmm(
  onset_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = yangambi_df,
  family = "gaussian",
  cov_ranef = list(species = yangambi_phy),
  REML = FALSE
)

summary(yangambi_onset_year1)


# summary(yangambi_onset_year1)
# Linear mixed model fit by maximum likelihood
# 
# Call:onset_doy ~ year_c
# 
# logLik    AIC    BIC 
# -24320  48657  48699 
# 
# Random effects:
#                       Variance   Std.Dev
# 1|species            7.203e+02   26.8381
# 1|species__          9.096e+00    3.0160
# 0 + year_c|species   1.189e+00    1.0902
# 0 + year_c|species__ 2.248e-01    0.4742
# 1|year               3.957e+07 6290.1018
# residual             3.715e+03   60.9536
# 
# Fixed effects:
#                 Value Std.Error  Zscore Pvalue
# (Intercept)  120.0648 2775.7615  0.0433 0.9655
# year_c        -3.4281  263.4638 -0.0130 0.9896

# 8. Onset WITHOUT year random intercept

yangambi_onset_no_year2 <- phyr::pglmm(
  onset_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__),
  data = yangambi_df,
  family = "gaussian",
  cov_ranef = list(species = yangambi_phy),
  REML = FALSE
)

summary(yangambi_onset_no_year2)

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


# 9. Termination WITH year random intercept

yangambi_end_year1 <- phyr::pglmm(
  end_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = yangambi_df,
  family = "gaussian",
  cov_ranef = list(species = yangambi_phy),
  REML = FALSE
)

summary(yangambi_end_year1)


# summary(yangambi_end_year1)
# Linear mixed model fit by maximum likelihood
# 
# Call:end_doy ~ year_c
# 
# logLik    AIC    BIC 
# -25358  50732  50774 
# 
# Random effects:
#                       Variance   Std.Dev
# 1|species            1.738e+03 4.169e+01
# 1|species__          1.495e+02 1.223e+01
# 0 + year_c|species   9.490e+00 3.081e+00
# 0 + year_c|species__ 1.015e-01 3.186e-01
# 1|year               1.073e+08 1.036e+04
# residual             5.515e+03 7.426e+01
# 
# Fixed effects:
#                 Value Std.Error Zscore Pvalue
# (Intercept)  283.5049 4571.3856  0.062 0.9505
# year_c        -2.6239  433.8959 -0.006 0.9952



# 10. Termination WITHOUT year random intercept

yangambi_end_no_year2 <- phyr::pglmm(
  end_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__),
  data = yangambi_df,
  family = "gaussian",
  cov_ranef = list(species = yangambi_phy),
  REML = FALSE
)

summary(yangambi_end_no_year2)

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


# 11. Duration WITH year random intercept

yangambi_duration_year1 <- phyr::pglmm(
  duration_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__) +
    (1 | year),
  data = yangambi_df,
  family = "gaussian",
  cov_ranef = list(species = yangambi_phy),
  REML = FALSE
)

summary(yangambi_duration_year1)

# summary(yangambi_duration_year1)
# Linear mixed model fit by maximum likelihood
# 
# Call:duration_doy ~ year_c
# 
# logLik    AIC    BIC 
# -25805  51627  51668 
# 
# Random effects:
#                       Variance   Std.Dev
# 1|species            2.266e+03    47.600
# 1|species__          1.947e+02    13.953
# 0 + year_c|species   1.203e+01     3.469
# 0 + year_c|species__ 1.115e+00     1.056
# 1|year               1.589e+08 12606.780
# residual             6.670e+03    81.669
# 
# Fixed effects:
#                 Value Std.Error Zscore Pvalue
# (Intercept)  159.2396 5563.2578 0.0286 0.9772
# year_c         1.1496  528.0408 0.0022 0.9983


# 12. Duration WITHOUT year random intercept

yangambi_duration_no_year2 <- phyr::pglmm(
  duration_doy ~ year_c +
    (1 | species__) +
    (0 + year_c | species__),
  data = yangambi_df,
  family = "gaussian",
  cov_ranef = list(species = yangambi_phy),
  REML = FALSE
)



summary(yangambi_duration_no_year2)


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



luki_onset_no_year2
luki_end_no_year2
luki_duration_no_year2
yangambi_onset_no_year2
yangambi_end_no_year2
yangambi_duration_no_year2
