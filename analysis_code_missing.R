# R script with supplemental analysis code reported in:
# 'Effects of time-varying parent input on child language outcomes differ for vocabulary and syntax'

# The following runs the analysis reported in Text S2, using only complete data

# get required packages
required_packages <- c("dplyr", "MASS", "lmtest", "sandwich", "mice", "ggplot2", "wiqid", "lm.beta")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# dplyr for making quantiles
library(dplyr)
# MASS for running ordinal regressions
library(MASS)
# lmtest and sandwich libraries for robust regressions
library(lmtest)
library(sandwich)
# mice for multiple imputation
library(mice)
# ggplot2 for plotting
library(ggplot2)
# wiqid for AICc
library(wiqid)
# lm.beta for returning standardised regression coefficients
library(lm.beta)

# define function for producing robust regression output
# (reproduces robust results from Stata regressions)
# takes a model x as argument
robust <- function(x) {
  coeftest(x, vcov = vcovHC(x, "HC1"))
}

# APPROACH USING ONLY COMPLETE CASES

# read in raw data
lang_data_raw <- read.table("fake_lang_data.txt", header=TRUE)

# add completeness indicators complete_v and complete_s

# if any missing values in predictors, covariates, or vocabulary outcome, complete_v = 0
# if any missing values in predictors, covariates, or syntax outcome, complete_s = 0

lang_data_raw$complete_v <- ifelse(rowSums(is.na(lang_data_raw[ , c(1:8,11,13,15,21,22,25)]))>0, 0,1)
lang_data_raw$complete_s <- ifelse(rowSums(is.na(lang_data_raw[ , c(1:8,11,13,15,21,22,26)]))>0, 0,1)


# gender and birth order are factors
lang_data_raw$c_gen <- factor(lang_data_raw$c_gen)
lang_data_raw$c_bo <- factor(lang_data_raw$c_bo)

# calculate compos26 and other scales

# composite measure of child language at 26 months: average z-score of word types and MLU
lang_data_raw$compos26 <- (scale(lang_data_raw$c_wt_26) + scale(lang_data_raw$c_mlu_26))/2

# Z1_syntax_scaled and Z2_syntax_scaled: raw scores * 100
lang_data_raw$Z1_syntax_rescale <- lang_data_raw$p_cps_14 * 100
lang_data_raw$Z2_syntax_rescale <- lang_data_raw$p_cps_30 * 100

# strata: 8 quantiles

lang_data_raw$Z1stratum_vocab <- factor(ntile(lang_data_raw$p_wt_14, 8))
lang_data_raw$Z2stratum_vocab <- factor(ntile(lang_data_raw$p_wt_30, 8))
lang_data_raw$Z1stratum_syntax <- factor(ntile(lang_data_raw$p_cps_14, 8))
lang_data_raw$Z2stratum_syntax <- factor(ntile(lang_data_raw$p_cps_30, 8))

# cumulative Z1 and Z2: Z1 + Z2

lang_data_raw$cumulative_vocab <- lang_data_raw$p_wt_14 + lang_data_raw$p_wt_30
lang_data_raw$cumulative_syntax <- lang_data_raw$Z1_syntax_rescale + lang_data_raw$Z2_syntax_rescale

# what proportion of data are complete for vocab? 0.75
mean(lang_data_raw$complete_v)

lang_data_raw$complete_v <- factor(lang_data_raw$complete_v)

# logistic regression model predicting missingness from all covariates in our main analysis for which we have complete data

miss_model_v <- glm(complete_v ~ p_wt_14 + c_bo + c_gen + c_wt_14 + c_gt_14 + p_inc +
                      p_educ, lang_data_raw, family = binomial)

# get predicted probabilities of response from this model

lang_data_raw$v_resp_prob <- predict(miss_model_v, type="response")

# calculate weights: numerator = average completeness for vocab, denominator = person-specific probability of completeness
# for subjects who are less likely than average to be complete, this will upweight them
# for subjects who are more likely than average to be complete, this will downweight them
lang_data_raw$v_resp_w <- 0.75/lang_data_raw$v_resp_prob


# now, we calculate covariate-based weights using same covariates in our models as we include in the main analysis

# Z1 vocab
# predict stratum assignment from covariates
denominator_Z1_vocab <- polr(Z1stratum_vocab ~ p_viq + p_inc + c_gen, lang_data_raw, na.action = na.exclude)

# save the probabilities of being in each stratum as a new data frame
denprobs_Z1_vocab <- data.frame(predict(denominator_Z1_vocab, type="p"))

# save the probability for each person of being in their own stratum
# first, make the column names match by assigning each stratum number to the appropriate column heading
colnames(denprobs_Z1_vocab) <- seq_len(8)
# add veridical stratum to this data frame
denprobs_Z1_vocab['Z1stratum_vocab'] <- lang_data_raw$Z1stratum_vocab
# personprob = probability for each person of being in their actual stratum as listed in the Zstratum column
lang_data_raw['denprob_Z1_vocab'] <- as.numeric(denprobs_Z1_vocab[cbind(seq_len(nrow(denprobs_Z1_vocab)), match(denprobs_Z1_vocab$Z1stratum_vocab, names(denprobs_Z1_vocab)))])

# calculate the weights
# numerator = marginal probability of being in each stratum, i.e. 1/8
# denominator = person-specific probability of being in their stratum
# qbw = quantile binning weight

lang_data_raw['Z1qbw_vocab'] <- 1/8/(lang_data_raw$denprob_Z1_vocab)

# multiply these weights by the missingness weights
lang_data_raw$Z1qbw_v_final <- lang_data_raw$Z1qbw_vocab * lang_data_raw$v_resp_w


# now we calculate covariate-based weights for vocab Z2

denominator_Z2_vocab <- polr(Z2stratum_vocab ~ p_wt_14 + compos26, lang_data_raw, na.action = na.exclude)

# save the probabilities of being in each stratum as a new data frame
denprobs_Z2_vocab <- data.frame(predict(denominator_Z2_vocab, type="p"))

# save the probability for each person of being in their own stratum
# first, make the column names match by assigning each stratum number to the appropriate column heading
colnames(denprobs_Z2_vocab) <- seq_len(8)
# add actual stratum to this data frame
denprobs_Z2_vocab['Z2stratum_vocab'] <- lang_data_raw$Z2stratum_vocab
# personprob = probability for each person of being in their actual stratum as listed in the Zstratum column
lang_data_raw['denprob_Z2_vocab'] <- as.numeric(denprobs_Z2_vocab[cbind(seq_len(nrow(denprobs_Z2_vocab)), match(denprobs_Z2_vocab$Z2stratum_vocab, names(denprobs_Z2_vocab)))])

# Now, we calculate the numerator of the Z2 weight.
# numerator = probability of being in Z2 stratum given Z1
# So it's an ordinal model like the denominator model, but the only predictor is Z1 (no covariates)

numerator_Z2_vocab <- polr(Z2stratum_vocab ~ p_wt_14, lang_data_raw, na.action = na.exclude)
# save the probabilities of being in each stratum as a new data frame
numprobs_Z2_vocab <- data.frame(predict(numerator_Z2_vocab, type="p"))
# save the probability for each person of being in their own stratum
# first, make the column names match by assigning each stratum number to the appropriate column heading
colnames(numprobs_Z2_vocab) <- seq_len(8)
# add actual stratum to this data frame
numprobs_Z2_vocab['Z2stratum_vocab'] <- lang_data_raw$Z2stratum_vocab
# numprob = probability for each person of being in their actual stratum as listed in the Zstratum column
lang_data_raw['numprob_Z2_vocab'] <- as.numeric(numprobs_Z2_vocab[cbind(seq_len(nrow(numprobs_Z2_vocab)), match(numprobs_Z2_vocab$Z2stratum_vocab, names(numprobs_Z2_vocab)))])

# divide the numerator probability by the denominator probability to get the weight

lang_data_raw['Z2qbw_vocab'] <- lang_data_raw$numprob_Z2_vocab/lang_data_raw$denprob_Z2_vocab

# combined weights (missingness weights already incorporated into Z1 final weights)
lang_data_raw$combqbw_v_final <- lang_data_raw$Z1qbw_v_final * lang_data_raw$Z2qbw_vocab


# Now we run outcome models for vocabulary

# Model with Z1 and Z2 as separate predictors (row 1 of Table S1)

v_sep <- lm(c_ppvt_74 ~ p_wt_14 + p_wt_30, weights = combqbw_v_final, lang_data_raw)

v_m2 <- robust(v_sep)

print(v_m2)

# confidence intervals

v_m2_z1_lowlim <- v_m2[2] - (qt(0.975, v_sep$df.residual) * v_m2[5])
v_m2_z1_uplim <- v_m2[2] + (qt(0.975, v_sep$df.residual) * v_m2[5])
print(v_m2_z1_lowlim)
print(v_m2_z1_uplim)

v_m2_z2_lowlim <- v_m2[3] - (qt(0.975, v_sep$df.residual) * v_m2[6])
v_m2_z2_uplim <- v_m2[3] + (qt(0.975, v_sep$df.residual) * v_m2[6])
print(v_m2_z2_lowlim)
print(v_m2_z2_uplim)

# AICc for this model
v_sep_AIC <- AICc(v_sep)

print(v_sep_AIC)


# hypothesis test: are effects of Z1 and Z2 equal?

# to answer this question, model cumulative + Z2 and look at coefficient for Z2

v_diff <- lm(c_ppvt_74 ~ cumulative_vocab + p_wt_30, weights = combqbw_v_final, lang_data_raw)

print(robust(v_diff))

coef_diff <- v_sep$coefficients[['p_wt_30']] - v_sep$coefficients[['p_wt_14']]
var_covar <- vcov(v_sep)
diff_se <- sqrt(var_covar[2,2] + var_covar[3,3] - 2 * var_covar[2,3])
lower_bound <- coef_diff - (diff_se * 2)
upper_bound <- coef_diff + (diff_se * 2)

print(lower_bound)
print(upper_bound)

# Constant effects model (row 3 of Table S1)

v_const <- lm(c_ppvt_74 ~ cumulative_vocab, weights = combqbw_v_final, lang_data_raw)

# save model summary as v_m3
v_m3 <- robust(v_const)
print(v_m3)


# confidence interval

v_m3_cumul_lowlim <- v_m3[2] - (qt(0.975, v_const$df.residual) * v_m3[4])
v_m3_cumul_uplim <- v_m3[2] + (qt(0.975, v_const$df.residual) * v_m3[4])
print(v_m3_cumul_lowlim)
print(v_m3_cumul_uplim)


# AICc for this model
v_const_AIC <- AICc(v_const)

print(v_const_AIC)


# First extreme model: Z1 > 0, Z2 = 0 (row 4 of Table S1)

v_Z1 <- lm(c_ppvt_74 ~ p_wt_14, weights = combqbw_v_final, lang_data_raw)

# Save model summary as v_m4
v_m4 <- robust(v_Z1)
print(v_m4)

# confidence interval

v_m4_z1_lowlim <- v_m4[2] - (qt(0.975, v_Z1$df.residual) * v_m4[4])
v_m4_z1_uplim <- v_m4[2] + (qt(0.975, v_Z1$df.residual) * v_m4[4])
print(v_m4_z1_lowlim)
print(v_m4_z1_uplim)


# AICc for this model
v_Z1_AIC <- AICc(v_Z1)

print(v_Z1_AIC)

# Second 'extreme' model: Z1 = 0, Z2 > 0 (row 5 in Table S1)

v_Z2 <- lm(c_ppvt_74 ~ p_wt_30, weights = combqbw_v_final, lang_data_raw)

# Save model summary as v_m5
v_m5 <- robust(v_Z2)
print(v_m5)

# confidence interval

v_m5_z2_lowlim <- v_m5[2] - (qt(0.975, v_Z2$df.residual) * v_m5[4])
v_m5_z2_uplim <- v_m5[2] + (qt(0.975, v_Z2$df.residual) * v_m5[4])
print(v_m5_z2_lowlim)
print(v_m5_z2_uplim)


# AICc for this model
v_Z2_AIC <- AICc(v_Z2)

print(v_Z2_AIC)


# Now we do the same process for syntax

# what proportion data complete for syntax? 0.6875
mean(lang_data_raw$complete_s)

lang_data_raw$complete_s <- factor(lang_data_raw$complete_s)

# response model

miss_model_s <- glm(complete_s ~ p_cps_14 + c_bo + c_gen + c_wt_14 + c_gt_14 + p_inc +
                      p_educ, lang_data_raw, family = "binomial")

lang_data_raw$s_resp_prob <- predict(miss_model_s, type="response")

# calculate weights by adding numerator of average v missingness
lang_data_raw$s_resp_w <- 0.6875/lang_data_raw$s_resp_prob


# Now we calculate covariate-based weights

# Z1 syntax
# predict stratum assignment from covariates

denominator_Z1_syntax <- polr(Z1stratum_syntax ~ c_wt_14 + c_bo + p_viq, lang_data_raw, na.action = na.exclude)

# save the probabilities of being in each stratum as a new data frame
denprobs_Z1_syntax <- data.frame(predict(denominator_Z1_syntax, type="p"))

# save the probability for each person of being in their own stratum
# first, make the column names match by assigning each stratum number to the appropriate column heading
colnames(denprobs_Z1_syntax) <- seq_len(8)
# add actual stratum to this data frame
denprobs_Z1_syntax['Z1stratum_syntax'] <- lang_data_raw$Z1stratum_syntax
# personprob = probability for each person of being in their actual stratum as listed in the Zstratum column
lang_data_raw['denprob_Z1_syntax'] <- as.numeric(denprobs_Z1_syntax[cbind(seq_len(nrow(denprobs_Z1_syntax)), match(denprobs_Z1_syntax$Z1stratum_syntax, names(denprobs_Z1_syntax)))])

# calculate the weights
# numerator = marginal probability of being in each stratum, i.e. 1/8
# denominator = person-specific probability of being in their stratum
# qbw = quantile binning weight
lang_data_raw['Z1qbw_syntax'] <- 1/8/(lang_data_raw$denprob_Z1_syntax)

# multiply these weights by the missingness weights
lang_data_raw$Z1qbw_s_final <- lang_data_raw$Z1qbw_syntax * lang_data_raw$s_resp_w


# Now we clculate covariate-based weights for syntax Z2

denominator_Z2_syntax <- polr(Z2stratum_syntax ~ p_viq + p_cps_14 + compos26 + p_educ + c_gen, lang_data_raw,
                               na.action = na.exclude)

# save the probabilities of being in each stratum as a new data frame
denprobs_Z2_syntax <- data.frame(predict(denominator_Z2_syntax, type="p"))

# save the probability for each person of being in their own stratum
# first, make the column names match by assigning each stratum number to the appropriate column heading
colnames(denprobs_Z2_syntax) <- seq_len(8)
# add actual stratum to this data frame
denprobs_Z2_syntax['Z2stratum_syntax'] <- lang_data_raw$Z2stratum_syntax
# personprob = probability for each person of being in their actual stratum as listed in the Zstratum column
lang_data_raw['denprob_Z2_syntax'] <- as.numeric(denprobs_Z2_syntax[cbind(seq_len(nrow(denprobs_Z2_syntax)), match(denprobs_Z2_syntax$Z2stratum_syntax, names(denprobs_Z2_syntax)))])

# numerator = probability of being in each Z2 stratum given Z1
numerator_Z2_syntax <- polr(Z2stratum_syntax ~ p_cps_14, lang_data_raw, na.action = na.exclude)
# save the probabilities of being in each stratum as a new data frame
numprobs_Z2_syntax <- data.frame(predict(numerator_Z2_syntax, type="p"))

# save the probability for each person of being in their own stratum
# first, make the column names match by assigning each stratum number to the appropriate column heading
colnames(numprobs_Z2_syntax) <- seq_len(8)
# add actual stratum to this data frame
numprobs_Z2_syntax['Z2stratum_syntax'] <- lang_data_raw$Z2stratum_syntax
# numprob = probability for each person of being in their actual stratum as listed in the Zstratum column
lang_data_raw['numprob_Z2_syntax'] <- as.numeric(numprobs_Z2_syntax[cbind(seq_len(nrow(numprobs_Z2_syntax)), match(numprobs_Z2_syntax$Z2stratum_syntax, names(numprobs_Z2_syntax)))])


# divide numerator by denominator to get the weights
# qbw = quantile binning weight
lang_data_raw['Z2qbw_syntax'] <- lang_data_raw$numprob_Z2_syntax/lang_data_raw$denprob_Z2_syntax

# combined weights - missingness already incorporated via Z1 weights
lang_data_raw$combqbw_s_final <- lang_data_raw$Z1qbw_s_final * lang_data_raw$Z2qbw_syntax



# Now we run syntax outcome models

# Model with Z1 and Z2 as separate predictors (Table S2)

s_sep <- lm(c_celf_rs ~ Z1_syntax_rescale + Z2_syntax_rescale, weights = combqbw_s_final, lang_data_raw)

# Save model summary as s_m2
s_m2 <- robust(s_sep)
print(s_m2)

# confidence intervals

s_m2_z1_lowlim <- s_m2[2] - (qt(0.975, s_sep$df.residual) * s_m2[5])
s_m2_z1_uplim <- s_m2[2] + (qt(0.975, s_sep$df.residual) * s_m2[5])
print(s_m2_z1_lowlim)
print(s_m2_z1_uplim)

s_m2_z2_lowlim <- s_m2[3] - (qt(0.975, s_sep$df.residual) * s_m2[6])
s_m2_z2_uplim <- s_m2[3] + (qt(0.975, s_sep$df.residual) * s_m2[6])
print(s_m2_z2_lowlim)
print(s_m2_z2_uplim)


# AICc for this model
s_sep_AIC <- AICc(s_sep)

print(s_sep_AIC)


# hypothesis test: are effects of syntax Z1 and Z2 equal?
# model cumulative + Z2 and look at coefficient for Z2

s_diff_2 <- lm(c_celf_rs ~ cumulative_syntax + Z2_syntax_rescale, weights = combqbw_s_final, lang_data_raw)

print(robust(s_diff_2))


# calculate confidence interval for difference
coef_diff <- s_sep$coefficients[['Z2_syntax_rescale']] - s_sep$coefficients[['Z1_syntax_rescale']]
var_covar <- vcov(s_sep)
diff_se <- sqrt(var_covar[2,2] + var_covar[3,3] - 2 * var_covar[2,3])
lower_bound <- coef_diff - (diff_se * 2)
upper_bound <- coef_diff + (diff_se * 2)

print(lower_bound)
print(upper_bound)
