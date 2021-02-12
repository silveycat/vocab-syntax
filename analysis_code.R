# R script with analysis code reported in:
# 'Effects of time-varying parent input on child language outcomes differ for vocabulary and syntax'

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

# 1) IMPUTATION (Text S2) AND WEIGHT GENERATION (Text S3)

# read in raw data
lang_data_raw <- read.table("fake_lang_data.txt", header=TRUE)

# gender and birth order are factors
lang_data_raw$c_gen <- factor(lang_data_raw$c_gen)
lang_data_raw$c_bo <- factor(lang_data_raw$c_bo)

# display missing data pattern
md.pattern(lang_data_raw, rotate.names = TRUE)

# run multiple imputation using mice
# method: predictive mean matching

# first, use passive imputation to calculate variables that depend on others

# create a vector of missing values (NA) for each of the variables we need to calculate

# X1, composite measure of child language at 26 months (depends on c_wt_26 and c_mlu_26)
compos26 <- NA

# scaled versions of syntax Z1 and Z2 (depend on p_cps_14 and p_cps_30)
Z1_syntax_rescale <- NA
Z2_syntax_rescale <- NA

# Z1 and Z2 strata for vocab and syntax (depend on p_wt_14, p_wt_30, p_cps_14 and p_cps_30)

Z1stratum_vocab <- NA
Z2stratum_vocab <- NA

Z1stratum_syntax <- NA
Z2stratum_syntax <- NA

# cumulative vocab and syntax (Z1 + Z2)

cumulative_vocab <- NA
cumulative_syntax <- NA

# weights (these are calculated for each imputation using same formulas)

Z1qbw_vocab <- NA
Z2qbw_vocab <- NA

Z1qbw_syntax <- NA
Z2qbw_syntax <- NA

combqbw_vocab <- NA
combqbw_syntax <- NA

# add all these to the data frame ready for imputation

lang_data_plus <- cbind(lang_data_raw, compos26, Z1_syntax_rescale, Z2_syntax_rescale,
                        Z1stratum_vocab, Z2stratum_vocab, Z1stratum_syntax, Z2stratum_syntax,
                        cumulative_vocab, cumulative_syntax,
                        Z1qbw_vocab, Z2qbw_vocab, Z1qbw_syntax, Z2qbw_syntax,
                        combqbw_vocab, combqbw_syntax)



# function for calculating quantile binning weights according to Naimi et al.'s method
# Calculates weights using the same models for each imputed dataset
# takes the following arguments:
# dataset = dataset for which weights need to be calculated (R data frame)
# exposure = timepoint of exposure we are weighting for, Z1 or Z2 (string)
# variable = which type of parent input, vocab or syntax (string)
# covariates = list of covariates to include in weight-generating model (vector of strings)
# When exposure = Z1, the numerator of the weights is the marginal probability, 1/8
# When exposure = Z2, the function selects the appropriate Z1 based on variable,
# and calculates the numerator as the probability of Z2 given Z1

weight_calc <- function(dataset, exposure, variable, covariates) {
  # construct appropriate stratum column name e.g. Z1_stratum_vocab
  stratum <- paste0(exposure, "stratum")
  stratum <- paste(stratum, variable, sep = "_")
  # also need a version that is a factor, for including in formula below
  stratum_call <- paste0("as.factor(", stratum, ")")
  # convert covariate list into string with +s for including in formula below
  covar_string <- paste(covariates, collapse = " + ")
  # construct regression formula, e.g. as.factor(Z1_stratum_vocab) ~ p_viq + p_inc
  formula <- paste(stratum_call, covar_string, sep = " ~ ")

  # use this formula in ordinal regression to predict stratum assignment from covariates
  denominator_model <- polr(formula = formula, data = dataset)

  # save the probabilities of being in each stratum as a new data frame, predprobs
  denprobs <- data.frame(predict(denominator_model, type="p"))

  # save the probability for each person of being in their own stratum
  # first, make the column names match by assigning each stratum number to the appropriate column heading
  colnames(denprobs) <- seq_len(8)
  # add actual stratum to this data frame
  denprobs[stratum] <- dataset[stratum]

  # personprob = probability for each person of being in their actual stratum
  personprob <- as.numeric(denprobs[cbind(seq_len(nrow(denprobs)),
                                           match(denprobs[stratum][,1], names(denprobs)))])

  # calculate and return the weights
  # create appropriate numerator for Z1 vs. Z2 weights
  # if exposure = Z1:
  if (exposure == "Z1") {
    # numerator = marginal probability of being in each stratum, i.e. 1/8
    # denominator = person-specific probability of being in their stratum
    weights <- 1/8/(personprob)
  } else if (exposure == "Z2") {
    # for Z2, numerator = probability of being in Z2 stratum given Z1
    # The following code selects the appropriate Z1 variable,
    # depending if we're looking at vocab or syntax
    if (variable == "vocab") {
      Z1var <- "p_wt_14"
    } else if (variable == "syntax") {
      Z1var <- "p_cps_14"
    }
    # numerator = probability of being in Z2 stratum given Z1
    formula <- paste(stratum_call, Z1var, sep = " ~ ")

    # use this formula in ordinal regression to predict stratum assignment from Z1
    numerator_model <- polr(formula = formula, data = dataset)

    # save the probabilities of being in each stratum as a new data frame
    numprobs <- data.frame(predict(numerator_model, type="p"))
    # save the probability for each person of being in their own stratum
    # first, make the column names match by assigning each stratum number to the appropriate column heading
    colnames(numprobs) <- seq_len(8)
    # add actual stratum to this data frame
    numprobs[stratum] <- dataset[stratum]
    # numprob = probability for each person of being in their actual stratum as listed in the Zstratum column
    numprob <- as.numeric(numprobs[cbind(seq_len(nrow(numprobs)), match(numprobs[stratum][,1], names(numprobs)))])


    # denominator = person-specific probability of being in their stratum given Z1, X1 and X0
    weights <- numprob/personprob
  }
  return(weights)
}

# now we have the weight-generating function, we use this in our multiple imputation process

# set up template imputation without running to get method vector
ini <- mice(lang_data_plus, maxit = 0)

# inspect the prediction matrix (which missing variables are being predicted from which other variables)
pred <- ini$pred

# quickpred function automatically excludes predictors of <0.1 correlation
pred <- quickpred(lang_data_plus)

# save method vector (by default, predictive mean matching for all variables except factors)
meth <- ini$meth

# calculation method for compos26: average of z-scores for c_wt_26 and c_mlu_26
meth['compos26'] <- "~ I((scale(c_wt_26) + scale(c_mlu_26))/2)"

# calculation method for Z1_syntax_scaled and Z2_syntax_scaled: raw scores * 100
meth['Z1_syntax_rescale'] <- "~ I(p_cps_14 * 100)"
meth['Z2_syntax_rescale'] <- "~ I(p_cps_30 * 100)"

# calculation method for strata: 8 quantiles

meth['Z1stratum_vocab'] <- "~I(ntile(p_wt_14, 8))"
meth['Z2stratum_vocab'] <- "~I(ntile(p_wt_30, 8))"
meth['Z1stratum_syntax'] <- "~I(ntile(p_cps_14, 8))"
meth['Z2stratum_syntax'] <- "~I(ntile(p_cps_30, 8))"

# calculation method for cumulative Z1 and Z2: Z1 + Z2

meth['cumulative_vocab'] <- "~I(p_wt_14 + p_wt_30)"
meth['cumulative_syntax'] <- "~I(Z1_syntax_rescale + Z2_syntax_rescale)"

# calculation method for weights
# involves calling the weight_calc function defined above,
# with the appropriate exposure, variable, and covariates
# covariates were assessed for inclusion using one imputed dataset,
# via the method in section 3) below (reported in Text S6)

meth['Z1qbw_vocab'] <- "~ I(weight_calc(data, \"Z1\", \"vocab\", c(\"p_viq\", \"p_inc\", \"c_gen\")))"
meth['Z2qbw_vocab'] <- "~ I(weight_calc(data, \"Z2\", \"vocab\", c(\"p_wt_14\", \"compos26\")))"
meth['Z1qbw_syntax'] <- "~ I(weight_calc(data, \"Z1\", \"syntax\", c(\"p_viq\", \"c_bo\", \"c_wt_14\")))"
meth['Z2qbw_syntax'] <- "~ I(weight_calc(data, \"Z2\", \"syntax\", c(\"p_cps_14\",\"compos26\", \"p_viq\", \"c_gen\", \"p_educ\")))"

# calculation method for combined weights: multiply Z1 * Z2 weights

meth['combqbw_vocab'] <- "~ I(Z1qbw_vocab * Z2qbw_vocab)"
meth['combqbw_syntax'] <- "~ I(Z1qbw_syntax * Z2qbw_syntax)"

# now everything is set up, run imputation
lang_data_imp <- mice(lang_data_plus, meth=meth, pred = pred)

# this generates 5 complete datasets


# 2) OUTCOME MODELS

# A) VOCABULARY

# 'Control' model replicating previous work: Z1 predictor with W1 weights. Reported in text (p.9)

v_control <- with(lang_data_imp, lm(c_ppvt_74 ~ p_wt_14, weights = Z1qbw_vocab))

# pool results according to Rubin's rules
pool.v_control <- pool(v_control)

# save model summary as v_m1
v_m1 <- summary(pool.v_control)
print(v_m1)

# calculate confidence interval for effect of Z1

v_m1_z1_lowlim <- v_m1[2, 'estimate'] - (qt(0.975, v_m1[2, 'df']) * v_m1[2, 'std.error'])
v_m1_z1_uplim <- v_m1[2, 'estimate'] + (qt(0.975, v_m1[2, 'df']) * v_m1[2, 'std.error'])
print(v_m1_z1_lowlim)
print(v_m1_z1_uplim)

# Model with Z1 and Z2 as separate predictors (row 1 of Table 2)

v_sep <- with(lang_data_imp, lm(c_ppvt_74 ~ p_wt_14 + p_wt_30, weights = combqbw_vocab))

# pool results according to Rubin's rules
pool.v_sep <- pool(v_sep)

# save model summary as v_m2
v_m2 <- summary(pool.v_sep)
print(v_m2)

# confidence intervals

v_m2_z1_lowlim <- v_m2[2, 'estimate'] - (qt(0.975, v_m2[2, 'df']) * v_m2[2, 'std.error'])
v_m2_z1_uplim <- v_m2[2, 'estimate'] + (qt(0.975, v_m2[2, 'df']) * v_m2[2, 'std.error'])
print(v_m2_z1_lowlim)
print(v_m2_z1_uplim)

v_m2_z2_lowlim <- v_m2[3, 'estimate'] - (qt(0.975, v_m2[3, 'df']) * v_m2[3, 'std.error'])
v_m2_z2_uplim <- v_m2[3, 'estimate'] + (qt(0.975, v_m2[3, 'df']) * v_m2[3, 'std.error'])
print(v_m2_z2_lowlim)
print(v_m2_z2_uplim)

# AICc for this model
# averaged over imputed datasets
v_sep_AIC <- NULL
for (i in 1:5) {
  v_sep_AIC[i] <- AICc(v_sep$analyses[[i]])
}

print(mean(v_sep_AIC))

# hypothesis test: are effects of Z1 and Z2 equal?

# to answer this question, model cumulative + Z2 and look at coefficient for Z2

v_diff <- with(lang_data_imp, lm(c_ppvt_74 ~ cumulative_vocab + p_wt_30, weights = combqbw_vocab))

# pool results according to Rubin's rules
pool.v_diff <- pool(v_diff)
summary(pool.v_diff)

# calculate confidence intervals for difference and average across datasets
v_diff_lower_bounds <- NULL
v_diff_upper_bounds <- NULL
for (i in 1:5) {
  model_1 <- v_sep$analyses[[i]]
  coef_diff <- model_1$coefficients[['p_wt_30']] - model_1$coefficients[['p_wt_14']]
  var_covar <- vcov(model_1)
  diff_se <- sqrt(var_covar[2,2] + var_covar[3,3] - 2 * var_covar[2,3])
  lower_bound <- coef_diff - (diff_se * 2)
  upper_bound <- coef_diff + (diff_se * 2)
  v_diff_lower_bounds <- c(v_diff_lower_bounds, lower_bound)
  v_diff_upper_bounds <- c(v_diff_upper_bounds, upper_bound)
}
print(mean(v_diff_lower_bounds))
print(mean(v_diff_upper_bounds))

# Constant effects model (row 3 of Table 2)

v_const <- with(lang_data_imp, lm(c_ppvt_74 ~ cumulative_vocab, weights = combqbw_vocab))

# Pool results according to Rubin's rules
pool.v_const <- pool(v_const)

# save model summary as v_m3
v_m3 <- summary(pool.v_const)
print(v_m3)

# confidence interval

v_m3_cumul_lowlim <- v_m3[2, 'estimate'] - (qt(0.975, v_m3[2, 'df']) * v_m3[2, 'std.error'])
v_m3_cumul_uplim <- v_m3[2, 'estimate'] + (qt(0.975, v_m3[2, 'df']) * v_m3[2, 'std.error'])
print(v_m3_cumul_lowlim)
print(v_m3_cumul_uplim)


# AICc for this model
# averaged over imputed datasets
v_const_AIC <- NULL
for (i in 1:5) {
  v_const_AIC[i] <- AICc(v_const$analyses[[i]])
}

print(mean(v_const_AIC))

# First extreme model: Z1 > 0, Z2 = 0 (row 4 of Table 2)

v_Z1 <- with(lang_data_imp, lm(c_ppvt_74 ~ p_wt_14, weights = combqbw_vocab))

# Pool results according to Rubin's rules
pool.v_Z1 <- pool(v_Z1)

# Save model summary as v_m4
v_m4 <- summary(pool.v_Z1)
print(v_m4)

# confidence interval

v_m4_z1_lowlim <- v_m4[2, 'estimate'] - (qt(0.975, v_m4[2, 'df']) * v_m4[2, 'std.error'])
v_m4_z1_uplim <- v_m4[2, 'estimate'] + (qt(0.975, v_m4[2, 'df']) * v_m4[2, 'std.error'])
print(v_m4_z1_lowlim)
print(v_m4_z1_uplim)


# AICc for this model
# averaged over imputed datasets
v_Z1_AIC <- NULL
for (i in 1:5) {
  v_Z1_AIC[i] <- AICc(v_Z1$analyses[[i]])
}

print(mean(v_Z1_AIC))

# Second 'extreme' model: Z1 = 0, Z2 > 0 (row 5 of Table 2)

v_Z2 <- with(lang_data_imp, lm(c_ppvt_74 ~ p_wt_30, weights = combqbw_vocab))

# Pool results according to Rubin's rules
pool.v_Z2 <- pool(v_Z2)

# Save model summary as v_m5
v_m5 <- summary(pool.v_Z2)
print(v_m5)

# confidence interval

v_m5_z2_lowlim <- v_m5[2, 'estimate'] - (qt(0.975, v_m5[2, 'df']) * v_m5[2, 'std.error'])
v_m5_z2_uplim <- v_m5[2, 'estimate'] + (qt(0.975, v_m5[2, 'df']) * v_m5[2, 'std.error'])
print(v_m5_z2_lowlim)
print(v_m5_z2_uplim)


# AICc for this model
# averaged over imputed datasets
v_Z2_AIC <- NULL
for (i in 1:5) {
  v_Z2_AIC[i] <- AICc(v_Z2$analyses[[i]])
}

print(mean(v_Z2_AIC))

# B) SYNTAX

# 'Control' model replicating previous work: Z1 predictor with W1 weights. Reported in text (p.12)

s_control <- with(lang_data_imp, lm(c_celf_rs ~ Z1_syntax_rescale, weights = Z1qbw_syntax))

# Pool results according to Rubin's rules
pool.s_control <- pool(s_control)

# Save model summary as s_m1
s_m1 <- summary(pool.s_control)
print(s_m1)

# confidence interval

s_m1_z1_lowlim <- s_m1[2, 'estimate'] - (qt(0.975, s_m1[2, 'df']) * s_m1[2, 'std.error'])
s_m1_z1_uplim <- s_m1[2, 'estimate'] + (qt(0.975, s_m1[2, 'df']) * s_m1[2, 'std.error'])
print(s_m1_z1_lowlim)
print(s_m1_z1_uplim)


# Model with Z1 and Z2 as separate predictors (Table 4)

s_sep <- with(lang_data_imp, lm(c_celf_rs ~ Z1_syntax_rescale + Z2_syntax_rescale, weights = combqbw_syntax))

# Pool results according to Rubin's rules
pool.s_sep <- pool(s_sep)

# Save model summary as s_m2
s_m2 <- summary(pool.s_sep)
print(s_m2)

# confidence intervals

s_m2_z1_lowlim <- s_m2[2, 'estimate'] - (qt(0.975, s_m2[2, 'df']) * s_m2[2, 'std.error'])
s_m2_z1_uplim <- s_m2[2, 'estimate'] + (qt(0.975, s_m2[2, 'df']) * s_m2[2, 'std.error'])
print(s_m2_z1_lowlim)
print(s_m2_z1_uplim)

s_m2_z2_lowlim <- s_m2[3, 'estimate'] - (qt(0.975, s_m2[3, 'df']) * s_m2[3, 'std.error'])
s_m2_z2_uplim <- s_m2[3, 'estimate'] + (qt(0.975, s_m2[3, 'df']) * s_m2[3, 'std.error'])
print(s_m2_z2_lowlim)
print(s_m2_z2_uplim)


# AICc for this model
# averaged over imputed datasets
s_sep_AIC <- NULL
for (i in 1:5) {
  s_sep_AIC[i] <- AICc(s_sep$analyses[[i]])
}

print(mean(s_sep_AIC))


# hypothesis test: are effects of syntax Z1 and Z2 equal?
# model cumulative + Z2 and look at coefficient for Z2

s_diff_2 <- with(lang_data_imp, lm(c_celf_rs ~ cumulative_syntax + Z2_syntax_rescale, weights = combqbw_syntax))

# pool results according to Rubin's rules
pool.s_diff_2 <- pool(s_diff_2)
summary(pool.s_diff_2)


# calculate confints for difference and average across datasets
s_diff_lower_bounds <- NULL
s_diff_upper_bounds <- NULL
for (i in 1:5) {
  model_1 <- s_sep$analyses[[i]]
  coef_diff <- model_1$coefficients[['Z2_syntax_rescale']] - model_1$coefficients[['Z1_syntax_rescale']]
  var_covar <- vcov(model_1)
  diff_se <- sqrt(var_covar[2,2] + var_covar[3,3] - 2 * var_covar[2,3])
  lower_bound <- coef_diff - (diff_se * 2)
  upper_bound <- coef_diff + (diff_se * 2)
  s_diff_lower_bounds <- c(s_diff_lower_bounds, lower_bound)
  s_diff_upper_bounds <- c(s_diff_upper_bounds, upper_bound)
}

print(mean(s_diff_lower_bounds))
print(mean(s_diff_upper_bounds))



# 3) COVARIATE SELECTION, BALANCE, AND COMMON SUPPORT
# (Text S6, Tables S3-S6, Text S7, Figures S3-S6)

# select one complete dataset
lang_data_balance <- complete(lang_data_imp)

# for the purposes of balance checking, we use continuous versions of gender and birth order
# this makes it easier to compare t-values and standardised regression coefficients
lang_data_balance$c_gen_cont <- as.numeric(lang_data_balance$c_gen)
lang_data_balance$c_bo_cont <- as.numeric(lang_data_balance$c_bo)

# I. Covariate selection and balance checking for Z1 vocabulary input (Table S3)

# Check balance of covariates X0 on Z1

# stratify the sample on Z1 (parent word types at 14 months) into 8 strata
lang_data_balance['Z1stratum_vocab'] <- factor(ntile(lang_data_balance$p_wt_14, 8))

# Assess covariates for inclusion
# First, run all unweighted regressions below predicting each covariate from Z1

# Add covariate where Z1 has largest t value (where t > 1.67) to ordinal model
# Generate weights
# Run all weighted regressions predicting each covariate X0 from Z1
# If t-value for Z1 coefficient in any model is > 1.67,
# add covariate where Z1 has largest t value to ordinal model
# Repeat until balance achieved on all covariates
# (i.e. t-value for Z1 coefficient is < 1.67 > -1.67 in all weighted regressions)

# In the final balance table, report standardised beta and t-value for each covariate,
# before and after weighting

# Covariates justified for inclusion for Z1 vocab were parent verbal IQ, household income, & gender

# predict stratum assignment from covariates X0
denominator_Z1_vocab <- polr(Z1stratum_vocab ~ p_viq + p_inc + c_gen, lang_data_balance)

# save the probabilities of being in each stratum as a new data frame
denprobs_Z1_vocab <- data.frame(predict(denominator_Z1_vocab, type="p"))

# save the probability for each person of being in their own stratum
# first, make the column names match by assigning each stratum number to the appropriate column heading
colnames(denprobs_Z1_vocab) <- seq_len(8)
# add veridical stratum to this data frame
denprobs_Z1_vocab['Z1stratum_vocab'] <- lang_data_balance$Z1stratum_vocab
# denprob = probability for each person of being in their actual stratum as listed in the Zstratum column
lang_data_balance['denprob_Z1_vocab'] <- as.numeric(denprobs_Z1_vocab[cbind(seq_len(nrow(denprobs_Z1_vocab)), match(denprobs_Z1_vocab$Z1stratum_vocab, names(denprobs_Z1_vocab)))])

# calculate the weights
# here, the numerator = marginal probability of being in each stratum, i.e. 1/8
# denominator = person-specific probability of being in their stratum
# qbw = quantile binning weight

lang_data_balance['Z1qbw_vocab'] <- 1/8/(lang_data_balance$denprob_Z1_vocab)

# balance checking regressions

# birth order
# unweighted
c_bo_unweighted <- lm(c_bo_cont ~ p_wt_14, lang_data_balance)
print(robust(c_bo_unweighted))
print(lm.beta(c_bo_unweighted))

# weighted
c_bo_weighted <- lm(c_bo_cont ~ p_wt_14, weights = Z1qbw_vocab, data = lang_data_balance)
print(robust(c_bo_weighted))
print(lm.beta(c_bo_weighted))

# child gender
# unweighted
c_gen_unweighted <- lm(c_gen_cont ~ p_wt_14, lang_data_balance)
print(robust(c_gen_unweighted))
print(lm.beta(c_gen_unweighted))

# weighted
c_gen_weighted <- lm(c_gen_cont ~ p_wt_14, weights = Z1qbw_vocab, data = lang_data_balance)
print(robust(c_gen_weighted))
print(lm.beta(c_gen_weighted))

# child word types 14 months
# unweighted
c_wt_14_unweighted <- lm(c_wt_14 ~ p_wt_14, lang_data_balance)
print(robust(c_wt_14_unweighted))
print(lm.beta(c_wt_14_unweighted))

# weighted
c_wt_14_weighted <- lm(c_wt_14 ~ p_wt_14, weights = Z1qbw_vocab, data = lang_data_balance)
print(robust(c_wt_14_weighted))
print(lm.beta(c_wt_14_weighted))

# child gesture types 14 months
# unweighted
c_gt_14_unweighted <- lm(c_gt_14 ~ p_wt_14, lang_data_balance)
print(robust(c_gt_14_unweighted))
print(lm.beta(c_gt_14_unweighted))

# weighted
c_gt_14_weighted <- lm(c_gt_14 ~ p_wt_14, weights = Z1qbw_vocab, data = lang_data_balance)
print(robust(c_gt_14_weighted))
print(lm.beta(c_gt_14_weighted))

# parent IQ
# unweighted
p_viq_unweighted <- lm(p_viq ~ p_wt_14, lang_data_balance)
print(robust(p_viq_unweighted))
print(lm.beta(p_viq_unweighted))

# weighted
p_viq_weighted <- lm(p_viq ~ p_wt_14, weights = Z1qbw_vocab, data = lang_data_balance)
print(robust(p_viq_weighted))
print(lm.beta(p_viq_weighted))

# parent education
# unweighted
p_educ_unweighted <- lm(p_educ ~ p_wt_14, lang_data_balance)
print(robust(p_educ_unweighted))
print(lm.beta(p_educ_unweighted))

# weighted
p_educ_weighted <- lm(p_educ ~ p_wt_14, weights = Z1qbw_vocab, data = lang_data_balance)
print(robust(p_educ_weighted))
print(lm.beta(p_educ_weighted))

# parent income
# unweighted
p_inc_unweighted <- lm(p_inc ~ p_wt_14, lang_data_balance)
print(robust(p_inc_unweighted))
print(lm.beta(p_inc_unweighted))

# weighted
p_inc_weighted <- lm(p_inc ~ p_wt_14, weights = Z1qbw_vocab, data = lang_data_balance)
print(robust(p_inc_weighted))
print(lm.beta(p_inc_weighted))


# II. Common support for vocabulary Z1 (Figure S3)

# use denprobs frame to figure out which was highest probability stratum for each case
denprobs_Z1_vocab$max_strat <- colnames(denprobs_Z1_vocab[,-9])[max.col(denprobs_Z1_vocab[,-9],ties.method="first")]

# convert these to 4 categories (1 per 2 strata) and then show histograms of which stratum actually in

denprobs_Z1_vocab$category <- NA

# if max_strat is 1 or 2, designate category as 1
denprobs_Z1_vocab$category[denprobs_Z1_vocab$max_strat == 1 | denprobs_Z1_vocab$max_strat == 2] <- 1

# if max_strat is 3 or 4, designate category as 2
denprobs_Z1_vocab$category[denprobs_Z1_vocab$max_strat == 3 | denprobs_Z1_vocab$max_strat == 4] <- 2

# if max_strat is 5 or 6, designate category as 3
denprobs_Z1_vocab$category[denprobs_Z1_vocab$max_strat == 5 | denprobs_Z1_vocab$max_strat == 6] <- 3

# if max_strat is 7 or 8, designate category as 4
denprobs_Z1_vocab$category[denprobs_Z1_vocab$max_strat == 7 | denprobs_Z1_vocab$max_strat == 8] <- 4

# plot: histogram faceted by category, with actual stratum as x

# New facet label names
facet.labs <- c("quantile 1-2 predicted", "quantile 3-4 predicted", "quantile 5-6 predicted",
                "quantile 7-8 predicted")
names(facet.labs) <- c("1", "2", "3", "4")

Z1_cs_plot_vocab <- ggplot(denprobs_Z1_vocab, aes(x = Z1stratum_vocab, group = category))
Z1_cs_plot_vocab + geom_histogram(stat = "count") + facet_grid(. ~ category,
                                                         labeller = labeller(category = facet.labs)) +
  theme_bw() + scale_x_discrete("\nQuantile observed") + scale_y_continuous("Frequency\n")

ggsave("Common support Z1 vocab final.png")


# III. Covariate selection for Z2 vocabulary input (Table S4)

# first we stratify the sample on Z2 (parent word types at 30 months) into 8 strata
lang_data_balance['Z2stratum_vocab'] <- factor(ntile(lang_data_balance$p_wt_30, 8))

# For Z2, we run 'combined' balance checks, using combined weights.
# To get combined weights, weneed to generate Z2 weights first.

# first, we calculate the denominator:
# the person-specific probability of being in their Z2 stratum given Z1, X1 and X0
# At the beginning of the balance check, we only include Z1 in the denominator model.
# We will only add covariates to this model if it turns out they are not balanced with respect to Z2.
# (in our balance checks we found X1 was not balanced, which is why it is in the model below)

# predict stratum assignment from Z1 and covariates
denominator_Z2_vocab <- polr(Z2stratum_vocab ~ p_wt_14 + compos26, lang_data_balance)

# save the probabilities of being in each stratum as a new data frame
denprobs_Z2_vocab <- data.frame(predict(denominator_Z2_vocab, type="p"))

# save the probability for each person of being in their own stratum
# first, make the column names match by assigning each stratum number to the appropriate column heading
colnames(denprobs_Z2_vocab) <- seq_len(8)
# add actual stratum to this data frame
denprobs_Z2_vocab['Z2stratum_vocab'] <- lang_data_balance$Z2stratum_vocab
# personprob = probability for each person of being in their actual stratum as listed in the Zstratum column
lang_data_balance['denprob_Z2_vocab'] <- as.numeric(denprobs_Z2_vocab[cbind(seq_len(nrow(denprobs_Z2_vocab)), match(denprobs_Z2_vocab$Z2stratum_vocab, names(denprobs_Z2_vocab)))])

# Now, we calculate the numerator of the Z2 weights
# numerator = probability of being in Z2 stratum given Z1

numerator_Z2_vocab <- polr(Z2stratum_vocab ~ p_wt_14, lang_data_balance)
# save the probabilities of being in each stratum as a new data frame
numprobs_Z2_vocab <- data.frame(predict(numerator_Z2_vocab, type="p"))
# save the probability for each person of being in their own stratum
# first, make the column names match by assigning each stratum number to the appropriate column heading
colnames(numprobs_Z2_vocab) <- seq_len(8)
# add actual stratum to this data frame
numprobs_Z2_vocab['Z2stratum_vocab'] <- lang_data_balance$Z2stratum_vocab
# numprob = probability for each person of being in their actual stratum as listed in the Zstratum column
lang_data_balance['numprob_Z2_vocab'] <- as.numeric(numprobs_Z2_vocab[cbind(seq_len(nrow(numprobs_Z2_vocab)), match(numprobs_Z2_vocab$Z2stratum_vocab, names(numprobs_Z2_vocab)))])


# divide the numerator probability by the denominator probability to get the weight
# Note that at the beginning of the balance checking process, the numerator and denominator are the same
# so Z2 weights are initially set to 1.
# This will change if & when we add covariates to the denominator model.

lang_data_balance['Z2qbw_vocab_check'] <- lang_data_balance$numprob_Z2_vocab/lang_data_balance$denprob_Z2_vocab


# Multiply Z1 and Z2 weights to create combined weights

lang_data_balance['combqbw_vocab'] <- lang_data_balance$Z1qbw_vocab * lang_data_balance$Z2qbw_vocab

# Balance checking regressions, using combined weights
# Run the weighted regressions to check balance.
# The unweighted versions are there so we can report associations before and after weighting.

# Check balance on X0 and X1. Here, both Z1 and Z2 are predictors in the balance checking models;
# however, we only need to pay attention to and report the coefficient for Z2.
# If t > 1.67 or t < -1.67, balance is not achieved.
# If this happens, add the covariate with the highest t-value to our denominator model for Z2,
# recalculate the Z2 and combined weights,
# and run the weighted balance checking regressions again.
# Repeat these steps until the t for Z2 is < 1.67 > -1.67 in all the weighted models.

# X1 - child language at 26 months
# unweighted - don't use this to check balance, only to report associations before weighting
X1_comb_unweighted <- lm(compos26 ~ p_wt_14 + p_wt_30, lang_data_balance)
print(robust(X1_comb_unweighted))
print(lm.beta(X1_comb_unweighted))

# weighted - use this to check balance
X1_comb_weighted <- lm(compos26 ~ p_wt_14 + p_wt_30, weights = combqbw_vocab, lang_data_balance)
print(robust(X1_comb_weighted))
print(lm.beta(X1_comb_weighted))

# birth order
# unweighted
c_bo_comb_unweighted <- lm(c_bo_cont ~ p_wt_14 + p_wt_30, lang_data_balance)
print(robust(c_bo_comb_unweighted))
print(lm.beta(c_bo_comb_unweighted))

# weighted
c_bo_comb_weighted <- lm(c_bo_cont ~ p_wt_14 + p_wt_30, weights = combqbw_vocab, lang_data_balance)
print(robust(c_bo_comb_weighted))
print(lm.beta(c_bo_comb_weighted))

# gender
# unweighted
c_gen_comb_unweighted <- lm(c_gen_cont ~ p_wt_14 + p_wt_30, lang_data_balance)
print(robust(c_gen_comb_unweighted))
print(lm.beta(c_gen_comb_unweighted))

# weighted
c_gen_comb_weighted <- lm(c_gen_cont ~ p_wt_14 + p_wt_30, weights = combqbw_vocab, lang_data_balance)
print(robust(c_gen_comb_weighted))
print(lm.beta(c_gen_comb_weighted))

# child word types at 14 m
# unweighted
c_wt_14_comb_unweighted <- lm(c_wt_14 ~ p_wt_14 + p_wt_30, lang_data_balance)
print(robust(c_wt_14_comb_unweighted))
print(lm.beta(c_wt_14_comb_unweighted))

# weighted
c_wt_14_comb_weighted <- lm(c_wt_14 ~ p_wt_14 + p_wt_30, weights = combqbw_vocab, lang_data_balance)
print(robust(c_wt_14_comb_weighted))
print(lm.beta(c_wt_14_comb_weighted))

# child gesture types at 14 m
# unweighted
c_gt_14_comb_unweighted <- lm(c_gt_14 ~ p_wt_14 + p_wt_30, lang_data_balance)
print(robust(c_gt_14_comb_unweighted))
print(lm.beta(c_gt_14_comb_unweighted))

# weighted
c_gt_14_comb_weighted <- lm(c_gt_14 ~ p_wt_14 + p_wt_30, weights = combqbw_vocab, lang_data_balance)
print(robust(c_gt_14_comb_weighted))
print(lm.beta(c_gt_14_comb_weighted))

# parent IQ
# unweighted
p_viq_comb_unweighted <- lm(p_viq ~ p_wt_14 + p_wt_30, lang_data_balance)
print(robust(p_viq_comb_unweighted))
print(lm.beta(p_viq_comb_unweighted))

# weighted
p_viq_comb_weighted <- lm(p_viq ~ p_wt_14 + p_wt_30, weights = combqbw_vocab, lang_data_balance)
print(robust(p_viq_comb_weighted))
print(lm.beta(p_viq_comb_weighted))

# parent education
# unweighted
p_educ_comb_unweighted <- lm(p_educ ~ p_wt_14 + p_wt_30, lang_data_balance)
print(robust(p_educ_comb_unweighted))
print(lm.beta(p_educ_comb_unweighted))

# weighted
p_educ_comb_weighted <- lm(p_educ ~ p_wt_14 + p_wt_30, weights = combqbw_vocab, lang_data_balance)
print(robust(p_educ_comb_weighted))
print(lm.beta(p_educ_comb_weighted))

# parent income
# unweighted
p_inc_unweighted <- lm(p_inc ~ p_wt_14 + p_wt_30, lang_data_balance)
print(robust(p_inc_unweighted))
print(lm.beta(p_inc_unweighted))

# weighted
p_inc_weighted <- lm(p_inc ~ p_wt_14 + p_wt_30, weights = combqbw_vocab, lang_data_balance)
print(robust(p_inc_weighted))
print(lm.beta(p_inc_weighted))



# IV. Common support for vocabulary Z2 (Figure S4)

# use denprobs frame to figure out which was highest probability stratum for each case
denprobs_Z2_vocab$max_strat <- colnames(denprobs_Z2_vocab[,-9])[max.col(denprobs_Z2_vocab[,-9],ties.method="first")]

# convert these to 4 categories (1 per 2 strata) and then show histograms of which stratum actually in

denprobs_Z2_vocab$category <- NA

# if max_strat is 1 or 2, designate category as 1
denprobs_Z2_vocab$category[denprobs_Z2_vocab$max_strat == 1 | denprobs_Z2_vocab$max_strat == 2] <- 1

# if max_strat is 3 or 4, designate category as 2
denprobs_Z2_vocab$category[denprobs_Z2_vocab$max_strat == 3 | denprobs_Z2_vocab$max_strat == 4] <- 2

# if max_strat is 5 or 6, designate category as 3
denprobs_Z2_vocab$category[denprobs_Z2_vocab$max_strat == 5 | denprobs_Z2_vocab$max_strat == 6] <- 3

# if max_strat is 7 or 8, designate category as 4
denprobs_Z2_vocab$category[denprobs_Z2_vocab$max_strat == 7 | denprobs_Z2_vocab$max_strat == 8] <- 4

# plot: histogram faceted by category, with actual stratum as x

# New facet label names
facet.labs <- c("quantile 1-2 predicted", "quantile 3-4 predicted", "quantile 5-6 predicted",
                "quantile 7-8 predicted")
names(facet.labs) <- c("1", "2", "3", "4")

Z2_cs_plot_vocab <- ggplot(denprobs_Z2_vocab, aes(x = Z2stratum_vocab, group = category))
Z2_cs_plot_vocab + geom_histogram(stat = "count") + facet_grid(. ~ category,
                                                         labeller = labeller(category = facet.labs)) +
  theme_bw() + scale_x_discrete("\nQuantile observed") + scale_y_continuous("Frequency\n")

ggsave("Common support Z2 vocab final.png")


# V. Covariate selection and balance checking for Z1 syntax input (Table S5)

# Z1: parent average clauses per sentence at 14 months

# same approach as above

# stratify the sample on Z1 (parent cps at 14 months) into 8 strata
lang_data_balance['Z1stratum_syntax'] <- factor(ntile(lang_data_balance$p_cps_14, 8))

# assess covariates for inclusion
# same process as above
# justified for inclusion: p_viq, c_wt_14, c_bo

# predict stratum assignment from these covariates
denominator_Z1_syntax <- polr(Z1stratum_syntax ~ p_viq + c_bo + c_wt_14, lang_data_balance)

# save the probabilities of being in each stratum as a new data frame
denprobs_Z1_syntax <- data.frame(predict(denominator_Z1_syntax, type="p"))

# save the probability for each person of being in their own stratum
# first, make the column names match by assigning each stratum number to the appropriate column heading
colnames(denprobs_Z1_syntax) <- seq_len(8)
# add actual stratum to this data frame
denprobs_Z1_syntax['Z1stratum_syntax'] <- lang_data_balance$Z1stratum_syntax
# personprob = probability for each person of being in their actual stratum as listed in the Zstratum column
lang_data_balance['denprob_Z1_syntax'] <- as.numeric(denprobs_Z1_syntax[cbind(seq_len(nrow(denprobs_Z1_syntax)), match(denprobs_Z1_syntax$Z1stratum_syntax, names(denprobs_Z1_syntax)))])

# calculate the weights
# here, numerator = marginal probability of being in each stratum, i.e. 1/8
# denominator = person-specific probability of being in their stratum given covariates
# qbw = quantile binning weight
lang_data_balance['Z1qbw_syntax_check'] <- 1/8/(lang_data_balance$denprob_Z1_syntax)

# balance checking regressions

# birth order
# unweighted
c_bo_unweighted <- lm(c_bo_cont ~ p_cps_14, lang_data_balance)
print(robust(c_bo_unweighted))
print(lm.beta(c_bo_unweighted))

# weighted
c_bo_weighted <- lm(c_bo_cont ~ p_cps_14, weights = Z1qbw_syntax, data = lang_data_balance)
print(robust(c_bo_weighted))
print(lm.beta(c_bo_weighted))

# child gender
# unweighted
c_gen_unweighted <- lm(c_gen_cont ~ p_cps_14, lang_data_balance)
print(robust(c_gen_unweighted))
print(lm.beta(c_gen_unweighted))

# weighted
c_gen_weighted <- lm(c_gen_cont ~ p_cps_14, weights = Z1qbw_syntax, data = lang_data_balance)
print(robust(c_gen_weighted))
print(lm.beta(c_gen_weighted))

# child word types 14 months
# unweighted
c_wt_14_unweighted <- lm(c_wt_14 ~ p_cps_14, lang_data_balance)
print(robust(c_wt_14_unweighted))
print(lm.beta(c_wt_14_unweighted))

# weighted
c_wt_14_weighted <- lm(c_wt_14 ~ p_cps_14, weights = Z1qbw_syntax, data = lang_data_balance)
print(robust(c_wt_14_weighted))
print(lm.beta(c_wt_14_weighted))

# child gesture types 14 months
# unweighted
c_gt_14_unweighted <- lm(c_gt_14 ~ p_cps_14, lang_data_balance)
print(robust(c_gt_14_unweighted))
print(lm.beta(c_gt_14_unweighted))

# weighted
c_gt_14_weighted <- lm(c_gt_14 ~ p_cps_14, weights = Z1qbw_syntax, data = lang_data_balance)
print(robust(c_gt_14_weighted))
print(lm.beta(c_gt_14_weighted))

# parent IQ
# unweighted
p_viq_unweighted <- lm(p_viq ~ p_cps_14, lang_data_balance)
print(robust(p_viq_unweighted))
print(lm.beta(p_viq_unweighted))

# weighted
p_viq_weighted <- lm(p_viq ~ p_cps_14, weights = Z1qbw_syntax, data = lang_data_balance)
print(robust(p_viq_weighted))
print(lm.beta(p_viq_weighted))

# parent education
# unweighted
p_educ_unweighted <- lm(p_educ ~ p_cps_14, lang_data_balance)
print(robust(p_educ_unweighted))
print(lm.beta(p_educ_unweighted))

# weighted
p_educ_weighted <- lm(p_educ ~ p_cps_14, weights = Z1qbw_syntax, data = lang_data_balance)
print(robust(p_educ_weighted))
print(lm.beta(p_educ_weighted))

# parent income
# unweighted
p_inc_unweighted <- lm(p_inc ~ p_cps_14, lang_data_balance)
print(robust(p_inc_unweighted))
print(lm.beta(p_inc_unweighted))

# weighted
p_inc_weighted <- lm(p_inc ~ p_cps_14, weights = Z1qbw_syntax, data = lang_data_balance)
print(robust(p_inc_weighted))
print(lm.beta(p_inc_weighted))


# VI. Common support for syntax Z1 (Figure S5)

# use denprobs frame to figure out which was highest probability stratum for each case
denprobs_Z1_syntax$max_strat <- colnames(denprobs_Z1_syntax[,-9])[max.col(denprobs_Z1_syntax[,-9],ties.method="first")]

# convert these to 4 categories (1 per 2 strata) and then show histograms of which stratum actually in

denprobs_Z1_syntax$category <- NA

# if max_strat is 1 or 2, designate category as 1
denprobs_Z1_syntax$category[denprobs_Z1_syntax$max_strat == 1 | denprobs_Z1_syntax$max_strat == 2] <- 1

# if max_strat is 3 or 4, designate category as 2
denprobs_Z1_syntax$category[denprobs_Z1_syntax$max_strat == 3 | denprobs_Z1_syntax$max_strat == 4] <- 2

# if max_strat is 5 or 6, designate category as 3
denprobs_Z1_syntax$category[denprobs_Z1_syntax$max_strat == 5 | denprobs_Z1_syntax$max_strat == 6] <- 3

# if max_strat is 7 or 8, designate category as 4
denprobs_Z1_syntax$category[denprobs_Z1_syntax$max_strat == 7 | denprobs_Z1_syntax$max_strat == 8] <- 4

# plot: histogram faceted by category, with actual stratum as x

# New facet label names
facet.labs <- c("quantile 1-2 predicted", "quantile 3-4 predicted", "quantile 5-6 predicted",
                "quantile 7-8 predicted")
names(facet.labs) <- c("1", "2", "3", "4")

Z1_cs_plot_syntax <- ggplot(denprobs_Z1_syntax, aes(x = Z1stratum_syntax, group = category))
Z1_cs_plot_syntax + geom_histogram(stat = "count") + facet_grid(. ~ category,
                                                         labeller = labeller(category = facet.labs)) +
  theme_bw() + scale_x_discrete("\nQuantile observed") + scale_y_continuous("Frequency\n")

ggsave("Common support Z1 syntax final.png")


# VII. Covariate selection and balance checking for Z2 syntax input (Table S6)

# first we stratify the sample on Z2 (parent cps at 30 months) into 8 strata
lang_data_balance['Z2stratum_syntax'] <- factor(ntile(lang_data_balance$p_cps_30, 8))


# covariates for inclusion

# We use the same process as we used for vocabulary Z2
# denominator = person-specific probability of being in their stratum given Z1, X1 and X0
# initially we only put Z1 in the denominator model, so initially Z2 weights = 1

# our balance checks lead to p_viq, compos26 (X1), p_educ, and gender being included
denominator_Z2_syntax <- polr(Z2stratum_syntax ~ p_cps_14 + p_viq + compos26 + c_gen + p_educ, lang_data_balance)

# save the probabilities of being in each stratum as a new data frame
denprobs_Z2_syntax <- data.frame(predict(denominator_Z2_syntax, type="p"))

# save the probability for each person of being in their own stratum
# first, make the column names match by assigning each stratum number to the appropriate column heading
colnames(denprobs_Z2_syntax) <- seq_len(8)
# add actual stratum to this data frame
denprobs_Z2_syntax['Z2stratum_syntax'] <- lang_data_balance$Z2stratum_syntax
# personprob = probability for each person of being in their actual stratum as listed in the Zstratum column
lang_data_balance['denprob_Z2_syntax'] <- as.numeric(denprobs_Z2_syntax[cbind(seq_len(nrow(denprobs_Z2_syntax)), match(denprobs_Z2_syntax$Z2stratum_syntax, names(denprobs_Z2_syntax)))])

# calculate the weights
# as for vocab Z2, we now need a numerator model
# numerator = probability of being in each Z2 stratum given Z1
numerator_Z2_syntax <- polr(Z2stratum_syntax ~ p_cps_14, lang_data_balance)
# save the probabilities of being in each stratum as a new data frame
numprobs_Z2_syntax <- data.frame(predict(numerator_Z2_syntax, type="p"))

# save the probability for each person of being in their own stratum
# first, make the column names match by assigning each stratum number to the appropriate column heading
colnames(numprobs_Z2_syntax) <- seq_len(8)
# add actual stratum to this data frame
numprobs_Z2_syntax['Z2stratum_syntax'] <- lang_data_balance$Z2stratum_syntax
# numprob = probability for each person of being in their actual stratum as listed in the Zstratum column
lang_data_balance['numprob_Z2_syntax'] <- as.numeric(numprobs_Z2_syntax[cbind(seq_len(nrow(numprobs_Z2_syntax)), match(numprobs_Z2_syntax$Z2stratum_syntax, names(numprobs_Z2_syntax)))])


# divide numerator by denominator to get the weights
# qbw = quantile binning weight
lang_data_balance['Z2qbw_syntax_check'] <- lang_data_balance$numprob_Z2_syntax/lang_data_balance$denprob_Z2_syntax

# Multiply weights

lang_data_balance['combqbw_syntax'] <- lang_data_balance$Z1qbw_syntax * lang_data_balance$Z2qbw_syntax


# balance checking regressions
# Again, we run combined models with both Z1 and Z2 as predictors, using combined weights
# However, we only pay attention to the coefficient for Z2 (here, p_cps_30)
# If in any model the t value for this coefficient > 1.67 or < -1.67, add the covariate with the highest
# t-value to the denominator model above, recalculate Z2 and combined weights, and repeat

# X1
# unweighted
X1_comb_unweighted <- lm(compos26 ~ p_cps_14 + p_cps_30, lang_data_balance)
print(robust(X1_comb_unweighted))
print(lm.beta(X1_comb_unweighted))

# weighted
X1_comb_weighted <- lm(compos26 ~ p_cps_14 + p_cps_30, weights = combqbw_syntax, lang_data_balance)
print(robust(X1_comb_weighted))
print(lm.beta(X1_comb_weighted))

# birth order
# unweighted
c_bo_comb_unweighted <- lm(c_bo_cont ~ p_cps_14 + p_cps_30, lang_data_balance)
print(robust(c_bo_comb_unweighted))
print(lm.beta(c_bo_comb_unweighted))

# weighted
c_bo_comb_weighted <- lm(c_bo_cont ~ p_cps_14 + p_cps_30, weights = combqbw_syntax, lang_data_balance)
print(robust(c_bo_comb_weighted))
print(lm.beta(c_bo_comb_weighted))

# gender
# unweighted
c_gen_comb_unweighted <- lm(c_gen_cont ~ p_cps_14 + p_cps_30, lang_data_balance)
print(robust(c_gen_comb_unweighted))
print(lm.beta(c_gen_comb_unweighted))

# weighted
c_gen_comb_weighted <- lm(c_gen_cont ~ p_cps_14 + p_cps_30, weights = combqbw_syntax, lang_data_balance)
print(robust(c_gen_comb_weighted))
print(lm.beta(c_gen_comb_weighted))

# child word types at 14 m
# unweighted
c_wt_14_comb_unweighted <- lm(c_wt_14 ~ p_cps_14 + p_cps_30, lang_data_balance)
print(robust(c_wt_14_comb_unweighted))
print(lm.beta(c_wt_14_comb_unweighted))

# weighted
c_wt_14_comb_weighted <- lm(c_wt_14 ~ p_cps_14 + p_cps_30, weights = combqbw_syntax, lang_data_balance)
print(robust(c_wt_14_comb_weighted))
print(lm.beta(c_wt_14_comb_weighted))

# child gesture types at 14 m
# unweighted
c_gt_14_comb_unweighted <- lm(c_gt_14 ~ p_cps_14 + p_cps_30, lang_data_balance)
print(robust(c_gt_14_comb_unweighted))
print(lm.beta(c_gt_14_comb_unweighted))

# weighted
c_gt_14_comb_weighted <- lm(c_gt_14 ~ p_cps_14 + p_cps_30, weights = combqbw_syntax, lang_data_balance)
print(robust(c_gt_14_comb_weighted))
print(lm.beta(c_gt_14_comb_weighted))

# parent IQ
# unweighted
p_viq_comb_unweighted <- lm(p_viq ~ p_cps_14 + p_cps_30, lang_data_balance)
print(robust(p_viq_comb_unweighted))
print(lm.beta(p_viq_comb_unweighted))

# weighted
p_viq_comb_weighted <- lm(p_viq ~ p_cps_14 + p_cps_30, weights = combqbw_syntax, lang_data_balance)
print(robust(p_viq_comb_weighted))
print(lm.beta(p_viq_comb_weighted))

# parent education
# unweighted
p_educ_comb_unweighted <- lm(p_educ ~ p_cps_14 + p_cps_30, lang_data_balance)
print(robust(p_educ_comb_unweighted))
print(lm.beta(p_educ_comb_unweighted))

# weighted
p_educ_comb_weighted <- lm(p_educ ~ p_cps_14 + p_cps_30, weights = combqbw_syntax, lang_data_balance)
print(robust(p_educ_comb_weighted))
print(lm.beta(p_educ_comb_weighted))

# parent income
# unweighted
p_inc_comb_unweighted <- lm(p_inc ~ p_cps_14 + p_cps_30, lang_data_balance)
print(robust(p_inc_comb_unweighted))
print(lm.beta(p_inc_comb_unweighted))

# weighted
p_inc_comb_weighted <- lm(p_inc ~ p_cps_14 + p_cps_30, weights = combqbw_syntax, lang_data_balance)
print(robust(p_inc_comb_weighted))
print(lm.beta(p_inc_comb_weighted))


# VIII. Common support for Z2 syntax (Figure S6)

# use predprobs frame to figure out which was highest probability stratum for each case
denprobs_Z2_syntax$max_strat <- colnames(denprobs_Z2_syntax[,-9])[max.col(denprobs_Z2_syntax[,-9],ties.method="first")]

# convert these to 4 categories (1 per 2 strata) and then show histograms of which stratum actually in

denprobs_Z2_syntax$category <- NA

# if max_strat is 1 or 2, designate category as 1
denprobs_Z2_syntax$category[denprobs_Z2_syntax$max_strat == 1 | denprobs_Z2_syntax$max_strat == 2] <- 1

# if max_strat is 3 or 4, designate category as 2
denprobs_Z2_syntax$category[denprobs_Z2_syntax$max_strat == 3 | denprobs_Z2_syntax$max_strat == 4] <- 2

# if max_strat is 5 or 6, designate category as 3
denprobs_Z2_syntax$category[denprobs_Z2_syntax$max_strat == 5 | denprobs_Z2_syntax$max_strat == 6] <- 3

# if max_strat is 7 or 8, designate category as 4
denprobs_Z2_syntax$category[denprobs_Z2_syntax$max_strat == 7 | denprobs_Z2_syntax$max_strat == 8] <- 4

# plot: histogram faceted by category, with actual stratum as x

# New facet label names
facet.labs <- c("quantile 1-2 predicted", "quantile 3-4 predicted", "quantile 5-6 predicted",
                "quantile 7-8 predicted")
names(facet.labs) <- c("1", "2", "3", "4")

Z2_cs_plot_syntax <- ggplot(denprobs_Z2_syntax, aes(x = Z2stratum_syntax, group = category))
Z2_cs_plot_syntax + geom_histogram(stat = "count") + facet_grid(. ~ category,
                                                                labeller = labeller(category = facet.labs)) +
  theme_bw() + scale_x_discrete("\nQuantile observed") + scale_y_continuous("Frequency\n")

ggsave("Common support Z2 syntax final.png")


# 4) PREDICTED OUTCOME PLOTS (Figure 1)

# use raw, complete data (i.e. lang_data_raw)

# calculate cumulative vocab input
lang_data_raw$cumulative_vocab <- lang_data_raw$p_wt_14 + lang_data_raw$p_wt_30

# scale syntax input

lang_data_raw$Z1_syntax_rescale <- lang_data_raw$p_cps_14 * 100
lang_data_raw$Z2_syntax_rescale <- lang_data_raw$p_cps_30 * 100

# get predicted PPVT outcomes on basis of cumulative model

lang_data_raw$pred_ppvt <- v_m3$estimate[1] + v_m3$estimate[2] * lang_data_raw$cumulative_vocab


# add variables denoting if pred_ppvt in top or bottom 16

lang_data_raw <- mutate(lang_data_raw, voc_hi_16 = (min_rank(-pred_ppvt) <= 16),
                        voc_lo_16 = min_rank(-pred_ppvt) >= 46)

# colour code by these variables
lang_data_raw <- mutate(lang_data_raw, voc_color = ifelse(voc_hi_16, "a",
                                                      ifelse(voc_lo_16, "c", "b")))

# make vocab plot

vocab <- ggplot(lang_data_raw, aes(x = p_wt_14, y = p_wt_30, group = voc_color))

vocab + theme_bw() + geom_point(aes(color = voc_color)) +
  scale_color_manual(values = c("#67a9cf", "grey", "#ef8a62")) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  theme(legend.position = "none") +
  scale_x_continuous("\nParent word types at 14 months", limits = c(60, 750)) +
  scale_y_continuous("Parent word types at 30 months\n", limits = c(60, 750)) +
  coord_fixed(ratio = 1)
ggsave("VocPred.png", width=5, height=5)


# same for syntax

# get predicted CELF outcomes on basis of differing effects model

lang_data_raw$pred_celf <- s_m2$estimate[1] + s_m2$estimate[2] * lang_data_raw$Z1_syntax_rescale +
  s_m2$estimate[3] * lang_data_raw$Z2_syntax_rescale


# add variables denoting if pred_celf in top or bottom 16

lang_data_raw <- mutate(lang_data_raw, syn_hi_16 = (min_rank(-pred_celf) <= 16),
                        syn_lo_16 = min_rank(-pred_celf) >= 46)

# colour code by these variables
lang_data_raw <- mutate(lang_data_raw, syn_color = ifelse(syn_hi_16, "a",
                                                      ifelse(syn_lo_16, "c", "b")))

syntax <- ggplot(lang_data_raw, aes(x = Z1_syntax_rescale, y = Z2_syntax_rescale, group = syn_color))

syntax + theme_bw() + geom_point(aes(color = syn_color)) +
  scale_color_manual(values = c("#67a9cf", "grey", "#ef8a62")) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  theme(legend.position = "none") +
  scale_x_continuous("\nParent clauses per 100 sentences at 14 months", limits = c(100, 130)) +
  scale_y_continuous("Parent clauses per 100 sentences at 30 months\n", limits = c(100, 130)) +
  coord_fixed(ratio = 1)
ggsave("SynPred.png", width=5, height=5)


# 5) SENSITIVITY ANALYSIS (Text S8)

# goal: to estimate the biasing effect of an unobserved confounder
# as strong as any of our observed confounders (except parent IQ and household income)
# and to assess robustness to differing functional forms


# A) Vocabulary

# use one complete imputed data set

complete_data <- complete(lang_data_imp)

# for this analysis, treat child birth order and gender as continuous

complete_data$cbo_cont <- as.numeric(complete_data$c_bo)
complete_data$cgen_cont <- as.numeric(complete_data$c_gen)

# Regress cumulative input on parent IQ and income and save the residuals from this model

v_sens_model <- lm(cumulative_vocab ~ p_viq + p_inc, complete_data)
v_sens_residuals <- residuals(v_sens_model)

# Compute standard deviation of these residuals

v_sens_res_sd <- sd(v_sens_residuals)

# Regress Y on parent IQ and income and save the residuals from this model

v_sens_model_2 <- lm(c_ppvt_74 ~ p_viq + p_inc, complete_data)
v_sens_residuals_2 <- residuals(v_sens_model_2)

# compute the standard deviation of these residuals

v_sens_res_2_sd <- sd(v_sens_residuals_2)

# Regress each possible "unobserved" covariate u on parent IQ and income
# save the residuals, compute the SD

# child birth order

cbo_sens <- lm(cbo_cont ~ p_viq + p_inc, complete_data)
cbo_sens_res <- residuals(cbo_sens)
cbo_sens_res_sd <- sd(cbo_sens_res)

# child gender

cgen_sens <- lm(cgen_cont ~ p_viq + p_inc, complete_data)
cgen_sens_res <- residuals(cgen_sens)
cgen_sens_res_sd <- sd(cgen_sens_res)

# child word types at 14 months

cwt14_sens <- lm(c_wt_14 ~ p_viq + p_inc, complete_data)
cwt14_sens_res <- residuals(cwt14_sens)
cwt14_sens_res_sd <- sd(cwt14_sens_res)

# child gesture types at 14 months

cgt14_sens <- lm(c_gt_14 ~ p_viq + p_inc, complete_data)
cgt14_sens_res <- residuals(cgt14_sens)
cgt14_sens_res_sd <- sd(cgt14_sens_res)

# peduc

peduc_sens <- lm(p_educ ~ p_viq + p_inc, complete_data)
peduc_sens_res <- residuals(peduc_sens)
peduc_sens_res_sd <- sd(peduc_sens_res)

# Calculate the correlations among residuals

# correlation of ec and ey

v_ec_ey <- cor.test(v_sens_residuals, v_sens_residuals_2)$estimate

# correlation of ec and eu for each u

v_ec_cbo <- cor.test(v_sens_residuals, cbo_sens_res)$estimate
v_ec_cgen <- cor.test(v_sens_residuals, cgen_sens_res)$estimate
v_ec_cwt14 <- cor.test(v_sens_residuals, cwt14_sens_res)$estimate
v_ec_cgt14 <- cor.test(v_sens_residuals, cgt14_sens_res)$estimate
v_ec_peduc <- cor.test(v_sens_residuals, peduc_sens_res)$estimate

# correlation of eu and ey for each u

v_ey_cbo <- cor.test(v_sens_residuals_2, cbo_sens_res)$estimate
v_ey_cgen <- cor.test(v_sens_residuals_2, cgen_sens_res)$estimate
v_ey_cwt14 <- cor.test(v_sens_residuals_2, cwt14_sens_res)$estimate
v_ey_cgt14 <- cor.test(v_sens_residuals_2, cgt14_sens_res)$estimate
v_ey_peduc <- cor.test(v_sens_residuals_2, peduc_sens_res)$estimate

# for each of the 5 'omitted covariates' U, compute the bias

# child birth order

v_cbo_gamma <- ((v_ey_cbo - v_ec_cbo * v_ec_ey) * v_sens_res_2_sd) / ((1 - v_ec_cbo ^ 2) * cbo_sens_res_sd)

v_cbo_alpha <- (v_ec_cbo * cbo_sens_res_sd) / v_sens_res_sd

v_cbo_bias <- v_cbo_gamma * v_cbo_alpha

print(v_cbo_bias)

# child gender

v_cgen_gamma <- ((v_ey_cgen - v_ec_cgen * v_ec_ey) * v_sens_res_2_sd) / ((1 - v_ec_cgen ^ 2) * cgen_sens_res_sd)

v_cgen_alpha <- (v_ec_cgen * cgen_sens_res_sd) / v_sens_res_sd

v_cgen_bias <- v_cgen_gamma * v_cgen_alpha

print(v_cgen_bias)

# child word types at 14 months

v_cwt14_gamma <- ((v_ey_cwt14 - v_ec_cwt14 * v_ec_ey) * v_sens_res_2_sd) / ((1 - v_ec_cwt14 ^ 2) * cwt14_sens_res_sd)

v_cwt14_alpha <- (v_ec_cwt14 * cwt14_sens_res_sd) / v_sens_res_sd

v_cwt14_bias <- v_cwt14_gamma * v_cwt14_alpha

print(v_cwt14_bias)

# child gesture types at 14 months

v_cgt14_gamma <- ((v_ey_cgt14 - v_ec_cgt14 * v_ec_ey) * v_sens_res_2_sd) / ((1 - v_ec_cgt14 ^ 2) * cgt14_sens_res_sd)

v_cgt14_alpha <- (v_ec_cgt14 * cgt14_sens_res_sd) / v_sens_res_sd

v_cgt14_bias <- v_cgt14_gamma * v_cgt14_alpha

print(v_cgt14_bias)

# parent education

v_peduc_gamma <- ((v_ey_peduc - v_ec_peduc * v_ec_ey) * v_sens_res_2_sd) / ((1 - v_ec_peduc ^ 2) * peduc_sens_res_sd)

v_peduc_alpha <- (v_ec_peduc * peduc_sens_res_sd) / v_sens_res_sd

v_peduc_bias <- v_peduc_gamma * v_peduc_alpha

print(v_peduc_bias)

# retrieve estimate from constant effects model v_m3 above
vocab_estimate <- v_m3[2, 'estimate']

v_cbo_estimate <- vocab_estimate - v_cbo_bias
v_cgen_estimate <- vocab_estimate - v_cgen_bias
v_cwt14_estimate <- vocab_estimate - v_cwt14_bias
v_cgt14_estimate <- vocab_estimate - v_cgt14_bias
v_peduc_estimate <- vocab_estimate - v_peduc_bias

# here, the strongest bias is for cgt14
# ratio of estimate ignoring the confounder and estimate taking the confounder into account
vocab_estimate / v_cgt14_estimate



# B) Syntax

# Regress Z1 syntax input on parent IQ and income and save the residuals from this model

s_sens_model_Z1 <- lm(Z1_syntax_rescale ~ p_viq + p_inc, complete_data)
s_sens_residuals_Z1 <- residuals(s_sens_model_Z1)

# Compute standard deviation of these residuals

s_sens_res_Z1_sd <- sd(s_sens_residuals_Z1)

# Regress Z2 syntax input on parent IQ and income and save the residuals from this model

s_sens_model_Z2 <- lm(Z2_syntax_rescale ~ p_viq + p_inc, complete_data)
s_sens_residuals_Z2 <- residuals(s_sens_model_Z2)

# Compute standard deviation of these residuals

s_sens_res_Z2_sd <- sd(s_sens_residuals_Z2)

# Regress Y on parent IQ and income and save the residuals from this model

s_sens_model_Y <- lm(c_celf_rs ~ p_viq + p_inc, complete_data)
s_sens_residuals_Y <- residuals(s_sens_model_Y)

# compute the standard deviation of these residuals

s_sens_res_Y_sd <- sd(s_sens_residuals_Y)

# We already have residuals and SDs from each 'unobserved' covariate, calculated above

# Calculate the correlations among residuals

# correlation of ez1 and ey

s_ez1_ey <- cor.test(s_sens_residuals_Z1, s_sens_residuals_Y)$estimate

# correlation of ez2 and ey

s_ez2_ey <- cor.test(s_sens_residuals_Z2, s_sens_residuals_Y)$estimate

# i) Z1

# correlation of ez1 and eu for each u

s_ez1_cbo <- cor.test(s_sens_residuals_Z1, cbo_sens_res)$estimate
s_ez1_cgen <- cor.test(s_sens_residuals_Z1, cgen_sens_res)$estimate
s_ez1_cwt14 <- cor.test(s_sens_residuals_Z1, cwt14_sens_res)$estimate
s_ez1_cgt14 <- cor.test(s_sens_residuals_Z1, cgt14_sens_res)$estimate
s_ez1_peduc <- cor.test(s_sens_residuals_Z1, peduc_sens_res)$estimate

# correlation of eu and ey for each u

s_ey_cbo <- cor.test(s_sens_residuals_Y, cbo_sens_res)$estimate
s_ey_cgen <- cor.test(s_sens_residuals_Y, cgen_sens_res)$estimate
s_ey_cwt14 <- cor.test(s_sens_residuals_Y, cwt14_sens_res)$estimate
s_ey_cgt14 <- cor.test(s_sens_residuals_Y, cgt14_sens_res)$estimate
s_ey_peduc <- cor.test(s_sens_residuals_Y, peduc_sens_res)$estimate

# correlation between Z1 and Z2 residuals

s_Z1_Z2 <- cor.test(s_sens_residuals_Z1, s_sens_residuals_Z2)$estimate

# for each of the 5 'omitted covariates' U, compute the bias

# child birth order

s_z1_cbo_gamma <- (((1 - (s_Z1_Z2 ^ 2)) * s_ey_cbo - s_ez1_cbo * s_ez1_ey - s_ez2_ey * s_ez1_cbo * s_Z1_Z2) * (s_sens_res_Y_sd/cbo_sens_res_sd)) /
  (1 - (s_Z1_Z2 ^ 2) - (s_ez1_cbo ^ 2))

s_z1_cbo_alpha <- (s_ez1_cbo * cbo_sens_res_sd) / s_sens_res_Z1_sd

s_z1_cbo_bias <- s_z1_cbo_gamma * s_z1_cbo_alpha

print(s_z1_cbo_bias)

# child gender

s_z1_cgen_gamma <- (((1 - (s_Z1_Z2 ^ 2)) * s_ey_cgen - s_ez1_cgen * s_ez1_ey - s_ez2_ey * s_ez1_cgen * s_Z1_Z2) * (s_sens_res_Y_sd/cgen_sens_res_sd)) /
  (1 - (s_Z1_Z2 ^ 2) - (s_ez1_cgen ^ 2))

s_z1_cgen_alpha <- (s_ez1_cgen * cgen_sens_res_sd) / s_sens_res_Z1_sd

s_z1_cgen_bias <- s_z1_cgen_gamma * s_z1_cgen_alpha

print(s_z1_cgen_bias)

# child word types at 14 months

s_z1_cwt14_gamma <- (((1 - (s_Z1_Z2 ^ 2)) * s_ey_cwt14 - s_ez1_cwt14 * s_ez1_ey - s_ez2_ey * s_ez1_cwt14 * s_Z1_Z2) * (s_sens_res_Y_sd/cwt14_sens_res_sd)) /
  (1 - (s_Z1_Z2 ^ 2) - (s_ez1_cwt14 ^ 2))

s_z1_cwt14_alpha <- (s_ez1_cwt14 * cwt14_sens_res_sd) / s_sens_res_Z1_sd

s_z1_cwt14_bias <- s_z1_cwt14_gamma * s_z1_cwt14_alpha

print(s_z1_cwt14_bias)

# child gesture types at 14 months

s_z1_cgt14_gamma <- (((1 - (s_Z1_Z2 ^ 2)) * s_ey_cgt14 - s_ez1_cgt14 * s_ez1_ey - s_ez2_ey * s_ez1_cgt14 * s_Z1_Z2) * (s_sens_res_Y_sd/cgt14_sens_res_sd)) /
  (1 - (s_Z1_Z2 ^ 2) - (s_ez1_cgt14 ^ 2))

s_z1_cgt14_alpha <- (s_ez1_cgt14 * cgt14_sens_res_sd) / s_sens_res_Z1_sd

s_z1_cgt14_bias <- s_z1_cgt14_gamma * s_z1_cgt14_alpha

print(s_z1_cgt14_bias)

# parent education

s_z1_peduc_gamma <- (((1 - (s_Z1_Z2 ^ 2)) * s_ey_peduc - s_ez1_peduc * s_ez1_ey - s_ez2_ey * s_ez1_peduc * s_Z1_Z2) * (s_sens_res_Y_sd/peduc_sens_res_sd)) /
  (1 - (s_Z1_Z2 ^ 2) - (s_ez1_peduc ^ 2))

s_z1_peduc_alpha <- (s_ez1_peduc * peduc_sens_res_sd) / s_sens_res_Z1_sd

s_z1_peduc_bias <- s_z1_peduc_gamma * s_z1_peduc_alpha

print(s_z1_peduc_bias)

# get syntax Z1 effect estimate from model s_m2 above
syntax_Z1_estimate <- s_m2[2, 'estimate']

s_z1_cbo_estimate <- syntax_Z1_estimate - s_z1_cbo_bias
s_z1_cgen_estimate <- syntax_Z1_estimate - s_z1_cgen_bias
s_z1_cwt14_estimate <- syntax_Z1_estimate - s_z1_cwt14_bias
s_z1_cgt14_estimate <- syntax_Z1_estimate - s_z1_cgt14_bias
s_z1_peduc_estimate <- syntax_Z1_estimate - s_z1_peduc_bias

# here, the strongest bias is for cwt14
# ratio of estimate ignoring the confounder and estimate taking the confounder into account
syntax_Z1_estimate / s_z1_cwt14_estimate



# ii) Z2

# correlation of ez2 and eu for each u

s_ez2_cbo <- cor.test(s_sens_residuals_Z2, cbo_sens_res)$estimate
s_ez2_cgen <- cor.test(s_sens_residuals_Z2, cgen_sens_res)$estimate
s_ez2_cwt14 <- cor.test(s_sens_residuals_Z2, cwt14_sens_res)$estimate
s_ez2_cgt14 <- cor.test(s_sens_residuals_Z2, cgt14_sens_res)$estimate
s_ez2_peduc <- cor.test(s_sens_residuals_Z2, peduc_sens_res)$estimate

# for each of the 5 'omitted covariates' U, compute the bias

# child birth order

s_z2_cbo_gamma <- (((1 - (s_Z1_Z2 ^ 2)) * s_ey_cbo - s_ez2_cbo * s_ez2_ey - s_ez1_ey * s_ez2_cbo * s_Z1_Z2) * (s_sens_res_Y_sd/cbo_sens_res_sd)) /
  (1 - (s_Z1_Z2 ^ 2) - (s_ez2_cbo ^ 2))

s_z2_cbo_alpha <- (s_ez2_cbo * cbo_sens_res_sd) / s_sens_res_Z2_sd

s_z2_cbo_bias <- s_z2_cbo_gamma * s_z2_cbo_alpha

print(s_z2_cbo_bias)

# child gender

s_z2_cgen_gamma <- (((1 - (s_Z1_Z2 ^ 2)) * s_ey_cgen - s_ez2_cgen * s_ez2_ey - s_ez1_ey * s_ez2_cgen * s_Z1_Z2) * (s_sens_res_Y_sd/cgen_sens_res_sd)) /
  (1 - (s_Z1_Z2 ^ 2) - (s_ez2_cgen ^ 2))

s_z2_cgen_alpha <- (s_ez2_cgen * cgen_sens_res_sd) / s_sens_res_Z2_sd

s_z2_cgen_bias <- s_z2_cgen_gamma * s_z2_cgen_alpha

print(s_z2_cgen_bias)

# child word types at 14 months

s_z2_cwt14_gamma <- (((1 - (s_Z1_Z2 ^ 2)) * s_ey_cwt14 - s_ez2_cwt14 * s_ez2_ey - s_ez1_ey * s_ez2_cwt14 * s_Z1_Z2) * (s_sens_res_Y_sd/cwt14_sens_res_sd)) /
  (1 - (s_Z1_Z2 ^ 2) - (s_ez2_cwt14 ^ 2))

s_z2_cwt14_alpha <- (s_ez2_cwt14 * cwt14_sens_res_sd) / s_sens_res_Z2_sd

s_z2_cwt14_bias <- s_z2_cwt14_gamma * s_z2_cwt14_alpha

print(s_z2_cwt14_bias)

# child gesture types at 14 months

s_z2_cgt14_gamma <- (((1 - (s_Z1_Z2 ^ 2)) * s_ey_cgt14 - s_ez2_cgt14 * s_ez2_ey - s_ez1_ey * s_ez2_cgt14 * s_Z1_Z2) * (s_sens_res_Y_sd/cgt14_sens_res_sd)) /
  (1 - (s_Z1_Z2 ^ 2) - (s_ez2_cgt14 ^ 2))

s_z2_cgt14_alpha <- (s_ez2_cgt14 * cgt14_sens_res_sd) / s_sens_res_Z2_sd

s_z2_cgt14_bias <- s_z2_cgt14_gamma * s_z2_cgt14_alpha

print(s_z2_cgt14_bias)

# parent education

s_z2_peduc_gamma <- (((1 - (s_Z1_Z2 ^ 2)) * s_ey_peduc - s_ez2_peduc * s_ez2_ey - s_ez1_ey * s_ez2_peduc * s_Z1_Z2) * (s_sens_res_Y_sd/peduc_sens_res_sd)) /
  (1 - (s_Z1_Z2 ^ 2) - (s_ez2_peduc ^ 2))

s_z2_peduc_alpha <- (s_ez2_peduc * peduc_sens_res_sd) / s_sens_res_Z2_sd

s_z2_peduc_bias <- s_z2_peduc_gamma * s_z2_peduc_alpha

print(s_z2_peduc_bias)


# get syntax Z2 effect estimate from model s_m2 above

syntax_Z2_estimate <- s_m2[3, 'estimate']

s_z2_cbo_estimate <- syntax_Z2_estimate - s_z2_cbo_bias
s_z2_cgen_estimate <- syntax_Z2_estimate - s_z2_cgen_bias
s_z2_cwt14_estimate <- syntax_Z2_estimate - s_z2_cwt14_bias
s_z2_cgt14_estimate <- syntax_Z2_estimate - s_z2_cgt14_bias
s_z2_peduc_estimate <- syntax_Z2_estimate - s_z2_peduc_bias

# here, the strongest bias is for peduc
# ratio of estimate ignoring the confounder and estimate taking the confounder into account
syntax_Z2_estimate / s_z2_peduc_estimate


# We also want to test our models' robustness to violations of the linearity assumption
# To do this, we compare model graphs for linear, quadratic, and log

# again, do this with one complete dataset

complete_data <- complete(lang_data_imp)

# A) Vocabulary

# retrieve linear model of cumulative input using this dataset

v_lin <- summary(v_const$analyses[[1]])

v_linear_intercept <- v_lin$coefficients["(Intercept)", "Estimate"]
v_linear_slope <- v_lin$coefficients["cumulative_vocab", "Estimate"]

# define function for plotting linear model graph

v_linearFun <- function(x) {
  v_linear_intercept + v_linear_slope * x
}

# quadratic model

# center and square cumulative vocabulary input

complete_data["cumul_center"] <- complete_data$cumulative_vocab - (mean(complete_data$cumulative_vocab))
complete_data["cumul_center_sq"] <- complete_data$cumul_center ^ 2

# cumulative model quadratic
cumul_model_quad <- lm(c_ppvt_74 ~ cumul_center + cumul_center_sq, weights = combqbw_vocab, complete_data)
v_quad <- summary(cumul_model_quad)

v_quad_intercept <- v_quad$coefficients["(Intercept)", "Estimate"]
v_quad_slope <- v_quad$coefficients["cumul_center", "Estimate"]
v_quad_accel <- v_quad$coefficients["cumul_center_sq", "Estimate"]

# define function for plotting quadratic model graph

v_quadFun <- function(x) {
  v_quad_intercept + v_quad_slope * (x - mean(complete_data$cumulative_vocab)) + v_quad_accel * ((x - mean(complete_data$cumulative_vocab)) ^ 2)
}

# log model

# log-transform cumulative vocab input

complete_data['log_cumulative'] <- log(complete_data$cumulative_vocab)

logcumul_model <- lm(c_ppvt_74 ~ log_cumulative, weights = combqbw_vocab, complete_data)
v_log <- summary(logcumul_model)

v_log_intercept <- v_log$coefficients["(Intercept)", "Estimate"]
v_log_slope <- v_log$coefficients["log_cumulative", "Estimate"]

# define function for plotting log model

v_logFun <- function(x) {
  v_log_intercept + v_log_slope * log(x)
}

# plot model graphs between 10th and 90th percentiles

cumul_10 <- quantile(complete_data$cumulative_vocab, .1)
cumul_90 <- quantile(complete_data$cumulative_vocab, .9)

v_models <- ggplot(data.frame(x = c(cumul_10, cumul_90)), aes(x = x)) +
  stat_function(fun = v_linearFun, aes(linetype = "Linear")) + stat_function(fun = v_quadFun, aes(linetype = "Quadratic")) + stat_function(fun = v_logFun, aes(linetype = "Logarithmic"))

v_models + theme_bw() + scale_x_continuous("\nCumulative vocabulary input") + scale_y_continuous("Predicted PPVT score\n", limits = c(0,140)) + scale_linetype_discrete("Model")
ggsave("FunctionalFormVocab_final.png")


# B. Syntax
# Here we plot Z1 and Z2 separately, holding the other constant at the mean

# Z1

# retrieve linear model of cumulative input using this dataset

s_lin <- summary(s_sep$analyses[[1]])

s_linear_intercept <- s_lin$coefficients["(Intercept)", "Estimate"]
s_linear_slope_Z1 <- s_lin$coefficients["Z1_syntax_rescale", "Estimate"]
s_linear_slope_Z2 <- s_lin$coefficients["Z2_syntax_rescale", "Estimate"]

# define function for plotting linear model graph

s_z1_linearFun <- function(x) {
  s_linear_intercept + s_linear_slope_Z1 * x + s_linear_slope_Z2 * mean(complete_data$Z2_syntax_rescale)
}

# quadratic model

# first center and square Z1 and Z2

complete_data["Z1_center"] <- complete_data$Z1_syntax_rescale - (mean(complete_data$Z1_syntax_rescale))
complete_data["Z1_center_sq"] <- complete_data$Z1_center ^ 2

complete_data["Z2_center"] <- complete_data$Z2_syntax_rescale - (mean(complete_data$Z2_syntax_rescale))
complete_data["Z2_center_sq"] <- complete_data$Z2_center ^ 2

sep_model_quad <- lm(c_celf_rs ~ Z1_center + Z1_center_sq + Z2_center + Z2_center_sq, weights = combqbw_syntax, complete_data)
s_quad <- summary(sep_model_quad)

s_quad_intercept <- s_quad$coefficients["(Intercept)", "Estimate"]

s_quad_slope_Z1 <- s_quad$coefficients["Z1_center", "Estimate"]
s_quad_accel_Z1 <- s_quad$coefficients["Z1_center_sq", "Estimate"]

s_quad_slope_Z2 <- s_quad$coefficients["Z2_center", "Estimate"]
s_quad_accel_Z2 <- s_quad$coefficients["Z2_center_sq", "Estimate"]


# define function for plotting quadratic model graph

s_z1_quadFun <- function(x) {
  s_quad_intercept + s_quad_slope_Z1 * (x - mean(complete_data$Z1_syntax_rescale)) + s_quad_accel_Z1 * ((x - mean(complete_data$Z1_syntax_rescale)) ^ 2) +
    s_quad_slope_Z2 * mean(complete_data$Z2_center) + s_quad_accel_Z2 * (mean(complete_data$Z2_center) ^ 2)
}

# log model

# log-transform Z1 and Z2

complete_data['log_Z1'] <- log(complete_data$Z1_syntax_rescale)
complete_data['log_Z2'] <- log(complete_data$Z2_syntax_rescale)


logsep_model <- lm(c_celf_rs ~ log_Z1 + log_Z2, weights = combqbw_syntax, complete_data)
s_log <- summary(logsep_model)

s_log_intercept <- s_log$coefficients["(Intercept)", "Estimate"]
s_log_slope_Z1 <- s_log$coefficients["log_Z1", "Estimate"]
s_log_slope_Z2 <- s_log$coefficients["log_Z2", "Estimate"]

# define function for plotting log model

s_z1_logFun <- function(x) {
  s_log_intercept + s_log_slope_Z1 * log(x) + s_log_slope_Z2 * mean(complete_data$log_Z2)
}

# plot model graphs between 10th and 90th percentiles

Z1_10 <- quantile(complete_data$Z1_syntax_rescale, .1)
Z1_90 <- quantile(complete_data$Z1_syntax_rescale, .9)

s_z1_models <- ggplot(data.frame(x = c(Z1_10, Z1_90)), aes(x = x)) +
  stat_function(fun = s_z1_linearFun, aes(linetype = "Linear")) + stat_function(fun = s_z1_quadFun, aes(linetype = "Quadratic")) + stat_function(fun = s_z1_logFun, aes(linetype = "Logarithmic"))

s_z1_models + theme_bw() + scale_x_continuous("\nZ1 syntax input (Z2 held constant at mean)") + scale_y_continuous("Predicted CELF score\n", limits=c(0,16)) + scale_linetype_discrete("Model")
ggsave("FunctionalFormSyntaxZ1_final.png")


# Z2

# linear model is as before

# define function for plotting

s_z2_linearFun <- function(x) {
  s_linear_intercept + s_linear_slope_Z1 * mean(complete_data$Z1_syntax_rescale) + s_linear_slope_Z2 * x
}

# quadratic model is as before

# define function for plotting

s_z2_quadFun <- function(x) {
  s_quad_intercept + s_quad_slope_Z1 * mean(complete_data$Z1_center) + s_quad_accel_Z1 * (mean(complete_data$Z1_center) ^ 2) +
    s_quad_slope_Z2 * (x - mean(complete_data$Z2_syntax_rescale)) + s_quad_accel_Z2 * ((x - mean(complete_data$Z2_syntax_rescale)) ^ 2)
}

# log model is as before

# define function for plotting

s_z2_logFun <- function(x) {
  s_log_intercept + s_log_slope_Z1 * mean(complete_data$log_Z1) + s_log_slope_Z2 * log(x)
}

# plot model graphs between 10th and 90th percentiles

Z2_10 <- quantile(complete_data$Z2_syntax_rescale, .1)
Z2_90 <- quantile(complete_data$Z2_syntax_rescale, .9)


s_z2_models <- ggplot(data.frame(x = c(Z2_10, Z2_90)), aes(x = x)) +
  stat_function(fun = s_z2_linearFun, aes(linetype = "Linear")) + stat_function(fun = s_z2_quadFun, aes(linetype = "Quadratic")) + stat_function(fun = s_z2_logFun, aes(linetype = "Logarithmic"))

s_z2_models + theme_bw() + scale_x_continuous("\nZ2 syntax input (Z1 held constant at mean)") + scale_y_continuous("Predicted CELF score\n", limits=c(0,16)) + scale_linetype_discrete("Model")
ggsave("FunctionalFormSyntaxZ2_final.png")



