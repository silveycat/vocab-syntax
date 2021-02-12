# Simulation to determine whether quantile binning method returns robust results in small samples
# Reported in Text S5 of Supplemental Materials of:
# 'Effects of time-varying parent input on child language outcomes differ for vocabulary and syntax'

# get required packages
required_packages <- c("dplyr", "lmtest", "sandwich", "MASS", "ggplot2", "plotrix")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# use dplyr for summarising data and making quintiles
library(dplyr)
# use lmtest and sandwich libraries for robust regressions
library(lmtest)
library(sandwich)
# use MASS library for running ordinal regressions
library(MASS)
# ggplot2 for plotting
library(ggplot2)
# plotrix for SE
library(plotrix)

# define function for producing robust regression output (reproduces robust results from Stata regressions)
# takes a model x as argument
robust <- function(x) {
  coeftest(x, vcov = vcovHC(x, "HC1"))
}

# runs parameter - number of random datasets to generate & analyse
runs <- 1000

# set the number of observations - small (64) or large (500)
N <- c(64,500)

# set number of strata - previous simulations found best numbers are 8 for N=64 and 100 for N=500
numstrata <- c(8, 100)

# confounding parameter = association between X and Z. we assume same for effect of X on Y
confounding <- 0.5

# adjust residual variances accordingly to keep variance of Z and Y to 1
Zresvar <- 0.75
Yresvar <- 0.51


# main simulation function

# calculates:
# 1) bias = mean over all runs of error (= difference between estimated and true effect)
# 2) does average SE reported in models on each run correspond to empirical SE
# = SD of errors over all runs?

sim_run <- function(runs, N, numstrata, confounding) {
  # data frame for collecting results
  results <- data.frame((matrix(0, ncol=6, nrow=0)))
  names(sim_results) <- c("strata", "run", "estimate", "error", "weight_var")
  # variables for collecting the bias and SE
  overall_bias <- 0
  overall_se <- 0
  # main simulation loop. runs the number of times specified in variable 'runs'
  # generates data, analyzes, & writes data from each run to file
  for (i in 1:runs) {
    # print what run we're on
    print(i)
    # step 1: generate random variable X with N observations, mean 0, and variance 1
    X <- rnorm(N, 0, 1)
    # convert to data frame
    frame <- data.frame(X)
    # step 2: generate Z = cX + e, where c = confounding and e ~ N (0, Zresvar)
    # rnorm function takes mean and standard deviation as arguments, so we square root our variance
    frame['Z'] <- confounding * frame$X + rnorm(N, 0, sqrt(Zresvar))
    # step 3: generate Y = 0.4Z + 0.5X + f, where f ~ N (0, Yresvar)
    frame['Y'] <- 0.3 * frame$Z + 0.5 * frame$X + rnorm(N, 0, sqrt(Yresvar))
    # step 4: use quantile binning to do a weighted analysis
    # first stratify the sample on Z
    # into the given number of strata
    frame['Zstratum'] <- factor(ntile(frame$Z, numstrata))
    # predict stratum assignment from X
    ordinal_model <- polr(Zstratum ~ X, frame)
    # save the probabilities of being in each stratum as a new data frame, predprobs
    predprobs <- data.frame(predict(ordinal_model, type="p"))
    # save the probability for each person of being in their own stratum
    # first, make the column names match by assigning each stratum number to the appropriate column heading
    colnames(predprobs) <- seq_len(numstrata)
    # add actual stratum to this data frame
    predprobs['Zstratum'] <- frame$Zstratum
    # personprob = probability for each person of being in their actual stratum as listed in the Zstratum column
    frame['personprob'] <- as.numeric(predprobs[cbind(seq_len(nrow(predprobs)), match(predprobs$Zstratum, names(predprobs)))])
    # calculate the weights
    # numerator = marginal probability of being in each stratum, i.e. 1/8
    # denominator = person-specific probability of being in their stratum
    # qbw = quantile binning weight
    frame['qbw'] <- 1/numstrata/(frame$personprob)
    # record variance in the weights
    weight_var <- var(frame$qbw)
    # run the weighted outcome model
    outcome <- lm(Y ~ Z, weights = qbw, frame)
    outcome_robust <- robust(outcome)
    # get Z coefficient and standard error
    Zcoef <- outcome_robust["Z", 1]
    Zse <- outcome_robust["Z", 2]
    # get bias: estimated effect - true effect
    error <- Zcoef - 0.4
    # add to ongoing bias total
    overall_bias <- overall_bias + error
    # add SE to ongoing total for averaging
    overall_se <- overall_se + Zse
    # add to data frame
    results <- rbind(results, data.frame(N = N, strata = as.numeric(numstrata), run = i, estimate = Zcoef,
                                         error = error, se = Zse, weight_var = weight_var))
  }
  return(results)
}


sim_results <- data.frame((matrix(0, ncol=8, nrow=0)))
names(sim_results) <- c("N", "strata", "run", "estimate", "error", "se", "weight_var")


# run for small sample
small_result <- sim_run(runs, N[1], numstrata[1], confounding)
sim_results <- rbind(sim_results, small_result)
# run for large sample
large_result <- sim_run(runs, N[2], numstrata[2], confounding)
sim_results <- rbind(sim_results, large_result)


small_samp <- subset(sim_results, N == 64)
large_samp <- subset(sim_results, N == 500)

small_bias <- mean(small_samp$error)
large_bias <- mean(large_samp$error)

print(small_bias)
print(large_bias)

# assess how empirical SE compares to SE estimated in each model

# empirical SE = SD of errors across all runs
# estimated SE = average model SE across all runs

# for small sample
emp_SE_small <- sd(small_samp$error)
est_SE_small <- mean(small_samp$se)

print(emp_SE_small)
print(est_SE_small)

# for large sample
emp_SE_large <- sd(large_samp$error)
est_SE_large <- mean(large_samp$se)

print(emp_SE_large)
print(est_SE_large)