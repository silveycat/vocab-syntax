# Simulation to determine optimal number of strata for quantile binning analysis
# Reported in Text S4 of Supplemental Materials of:
# 'Effects of time-varying parent input on child language outcomes differ for vocabulary and syntax'

# get required packages
required_packages <- c("dplyr", "lmtest", "sandwich", "MASS", "ggplot2")

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

# define function for producing robust regression output (reproduces robust results from Stata regressions)
# takes a model x as argument
robust <- function(x) {
  coeftest(x, vcov = vcovHC(x, "HC1"))
}


# runs parameter - number of random datasets to generate & analyse
runs <- 1000

# set the number of observations to 64 (same as our sample)
N <- 64

# confounding parameter = association between X and Z
confounding <- 0.5

# adjust residual variances accordingly to keep variance of Z and Y to 1
Zresvar <- 0.75
Yresvar <- 0.51

# number of strata to test: from 3 to 64
stratum_values <- seq(3, 64)


# main simulation function

sim_run <- function(sim, runs, N, numstrata, confounding) {
  # data frame for collecting results
  results <- data.frame((matrix(0, ncol=6, nrow=0)))
  names(sim_results) <- c("sim", "strata", "run", "estimate", "bias", "weight_var")
  # variable for collecting the bias
  overall_bias <- 0
  # main simulation loop. runs the number of times specified in variable 'runs'
  # generates data, analyzes, & writes data from each run to file
  for (i in 1:runs) {
    # print what run we're on
    # print(i)
    # step 1: generate random variable X with N observations, mean 0, and variance 1
    X <- rnorm(N, 0, 1)
    # convert to data frame
    frame <- data.frame(X)
    # step 2: generate Z = cX + e, where c = confounding and e ~ N (0, Zresvar)
    # rnorm function takes mean and standard deviation as arguments, so we square root our variance
    frame['Z'] <- confounding * frame$X + rnorm(N, 0, sqrt(Zresvar))
    # step 3: generate Y = 0.4Z + 0.5X + f, where f ~ N (0, Yresvar)
    frame['Y'] <- 0.4 * frame$Z + 0.5 * frame$X + rnorm(N, 0, sqrt(Yresvar))
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
    # get Z coefficient
    Zcoef <- outcome_robust["Z", 1]
    # get bias: estimated effect - true effect
    bias <- Zcoef - 0.4
    # add to ongoing bias total
    overall_bias <- overall_bias + bias
    # add to data frame
    results <- rbind(results, data.frame(sim = sim, strata = as.numeric(numstrata), run = i, estimate = Zcoef,
                                         bias = bias, weight_var = weight_var))
  }
  return(results)
}


sim_results <- data.frame((matrix(0, ncol=6, nrow=0)))
names(sim_results) <- c("sim", "strata", "run", "estimate", "bias", "weight_var")

# main simulation loop

for (i in 1:5) {
  current_sim <- i
  print(paste0("simulation ", current_sim))
  for (s in stratum_values) {
    print (paste0(s, " strata"))
    result <- sim_run(current_sim, runs, N, s, confounding)
    sim_results <- rbind(sim_results, result)
  }
}


all_bias_means <- sim_results %>% group_by(sim, strata) %>% summarise(mean_bias = mean(bias))

together_bias <- ggplot(all_bias_means, aes(x=strata, y=mean_bias, group=sim))
together_bias + theme_bw() + geom_line(aes(group=sim), alpha = 0.5) +
  scale_x_continuous("\nNumber of quantiles",
                       breaks=c(8,16,24,32,40,48,56,64)) +
  scale_y_continuous("Mean bias\n")

ggsave("BiasAll.png")

all_weightvar_means <- sim_results %>% group_by(sim, strata) %>% summarise(mean_weightvar = mean(weight_var))

together_weightvar <- ggplot(all_weightvar_means, aes(x=strata, y=mean_weightvar, group=sim))
together_weightvar + theme_bw() + geom_line(aes(group=sim), alpha = 0.5) +
  scale_x_continuous("\nNumber of quantiles",
                     breaks=c(8,16,24,32,40,48,56,64)) +
  scale_y_continuous("Mean variance in weights\n")

ggsave("WeightvarAll.png")
