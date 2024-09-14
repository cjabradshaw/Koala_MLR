############################################
## MLR koala population projection model
############################################

## Authors: Corey JA Bradshaw, Katharina Peters & Frederik Saltre
## Date: January 2024

## Remove everything from the current R environment
rm(list = ls())  # Clears all objects from the workspace to avoid any conflicts

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Load necessary libraries
library(dplyr)   # For data manipulation
library(plotly)  # For creating interactive plots
source("matrixOperators.r")  # Sources external R script 'matrixOperators.r' that contains necessary functions


##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set functions

# Function to estimate the parameters of the beta distribution
# Inputs: mu (mean), var (variance)
# Outputs: alpha and beta parameters for the beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2  # Estimate alpha parameter
  beta <- alpha * (1 / mu - 1)  # Estimate beta parameter
  return(params = list(alpha = alpha, beta = beta))  # Return a list with alpha and beta
}

# Function to compute corrected Akaike Information Criterion (AICc) for multiple GLM models
# Inputs: One or more generalized linear models (GLMs)
# Outputs: AICc values for each model
AICc.glm <- function(...) {
  models <- list(...)  # Collect the input models into a list
  num.mod <- length(models)  # Count the number of models
  AICcs <- numeric(num.mod)  # Initialize a vector for AICc values
  ns <- numeric(num.mod)  # Initialize a vector for sample sizes
  ks <- numeric(num.mod)  # Initialize a vector for number of parameters
  AICc.vec <- rep(0,num.mod)  # Initialize the result vector for AICc values
  
  for (i in 1:num.mod) {  # Loop over each model
    # Calculate sample size (n) and number of parameters (k) for each model
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff)) + 1
    
    # Calculate AICc for the model
    AICcs[i] <- (-2 * logLik(models[[i]])) + ((2 * k * n) / (n - k - 1))  # AICc formula
    ns[i] <- n  # Store sample size
    ks[i] <- k  # Store number of parameters
    AICc.vec[i] <- AICcs[i]  # Store AICc for each model
  }
  return(AICc.vec)  # Return the vector of AICc values
}

# Helper function to extract the number of parameters (k) from a GLM
k.glm <- function(x) {
  if (length(x$df.residual) == 0) k <- sum(x$dims$ncol) else k <- (length(x$coeff) + 1)
}

delta.IC <- function(x) x - min(x)# Function to calculate the difference (delta) from the minimum value in a vector (x) of information criteria (IC)
weight.IC <- function(x) (exp(-0.5 * x)) / sum(exp(-0.5 * x))# Function to compute model weights from a vector (x) of delta ICs
chdev.glm <- function(x) (((as.numeric(x[12]) - as.numeric(x[10])) / as.numeric(x[12])) * 100)# Function to calculate percentage change in deviance for a GLM model

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Longevity parameter setup
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Defines parameters related to the koala's life span and age distribution for the population model. 
# 'longev' defines the maximum lifespan of koalas in years
longev <- 12  # Set koala longevity to 12 years (can adjust to 13 to see differences in predicted growth rates)
age.vec <- seq(0, longev, 1)  # Create an age vector from 0 to 'longev'
lage <- length(age.vec)  # Length of the age vector (number of age groups)
sex.ratio <- 0.5  # Assuming equal male-to-female ratio


##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Matrix construction
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Initializes the matrix used for modeling population transitions based on life stages.
stages <- lage  # Number of life stages corresponding to the age vector length
popmat <- matrix(0, nrow = stages, ncol = stages)  # Initialize a square matrix with dimensions based on life stages
colnames(popmat) <- age.vec[1:stages]  # Label matrix columns with age groups
rownames(popmat) <- age.vec[1:stages]  # Label matrix rows with age groups


##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Fertility data
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# based on Rhodes et al. (2011)
fert1.med <- c(0, 0.08, 0.47, 0.71, 0.71, 0.71, 0.71, 0.71, 0.71, 0.71, 0.71, 0.71, 0.71)  # Median fertility rates
fert1.lo <- c(0, 0.02, 0.30, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65)  # Lower bound fertility rates
fert1.up <- c(0, 0.19, 0.65, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75)  # Upper bound fertility rates
plot(age.vec, fert1.med, type = "l", xlab = "age", ylab = "proportion females with young", ylim = c(0, 1))# Plot fertility curves
lines(age.vec, fert1.lo, lty = 2, col = "red")  # Lower bounds (dashed red line)
lines(age.vec, fert1.up, lty = 2, col = "red")  # Upper bounds (dashed red line)


##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Survival data
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# female and male based on Rhodes et al. (2011)
survf1.med <- 1 - c(0.13, 0.20, 0.29, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25)  # Median survival rates for females
survf1.up <- 1 - c(0.08, 0.08, 0.16, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18)  # Upper bound survival rates for females
survf1.lo <- 1 - c(0.18, 0.36, 0.46, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34)  # Lower bound survival rates for females

survm1.med <- 1 - c(0.17, 0.23, 0.61, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34)  # Median survival rates for males
survm1.up <- 1 - c(0.12, 0.09, 0.38, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23)  # Upper bound survival rates for males
survm1.lo <- 1 - c(0.22, 0.44, 0.89, 0.46, 0.46, 0.46, 0.46, 0.46, 0.46, 0.46, 0.46, 0.46, 0.46)  # Lower bound survival rates for males

survf1.se.prop <- (((survf1.up - survf1.med) / survf1.med) + ((survf1.med - survf1.lo) / survf1.med)) / 2  # Proportional SE for females

# Plot survival curves for females and males
par(mfrow = c(1, 2))  # Set layout for 2 plots side by side
plot(age.vec, survf1.med, type = "l", xlab = "age", ylab = "survival", ylim = c(0, 1))
lines(age.vec, survf1.lo, lty = 2, col = "red")
lines(age.vec, survf1.up, lty = 2, col = "red")

plot(age.vec, survm1.med, type = "l", xlab = "age", ylab = "survival", ylim = c(0, 1))
lines(age.vec, survm1.lo, lty = 2, col = "red")
lines(age.vec, survm1.up, lty = 2, col = "red")
par(mfrow = c(1, 1))  # Reset to single plot layout


## Survival data for different koala populations from various studies (Oakey, Springsure, etc.)
# Penn et al. (2000) (note: SEs are all 10% of mean; Springsure males assumed to be same as females)
# based on Oakey
survf2.med <- 1 - c(0.33, 0.17, NA, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09)# Median survival rates for females
survf2.se <- c(0.033, 0.017, NA, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009)# SE survival rates for females
survm2.med <- 1 - c(0.20, 0.23, 0.23, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26) # Median survival rates for males
survm2.se <- c(0.02, 0.023, 0.023, 0.026, 0.026, 0.026, 0.026, 0.026, 0.026, 0.026, 0.026, 0.026, 0.026) # SE survival rates for males

survf2.se.prop <- (survf2.se/survf2.med)

# based on Springsure
survf3.med <- 1 - c(0.30, 0.16, NA, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08)# Median survival rates for females
survf3.se <- c(0.003, 0.002, NA, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008)# SE survival rates for females
survm3.med <- 1 - c(0.20, 0.23, 0.23, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26)# Median survival rates for males
survm3.se <- c(0.02, 0.023, 0.023, 0.026, 0.026, 0.026, 0.026, 0.026, 0.026, 0.026, 0.026, 0.026, 0.026)# SE survival rates for males

survf3.se.prop <- (survf3.se/survf3.med)

# Plot survival data for both populations
par(mfrow=c(2,2)) # Layout for 4 plots
plot(age.vec, survf2.med, type="l", xlab="age", ylab="survival", ylim=c(0,1))
lines(age.vec, survf2.med - 1.96*survf2.se, lty=2, col="red")
lines(age.vec, survf2.med + 1.96*survf2.se, lty=2, col="red")
plot(age.vec, survf3.med, type="l", xlab="age", ylab="survival", ylim=c(0,1))
lines(age.vec, survf3.med - 1.96*survf3.se, lty=2, col="red")
lines(age.vec, survf3.med + 1.96*survf3.se, lty=2, col="red")

plot(age.vec, survm2.med, type="l", xlab="age", ylab="survival", ylim=c(0,1))
lines(age.vec, survm2.med - 1.96*survm2.se, lty=2, col="red")
lines(age.vec, survm2.med + 1.96*survm2.se, lty=2, col="red")
plot(age.vec, survm3.med, type="l", xlab="age", ylab="survival", ylim=c(0,1))
lines(age.vec, survm3.med - 1.96*survm3.se, lty=2, col="red")
lines(age.vec, survm3.med + 1.96*survm3.se, lty=2, col="red")
par(mfrow=c(1,1))

# More survival data from Dique et al. (2003) and Lunney et al. (2007)
# Dique et al. (2003)
survf4.sed.med <- 1 - c(0.11, 0.06, 0.25, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
survf4.dsp.med <- 1 - c(1, 0.22, 0.5, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
survf4.bth.med <- (survf4.sed.med + survf4.dsp.med) / 2
  
survm4.sed.med <- 1 - c(0.17, 0.0, 0.0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
survm4.dsp.med <- 1 - c(1, 0.36, 0.0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
survm4.bth.med <- (survm4.sed.med + survm4.dsp.med) / 2

# Lunney et al. (2007)
survf5.med <- 1 - c(0.4, 0.4, NA, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23)
survf5.se <- c(0.04, 0.04, NA, 0.023, 0.023, 0.023, 0.023, 0.023, 0.023, 0.023, 0.023, 0.023, 0.023)
survm5.med <- 1 - c(0.4, 0.4, 0.4, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39)
survm5.se <- c(0.04, 0.004, 0.004, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039)

survf5.se.prop <- (survf5.se/survf5.med)

plot(age.vec, survf1.med, type="l", ylim=c(0,1))
lines(age.vec, survf2.med, lty=1)
lines(age.vec, survf3.med, lty=1)
lines(age.vec[1:3], survf4.bth.med[1:3], lty=1)
lines(age.vec, survf5.med, lty=1)

lines(age.vec, survf1.up, lty=2, col="red")
lines(age.vec, survf2.med+1.96*survf2.se, lty=2, col="red")
lines(age.vec, survf3.med + 1.96*survf3.se, lty=2, col="red")
lines(age.vec, survf5.med + 1.96*survf5.se, lty=2, col="red")

lines(age.vec, survf1.lo, lty=2, col="red")
lines(age.vec, survf2.med - 1.96*survf2.se, lty=2, col="red")
lines(age.vec, survf3.med - 1.96*survf3.se, lty=2, col="red")
lines(age.vec, survf5.med - 1.96*survf5.se, lty=2, col="red")


survfmn.dat <- data.frame(survf1.med, survf2.med, survf3.med, survf4.bth.med, survf5.med)
survf.mn <- apply(survfmn.dat, 1, median, na.rm=T)
survfsepr.dat <- data.frame(survf1.se.prop, survf2.se.prop, survf3.se.prop, survf5.se.prop)
survfsepr.mn <- apply(survfsepr.dat, 1, median, na.rm=T)

survflo.dat <- data.frame(survf1.lo, survf2.med - 1.96*survf2.se, survf3.med - 1.96*survf3.se, survf5.med - 1.96*survf5.se)
colnames(survflo.dat) <- c("s1lo", "s2lo", "s3lo", "s4lo")
survfup.dat <- data.frame(survf1.up, survf2.med + 1.96*survf2.se, survf3.med + 1.96*survf3.se, survf5.med + 1.96*survf5.se)
colnames(survfup.dat) <- c("s1up", "s2up", "s3up", "s4up")

plot(age.vec, survf.mn, type="l", ylim=c(0,1))
lines(age.vec, survf.mn - 1.96*survfsepr.mn, lty=2, col="red")
lines(age.vec, ifelse((survf.mn + 1.96*survfsepr.mn) > 1, 1, survf.mn + 1.96*survfsepr.mn), lty=2, col="red")



# fit function to mean survival to smooth
# Fits an exponential model to smooth the survival curve using the nls (nonlinear least squares) method.
# exponential association
# y = a * exp(-exp(b-c*x))
surv.dat <- data.frame(age.vec, c(0.70,0.80,0.84,0.85,0.86,0.87,0.89,rep(0.9,6)))
colnames(surv.dat) <- c("agemod", "survmod")
param.init <- c(8.6e-01, -1.84e00, 3.06e-01)
fit.expa <- nls(survmod ~ a1 * exp(-exp(b2 - c3 * agemod)), 
                data = surv.dat,
                algorithm = "default",
                start = c(a1 = param.init[1], b2 = param.init[2], c3 = param.init[3]),
                trace = TRUE,      
                nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))

# Plot smoothed survival curve
plot(age.vec, survf.mn, pch=19, xlab="age", ylab="survival", ylim=c(0.5, 0.95),col="white")
age.pred.vec <- seq(0,12,1)
pred.surv <- as.numeric(coef(fit.expa)[1]) * exp(-exp(as.numeric(coef(fit.expa)[2]) - as.numeric(coef(fit.expa)[3]) * age.pred.vec))
lines(age.pred.vec, pred.surv, lty=2, lwd=1, col="red")
points(age.vec, survfmn.dat[,1],pch=19,cex=0.5,col="light grey")
points(age.vec, survfmn.dat[,2],pch=19,cex=0.5,col="light grey")
points(age.vec, survfmn.dat[,3],pch=19,cex=0.5,col="light grey")
points(age.vec, survfmn.dat[,4],pch=19,cex=0.5,col="light grey")
points(age.vec, survfmn.dat[,5],pch=19,cex=0.5,col="light grey")
points(age.vec, survflo.dat[,1],pch=19,cex=0.5,col="pink")
points(age.vec, survflo.dat[,2],pch=19,cex=0.5,col="pink")
points(age.vec, survflo.dat[,3],pch=19,cex=0.5,col="pink")
points(age.vec, survflo.dat[,4],pch=19,cex=0.5,col="pink")
points(age.vec, survfup.dat[,1],pch=19,cex=0.5,col="light green")
points(age.vec, survfup.dat[,2],pch=19,cex=0.5,col="light green")
points(age.vec, survfup.dat[,3],pch=19,cex=0.5,col="light green")
points(age.vec, survfup.dat[,4],pch=19,cex=0.5,col="light green")


##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Beta-resampled survival simulation across ages
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Create survival standard deviation (SD) dataframe
# The SD is calculated by dividing the difference between upper and lower survival bounds by 1.96*2
# survfsd.dat calculates the standard deviation (SD) of survival rates for each age based on the difference between upper and lower bounds.
survfsd.dat <- (survfup.dat - survflo.dat) / (1.96*2)
survfsd.dat[3,2:4] <- survfsd.dat[3,1]
survfsd.dat# View the updated SD dataframe

# Interpolate missing values in `survfmn.dat` (mean survival rates) by taking the average of neighboring values
# Specifically, we interpolate the missing value in the 3rd row of columns 2 and 3
survfmn.dat[3,2] <- mean(c(survfmn.dat[2,2], survfmn.dat[4,2]))# Interpolate for column 2
survfmn.dat[3,3] <- mean(c(survfmn.dat[2,3], survfmn.dat[4,3]))# Interpolate for column 3
survfmn.dat# View the updated survival means dataframe

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Simulation setup
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The code sets up a survival simulation with 10,000 iterations. 
# A matrix surv.it.mat is initialized to store simulated survival rates for each iteration, 
# and vectors are prepared to store the parameters from fitting an exponential survival model (fit.expa).
surv.iter <- 10000  # Number of simulation iterations
itdivs <- surv.iter / 10  # For possible progress tracking or dividing iterations into chunks
survfmn.sim.dat <- survfmn.dat  # Copy of survival mean data to be used in the simulation
age.pred.vec <- seq(0, 12, 1)  # Predicted age sequence from 0 to 12 years

# Initialize a matrix to store survival rates for each iteration
# Rows represent iterations, and columns represent survival rates across ages
surv.it.mat <- matrix(data=0, nrow=surv.iter, ncol=dim(survfmn.sim.dat)[1])
fit.expa.list <- vector("list", surv.iter)# Initialize a list to store the exponential model fits from each iteration
fit.exp.a1.vec <- fit.exp.b1.vec <- fit.exp.c1.vec <- rep(0,surv.iter)# Vectors to store the fitted parameters (a1, b1, c1) from the nonlinear model across iterations

# Loop through iterations for beta-resampled survival simulation
# Key Points of the Loop:
#  1) Random Sampling: Each iteration samples survival standard deviations (surv.sd.ran) and survival means (surv.mn.ran) from the given data.
#                      Beta distribution parameters (alpha and beta) are estimated based on the sampled survival rates.
#  2) Beta Resampling:For each age group, the survival rate is sampled from a beta distribution with the estimated parameters.
#  3) Nonlinear Model Fitting: An exponential association model is fitted to the sampled survival data using nls.If the model fitting fails for an iteration, the parameters are set to NA.
#  4) Progress Monitoring: Every itdivs iterations (10% of total), the progress is printed to the console for tracking purposes.
for (k in 1:surv.iter) {
  # Randomly sample standard deviations for each age group
  surv.sd.ran <- apply(survfsd.dat, MARGIN=1, sample, 1, replace=T) # Sampling 1 SD per row
  survfmn.sim.dat[4:13, 4] <- apply(survfmn.dat[4:13, c(1:3, 5)], MARGIN = 1, sample, 1, replace = TRUE)  # Randomly fill in missing data for ages 4 to 13 in the 4th column by sampling from relevant columns
  survfmn.sim.dat[3, 5] <- sample(survfmn.sim.dat[3, 1:4], 1, replace = TRUE)  # Randomly sample for the missing value in row 3, column 5 from columns 1 to 4 in the same row
  survfmn.sim.dat# Display the updated simulated mean survival data (optional; you can remove this in production)
  
  # Estimate the alpha and beta parameters for the beta distribution using the sampled mean and SD
  s.alpha <- estBetaParams(surv.mn.ran, surv.sd.ran^2)$alpha  # Alpha for beta distribution
  s.beta <- estBetaParams(surv.mn.ran, surv.sd.ran^2)$beta    # Beta for beta distribution
  
  surv.it <- rep(0, dim(survfmn.sim.dat)[1])   # Initialize a vector to store the survival rate for each age group for this iteration
  # Loop through each age group and sample survival rates from the beta distribution
  for (x in 1:dim(survfmn.sim.dat)[1]) {
    surv.it[x] <- rbeta(1, s.alpha[x], s.beta[x])  # Sample from beta distribution
  }
  surv.it.mat[k, ] <- surv.it  # Store the sampled survival rates for this iteration in the matrix

  
  # Fit an exponential survival function to smooth the sampled survival rates
  # Model: y = a * exp(-exp(b - c * x))
  surv.dat2 <- data.frame(age.vec, surv.it.mat[k, ])  # Create dataframe with age and survival rates
  colnames(surv.dat2) <- c("agemod", "survmod")      # Rename columns for model fitting
  
  # Initial parameter estimates for the exponential model
  param.init <- c(8.6e-01, -1.84e00, 3.06e-01)
  
  # Fit the nonlinear model using the nls function, and store the result in the fit list
  fit.expa.list[[k]] <- try(fit.expa <- nls(
    survmod ~ a1 * exp(-exp(b2 - c3 * agemod)),  # Exponential model formula
    data = surv.dat2,                            # Data to fit
    algorithm = "default",                       # Default algorithm
    start = c(a1 = param.init[1], b2 = param.init[2], c3 = param.init[3]),  # Starting parameters
    trace = FALSE,                               # Disable trace
    nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1 / 1024)),  # Control parameters
    silent = TRUE)  # Run silently to avoid errors interrupting loop
  
  # Check if the model fit was successful (no error)
  if (is.null(attr(fit.expa.list[[k]], 'condition')[1]$message) == TRUE) {
    # If successful, extract the fitted coefficients (a1, b2, c3) and store in their respective vectors
    fit.exp.a1.vec[k] <- as.numeric(coef(fit.expa.list[[k]])[1])  # Store a1 parameter
    fit.exp.b1.vec[k] <- as.numeric(coef(fit.expa.list[[k]])[2])  # Store b2 parameter
    fit.exp.c1.vec[k] <- as.numeric(coef(fit.expa.list[[k]])[3])  # Store c3 parameter
  } else {
    # If model fit fails, store NA in the coefficient vectors for this iteration
    fit.exp.a1.vec[k] <- fit.exp.b1.vec[k] <- fit.exp.c1.vec[k] <- NA
  }
  
  # Print progress at every 'itdivs' interval (e.g., every 1000 iterations if 'surv.iter' = 10000)
  if (k %% itdivs == 0) print(k)
} # end k loop

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Analyzing and plotting the simulation results
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
coef.dat <- data.frame(fit.exp.a1.vec, fit.exp.b1.vec, fit.exp.c1.vec)# Create a dataframe to store the coefficients from the model fits across iterations
coef.mn <- apply(coef.dat, MARGIN = 2, median, na.rm = TRUE)# Calculate the median of the coefficients across all iterations (ignoring NA values)
surv.sim.pred <- coef.mn[1] * exp(-exp(coef.mn[2] - coef.mn[3] * age.pred.vec))# Use the median coefficients to predict survival across ages (using the fitted exponential model)
surv.mn.calc <- apply(surv.it.mat, MARGIN = 2, mean, na.rm = TRUE)# Calculate the mean survival rate across iterations for each age group
# Calculate the 2.5th and 97.5th percentiles to get the 95% credible interval for survival rates
surv.lo.calc <- apply(surv.it.mat, MARGIN = 2, quantile, probs = 0.025, na.rm = TRUE)  # Lower bound
surv.up.calc <- apply(surv.it.mat, MARGIN = 2, quantile, probs = 0.975, na.rm = TRUE)  # Upper bound
# Plot the mean simulated survival rates across ages
plot(age.vec, surv.mn.calc, type = "l", 
     xlab = "age (years)", ylab = "simulated survival", 
     ylim = c(min(surv.lo.calc), max(surv.up.calc)))  # Set Y-axis limits to the 95% interval
# Add the predicted survival line based on median coefficients (from the fitted model)
lines(age.vec, surv.sim.pred, lty = 2, lwd = 3)  # Dashed line for the predicted survival
# Add the 95% credible interval (shaded region)
lines(age.vec, surv.lo.calc, lty = 2, col = "pink")  # Lower bound (2.5th percentile)
lines(age.vec, surv.up.calc, lty = 2, col = "pink")  # Upper bound (97.5th percentile)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Population matrix construction and analysis
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set the diagonal of the population matrix (from row 2 onward) to the predicted survival rates
diag(popmat[2:(stages), ]) <- surv.sim.pred[-stages]  # For ages 1 to longev-1

# Set the survival rate of the oldest age group (last row and column of the matrix)
popmat[stages, stages] <- surv.sim.pred[stages]  # Survival rate for the last age class
popmat[1, ] <- fert1.up * sex.ratio# Set the first row (representing fertility) to the upper bound fertility rates multiplied by the sex ratio
popmat.orig <- popmat# Save the original population matrix (before modifications)
popmat# View the modified population matrix


##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Population matrix properties and demographic calculations
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Calculate and analyze key properties of the population matrix
# Calculate the maximum annual growth rate (lambda), also known as the dominant eigenvalue
max.lambda(popmat)  # 1-year lambda (annual population growth rate)
# Calculate the intrinsic rate of population increase (r), which is the natural log of lambda
max.r(popmat)  # Intrinsic rate of population growth (1-year r)
stable.stage.dist(popmat)# Calculate the stable stage distribution (the long-term proportion of individuals in each stage)
R.val(popmat, stages)# Calculate the reproductive value of each stage (expected future reproduction by an individual in each stage)
gen.l <- G.val(popmat, stages)# Calculate the mean generation length (average age of parents in a stable population)

# Calculate the probability of a catastrophe (Reed et al. 2003) based on generation length
cat.pr <- 0.14 / gen.l  # Probability of catastrophe
cat.mort <- 0.5  # Assumed mortality due to catastrophe

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Calculating rmax (maximum intrinsic growth rate) assuming perfect survival (Sx = 1)
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create a version of the population matrix where survival is set to 1 (maximum possible survival)
popmat.max <- popmat.orig  # Start with the original matrix
# Set survival to 1 for all stages except the first (fertility) row
diag(popmat.max[2:(stages), ]) <- 1
popmat.max[stages, stages] <- 1  # Set the survival for the last stage to 1
popmat.max# View the modified population matrix with maximum survival
max.r(popmat.max)# Calculate the intrinsic growth rate (r) assuming maximum survival

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Allometric prediction of rmax using Hennemann's formula (1983)
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Use Hennemann's allometric equation to predict rmax based on body mass
# Formula: r.max = 10^(0.6914 - 0.2622 * log10(M)), where M is body mass in grams
mass.koala <- 6.2  # Koala body mass in kg (reference: https://link.springer.com/article/10.1007/s00265-010-1136-4)
# Predict the maximum r (r.max) based on body mass
r.max.pred <- 10^(0.6914 - 0.2622 * (log10(mass.koala * 1000)))  # Convert kg to grams for the formula
r.max.pred# Display the predicted r.max value
r.max.pred / max.r(popmat.max)# Calculate the predicted maximum growth rate as a fraction of the intrinsic growth rate from the population matrix
10^r.max.pred / max.lambda(popmat.max)# Similarly, calculate the ratio of the predicted growth rate to the maximum lambda (annual growth rate)


##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Initial population vector setup and visualization
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# From the point pattern model, we estimated a median total population of 32,851 koalas (32231 - 38119).
# The current population size is split between males and females (sex ratio applied)
curr.pop <- 32851 * sex.ratio  # Multiply the population estimate by the sex ratio to get the female population
# Create the initial population vector using the stable stage distribution
init.vec <- curr.pop * stable.stage.dist(popmat)
# Plot the initial population distribution across ages
plot(age.vec, init.vec, xlab = "age (years)", ylab = "Nf", type = "l")  # Line plot of the female population across ages


#######################################################################################################################
## Population Projection Model ########################################################################################
#######################################################################################################################
# Set time parameters for the projection
yr.now <- 2016  # Current year (starting point of the projection)
yr.end <- 2040  # End year for the projection
t <- (yr.end - yr.now)  # Total number of years to project (e.g., 24 years)

# Use the original population matrix (before any modifications)
popmat <- popmat.orig

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Population Storage Setup
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create a matrix to store population counts for each stage over time
# Rows represent age stages, and columns represent time (years)
n.mat <- matrix(0, nrow = stages, ncol = (t + 1))  # Initialize with zeros

# Set the initial population vector (from stable stage distribution) as the first column in n.mat
n.mat[, 1] <- init.vec  # This represents the population at year 2016

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Population Projection Loop
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Loop through each year and project the population using matrix multiplication
for (i in 1:t) {
  n.mat[, i + 1] <- popmat %*% n.mat[, i]  # Multiply population matrix by current population vector
}

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Plotting the Projection Results
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## The resulting plot shows how the total population of females evolves over the projection period, 
## based on the fertility, survival rates, and population dynamics modeled in the population matrix.
# Create a vector of years for the x-axis of the plot
yrs <- seq(yr.now, yr.end, 1)  # Sequence of years from 2016 to 2040

# Plot the total population (sum across all age groups) over time
plot(yrs, as.vector(colSums(n.mat)), type = "l", 
     xlab = "year", ylab = "Nf", 
     xlim = c(yr.now, yr.end))  # Set x-axis limits to the projection period


#####################################################################################
## SCENARIO 0: 
## stochastic projection for UNMANAGED population SCENARIO
## no sterilisation scenarion 
## used as a baseline
#####################################################################################
# Key Concepts:
#  1) Stochasticity: This projection adds stochastic variability to fertility and survival rates and allows for possible catastrophic events (e.g., mass mortality due to environmental factors).
#  2) Density Dependence: Survival rates are adjusted based on population density, meaning that as population size increases, survival probability decreases.
#  3) Catastrophes: Catastrophic events can randomly occur, reducing survival rates significantly.

## Initial Setup
sex.ratio <- 0.5  # Assuming a 50% female population

start.pop <- curr.pop  # Start population size
init.vec <- start.pop * stable.stage.dist(popmat.orig)  # Initial population vector using stable stage distribution

catastrophe.used <- 1  # Set to 1 to enable stochastic catastrophes

yr.now <- 2016  # Start year
yr.end <- 2040  # End year for projection

t <- (yr.end - yr.now)  # Total number of years to project
yrs <- seq(yr.now, yr.end, 1)  # Create a year vector from 2016 to 2040

# Density feedback on survival (based on Todd et al. 2008)
K.pop <- 1.2 * curr.pop  # Carrying capacity (e.g., 1.2 times the current population)

# Survival multipliers for different population densities
surv.mult.start <- 1.0
surv.mult.lo <- 0.9975
surv.mult.mid <- 0.994
surv.mult.end <- 0.985

K.lo <- 1.05 * curr.pop
K.mid <- 1.1 * curr.pop

# Survival multipliers based on population sizes
K.vec <- c(curr.pop,K.lo,K.mid,K.pop)
surv.mult.vec <- c(surv.mult.start, surv.mult.lo, surv.mult.mid, surv.mult.end)
plot(K.vec, surv.mult.vec, pch=19)# Plot the relationship between population size and survival multipliers

## Fitting logistic function to survival multipliers based on population size
# Logistic function: y = a / (exp((b - X) / c) + 1)
library(deSolve)
DD.dat <- data.frame(K.vec, surv.mult.vec)  # Create data frame of carrying capacities and survival multipliers

# Initial parameter estimation for logistic model
param.init.est <- as.numeric(getInitial(surv.mult.vec ~ SSlogis(K.vec, a, b, c), dat=DD.dat))
param.init.est

# Fit a logistic function to model survival multipliers based on population size
fit.expd <- nls(
  surv.mult.vec ~ a / (exp((b - K.vec) / c) + 1), 
  data = DD.dat,
  algorithm = "port",
  start = c(a = param.init.est[1], b = param.init.est[2], c = param.init.est[3]),
  trace = TRUE,      
  nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024)
)

# Plot fitted logistic function for survival multipliers
plot(K.vec, surv.mult.vec, pch = 19, xlab = "Nf", ylab = "Reduction in survival")
K.pred.vec <- seq(curr.pop, K.pop, 10)  # Generate sequence of population sizes
pred.surv.mult <- coef(fit.expd)[1] / (exp((coef(fit.expd)[2] - K.pred.vec) / coef(fit.expd)[3]) + 1)
lines(K.pred.vec, pred.surv.mult, lty = 2, lwd = 1, col = "red")

## Population Storage Setup
n.mat <- matrix(0, nrow = stages, ncol = (t + 1))  # Initialize matrix for population counts
n.mat[, 1] <- init.vec  # Set initial population vector
popmat <- popmat.orig  # Use the original population matrix

## Stochastic Projection Loop
iter <- 1000  # Number of iterations (1000 simulations)
itdiv <- iter / 10  # For printing progress

# Storage for total population size and growth rates
n.sums.mat <- matrix(0, nrow = iter, ncol = (t + 1))  # Stores total population size per year for each iteration
r.mat <- matrix(0, nrow = iter, ncol = t)  # Stores stochastic growth rates (r) for each iteration

for (e in 1:iter) {
#  account for randomness in the population projections. 
#  This is typically referred to as Monte Carlo simulations in stochastic models. 
#  Each iteration (e) simulates one possible population trajectory over the projection period (from 2016 to 2040) 
#  by incorporating random variation (stochasticity) in fertility, survival, and potential catastrophes
  r.stoch <- rep(0, t)  # Reset stochastic growth rate vector for each iteration
  
  # Reset population matrix to original
  popmat <- popmat.orig
  
  # Projection Loop: project population year-by-year
  for (i in 1:t) {
    
    ## Fertility Beta Sampler
    fert.sd <- 0.03  # Assumed standard deviation for fertility
    fert.alpha <- estBetaParams(fert1.up, fert.sd^2)$alpha  # Beta distribution alpha for fertility
    fert.beta <- estBetaParams(fert1.up, fert.sd^2)$beta  # Beta distribution beta for fertility
    fert.alpha[1] <- fert.beta[1] <- 0  # Set fertility for age 0 to 0
    fert.stoch.calc <- sapply(1:stages, function(x) rbeta(1, fert.alpha[x], fert.beta[x]))  # Sample fertility
    fert.stoch <- ifelse(is.na(fert.stoch.calc), 0, fert.stoch.calc) * sex.ratio  # Adjust fertility by sex ratio
    fert.stoch[1] <- 0  # No reproduction for age 0
    
    ## Density-Dependent Survival Multiplier
    surv.mult.a <- coef(fit.expd)[1] / (exp((coef(fit.expd)[2] - sum(n.mat[, i])) / coef(fit.expd)[3]) + 1)
    surv.mult <- ifelse(round(sum(n.mat[, i]), 0) <= curr.pop, 1, surv.mult.a)  # Apply density feedback to survival
    surv.it <- surv.mult * surv.sim.pred  # Adjust survival rates by density
    
    ## Survival Beta Sampler
    Sx.sd <- 0.05  # Standard deviation for survival
    Sx.alpha <- estBetaParams(surv.it, Sx.sd^2)$alpha  # Beta distribution alpha for survival
    Sx.beta <- estBetaParams(surv.it, Sx.sd^2)$beta  # Beta distribution beta for survival
    Sx.stoch <- sapply(1:stages, function(x) rbeta(1, Sx.alpha[x], Sx.beta[x]))  # Sample survival

    ## Reconstruct Population Matrix with Stochastic Fertility and Survival
    popmat[1, ] <- fert.stoch  # Set fertility in the matrix
    diag(popmat[2:stages, ]) <- Sx.stoch[-stages]  # Set survival in the matrix
    
    ## Stochastic Catastrophes
    if (catastrophe.used == 1) {
      catastrophe <- rbinom(1, 1, cat.pr)  # Randomly determine if a catastrophe occurs
      if (catastrophe == 1) {
        cat.mort.stoch <- rbeta(1, estBetaParams(cat.mort, 0.05^2)$alpha, estBetaParams(cat.mort, 0.05^2)$beta)  # Catastrophic mortality rate
        diag(popmat[2:stages, ]) <- cat.mort.stoch * Sx.stoch[-stages]  # Adjust survival for catastrophe
      }
    }
    
    ## Project Population for the Next Year
    n.mat[, i + 1] <- popmat %*% n.mat[, i]  # Multiply population matrix by current population
    
    ## Save Stochastic Growth Rate (r)
    r.running <- log(sum(n.mat[, i + 1]) / sum(n.mat[, i]))  # Calculate growth rate (r)
    r.stoch[i] <- ifelse(r.running == -Inf, NA, r.running)  # Handle cases where r is -Inf
  }
  
  r.mat[e, ] <- r.stoch  # Store growth rates for this iteration
  n.sums.mat[e, ] <- colSums(n.mat)  # Store total population size for this iteration
  
  if (e %% itdiv == 0) print(e)  # Print progress every 10% of iterations
}


## Calculate Confidence Intervals for Population Size
# This loop calculates the median population size as well as the 95% confidence interval (2.5th and 97.5th percentiles) for each year of the projection.
n.up <- n.lo <- n.mn <- rep(0,(t+1))
for (q in 1:(t+1)) {
  n.mn[q] <- median(n.sums.mat[,q], na.rm=T) # assume equal sex ratio
  n.up[q] <- as.vector(quantile(n.sums.mat[,q], probs=0.975, na.rm=T))
  n.lo[q] <- as.vector(quantile(n.sums.mat[,q], probs=0.025, na.rm=T))
}
plot(yrs, n.mn, type="l", xlab="year",ylab="Nf",xlim=c(yr.now,yr.end),ylim=c(min(n.lo),max(n.up)))
lines(yrs, n.up, lty=2, col="red")
lines(yrs, n.lo, lty=2, col="red")

#The median, 2.5th, and 97.5th percentiles of the overall median population size are calculated and printed.
n.mn.mn <- median(n.mn, n.rm=T)
n.mn.lo <- quantile(n.mn, probs=0.025, n.rm=T)
n.mn.up <- quantile(n.mn, probs=0.975, n.rm=T)
print(n.mn.mn)# Print the median population size over time

# This section plots the population trajectories for the first 100 iterations to visualize the variability across different simulation runs.
# Each line represents a different simulation, showing the spread of possible population outcomes over time.
plot(yrs, n.sums.mat[1,],type="l",lty=2,lwd=0.2,xlab="",ylab="Nf", ylim=c(min(n.sums.mat), max(n.sums.mat)))
for (a in 2:100) {
  lines(yrs, n.sums.mat[a,],lty=2,lwd=0.2)
}
lines(yrs, n.mn, lwd=3, lty=2)

#The ndat.out data frame stores the years (yrs), median population (n.mn), and confidence intervals (n.lo and n.up) for further analysis or export.
ndat.out <- data.frame(yrs, n.mn, n.lo, n.up)
colnames(ndat.out) <- c("yrs","median","lowCI","upCI")
write.table(ndat.out, file="Unmanaged_Scenario(Median).csv",sep=",", dec = ".", row.names = F, col.names = T)

stats.ndat.out<-c(mean(n.mn, na.rm=T),quantile(n.mn, probs=0.025, na.rm=T),quantile(n.mn, probs=0.975, na.rm=T))
stats.ndat.out*2 #multiply by 2 for the full population


#####################################################################################
## Scenario 1 & 2
## stochastic projection for STERILIZED female ONLY (1) or FEMALE + DAUGHTER (2)
## incrementing proportion of individual sterilised
## with costs assumed from Prowse & Delean (2019)
#####################################################################################
# stochastic population projection model that includes sterilization as a management strategy. 
# It allows you to explore different sterilization rates, costs, and the effects on the koala population over time. 
# The cost of capturing and sterilizing koalas is modeled as a function of population density, and various results such as population size and costs are stored for analysis.
## Set assumed sex ratio (proportion of females in the population)
sex.ratio <- 0.5 # 50% female population

## Choose sterilization strategy:
## ster.which = 1 : only adult females sterilized
## ster.which = 2 : adult females + their female offspring sterilized
## ster.which = 0 : no sterilization
ster.which <- 1  # In this case, only adult females are sterilized

## Set Minimum Viable Population (MVP)
mvp <- 100 * sex.ratio  # MVP is 100 individuals, adjusted by the sex ratio (50 females)


## Initial population setup
start.pop <- curr.pop  # Start population size
init.vec <- start.pop * stable.stage.dist(popmat.orig)  # Initial population vector using stable stage distribution
popmat <- popmat.orig  # Use the original population matrix

## Enable catastrophic events (set to 0 if not used)
catastrophe.used <- 1  # Set to 1 to enable catastrophes in the projection

## Time projection setup
yr.now <- 2016  # Start year
yr.end <- 2040  # End year of the projection
t <- (yr.end - yr.now)  # Number of years in the projection (24 years)
yrs <- seq(yr.now, yr.end, 1)  # Create a sequence of years from 2016 to 2040

##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Cost Function Setup
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Area for the population (in hectares)
area <- 3080 * 100  # Set to 3080 ha (assuming MLR region)

## Set population density thresholds and corresponding search times
max.D <- K.pop / area  # Maximum population density (carrying capacity divided by area)
lolo.D <- 0.02  # Very low population density
lo.D <- 0.04  # Low population density
mid.D <- 0.08  # Medium population density
min.D <- mvp / area  # Minimum population density, based on MVP
srch.max <- 150  # Maximum search time (in minutes) at very low density
srch.lolo <- 90  # Search time at very low density
srch.lo <- 60  # Search time at low density
srch.mid <- 40  # Search time at medium density
srch.min <- 25  # Minimum search time at higher density

## Store the population density and search time values in vectors
D.vec <- c(min.D, lolo.D, lo.D, mid.D, max.D)  # Population densities
srch.vec <- c(srch.max, srch.lolo, srch.lo, srch.mid, srch.min)  # Search times in minutes

## Plot the relationship between population density and search time
plot(D.vec, srch.vec, pch = 19, xlab = "density", ylab = "search time (minutes)")

## Fit a nonlinear model to estimate the relationship between population density and search time
# Model: y = a * (x - b)^c
Dsrch.dat <- data.frame(D.vec, srch.vec)  # Create a data frame of population density and search time
colnames(Dsrch.dat) <- c("D", "srch")  # Rename columns for clarity
param.init <- c(0.25, -0.13, -3)  # Initial guesses for the parameters of the model
# Fit a nonlinear least squares (nls) model to the data
fit.sexpd <- nls(srch ~ a * (D.vec - b)^c, 
                 data = Dsrch.dat,
                 algorithm = "port",
                 start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1 / 1024))
## Plot the fitted model
plot(D.vec, srch.vec, pch = 19, xlab = "density", ylab = "search time (minutes)")
D.pred.vec <- seq(min.D, max.D * 1.25, 0.001)  # Create a sequence of population densities for prediction
pred.srch <- as.numeric(coef(fit.sexpd)[1]) * (D.pred.vec - as.numeric(coef(fit.sexpd)[2]))^as.numeric(coef(fit.sexpd)[3])
lines(D.pred.vec, pred.srch, lty = 2, lwd = 1, col = "red")  # Add fitted line to the plot


##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Cost Parameters Setup
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cph <- 30  # Cost per hour ($)
hpkc <- 0.83  # Hours per koala capture
spk <- 27  # Sterilization cost per koala

##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Iteration Setup for Population Projection
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
iter <- 10000  # Number of iterations (simulations) to run
itdiv <- iter / 10  # Print progress every 10% of iterations

## Set up a range of sterilization rates to explore (0% to 90% in 10% increments)
ster.vec <- seq(0, 0.9, 0.1)  # Sterilization rates to test (from 0% to 90%)

## Storage variables for population projection results
n.mn.mn <- n.mn.lo <- n.mn.up <- prop.mn <- prop.lo <- prop.up <- 
  min.mn <- min.lo <- min.up <- p.qext.max <- 
  nster.mn <- nster.up <- nster.lo <- cost.mn <- cost.up <- cost.lo <- 
  rep(0, length(ster.vec))  # Initialize vectors to store median, lower, and upper limits of population metrics


## Create matrices to store annual sterilization costs for each sterilization rate
cost.yr.mn.mat <- cost.yr.lo.mat <- cost.yr.up.mat <- matrix(0, nrow = length(ster.vec), ncol = t + 1)

for (s in 1:length(ster.vec)) {
  # This outer loop iterates over different sterilization proportions (ster.vec), testing how varying the proportion of sterilized females impacts the population.
  ## Initialize storage matrices for population and sterilization calculations
  n.mat <- matrix(0, nrow = stages, ncol = (t+1))  # Matrix for population projection over time
  n.mat[,1] <- init.vec  # Set initial population based on stable stage distribution
  popmat <- popmat.orig  # Use original population matrix for projection
  
  # Storage for population summaries and costs across iterations
  n.sums.mat <- nad.sums.mat <- njv.sums.mat <- cost.mat <- nad.ster.mat <- njv.ster.mat <- matrix(data = 0, nrow = iter, ncol = (t+1))
  
  ## Reset population matrix to original values before each iteration
  popmat <- popmat.orig
  
  ## Loop over each simulation iteration (stochastic projection)
  for (e in 1:iter) {
    #This middle loop runs iter (e.g., 1000) iterations of the population projection. Each iteration introduces stochastic variability in fertility, survival, and catastrophes.
    ## Reset population matrix to original values at the start of each iteration
    popmat <- popmat.orig
    
    ## Loop over each year in the projection
    for (i in 1:t) {
      #This inner loop projects the population for each year (t years), adjusting fertility and survival stochastically.
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Stochastic Fertility Sampling (Beta Distribution)
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      fert.sd <- 0.03  # Standard deviation for fertility sampling
      fert.alpha1 <- estBetaParams(fert1.up, fert.sd^2)$alpha  # Calculate alpha for beta distribution
      fert.alpha <- ifelse(is.na(fert.alpha1), 0, fert.alpha1)  # Handle missing values
      fert.beta1 <- estBetaParams(fert1.up, fert.sd^2)$beta  # Calculate beta for beta distribution
      fert.beta <- ifelse(is.na(fert.beta1), 0, fert.beta1)
      fert.alpha[1] <- fert.beta[1] <- 0  # Set fertility for age 0 to 0
      fert.stoch.calc <- rep(0, stages)  # Initialize stochastic fertility values
      
      
      # Loop over each stage and sample fertility
      for (x in 1:stages) {
        fert.stoch.if <- rbeta(1, fert.alpha[x], fert.beta[x])# Sample from beta distribution
        if (is.na(fert.stoch.if) == F) {# Handle missing values
          fert.stoch.calc[x] <- fert.stoch.if
        } # end if
        if (is.na(fert.stoch.if) == T) {
          fert.stoch.calc[x] <- 0
        }
      } # end x loop
      fert.stoch.calc[1] <- 0# Set fertility of newborns to 0
      fert.stoch.calc
      
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Stochastic Sterilization Proportion (Beta Distribution)
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (ster.vec[s] > 0 & ster.which > 0) {  # If sterilization is applied
        ster.sd <- 0.05  # Standard deviation for sterilization proportion
        ster.alpha <- estBetaParams(ster.vec[s], ster.sd^2)$alpha  # Alpha for beta distribution
        ster.beta <- estBetaParams(ster.vec[s], ster.sd^2)$beta  # Beta for beta distribution
        ster.stoch <- rbeta(1, ster.alpha, ster.beta)  # Sample sterilization proportion
        
        # Adjust sterilization for adult females and possibly their offspring
        if (ster.which == 2) {  # If sterilizing adult females + offspring
          ster.stoch.it1 <- ifelse((1 - (ster.stoch + ster.stoch * sex.ratio)) < 0, 0, 1 - (ster.stoch + ster.stoch * sex.ratio))
        }
        if (ster.which == 1) {  # If sterilizing only adult females
          ster.stoch.it1 <- ifelse((1 - ster.stoch) < 0, 0, 1 - ster.stoch)
        }
        ster.stoch.it <- ifelse(ster.stoch.it1 > 1, 1, ster.stoch.it1)
        fert.stoch.red <- fert.stoch.calc * ster.stoch.it  # Reduce fertility based on sterilization
        fert.stoch <- fert.stoch.red * sex.ratio  # Adjust fertility by sex ratio
      }
      
      # If no sterilization
      if (ster.vec[s] == 0) {
        ster.stoch <- 0  # No sterilization
        fert.stoch <- fert.stoch.calc * sex.ratio  # Set fertility based on sex ratio
      }
      
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Density-Dependent Survival Adjustment
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      surv.mult.a <- coef(fit.expd)[1] / ((exp((coef(fit.expd)[2] - sum(n.mat[, i])) / coef(fit.expd)[3])) + 1)
      surv.mult <- ifelse(sum(n.mat[, i]) <= curr.pop, 1, surv.mult.a)  # Adjust survival based on population density
      surv.it <- surv.mult * surv.sim.pred  # Adjust survival predictions
      
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Stochastic Survival Sampling (Beta Distribution)
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Sx.sd <- 0.05  # Standard deviation for survival
      Sx.alpha <- estBetaParams(surv.it, Sx.sd^2)$alpha  # Alpha for beta distribution
      Sx.beta <- estBetaParams(surv.it, Sx.sd^2)$beta  # Beta for beta distribution
      Sx.stoch <- rep(0, stages)  # Initialize stochastic survival values
      
      # Sample survival rates for each stage
      for (x in 1:stages) {
        Sx.stoch[x] <- rbeta(1, Sx.alpha[x], Sx.beta[x])  # Sample from beta distribution
      }
      
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Reconstruct Population Matrix with Stochastic Fertility and Survival
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      popmat[1, ] <- fert.stoch  # Set fertility in population matrix
      diag(popmat[2:stages, ]) <- Sx.stoch[-stages]  # Set survival in population matrix
      
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Catastrophic Mortality Adjustment (if used)
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      cat.mort.sd <- 0.05  # Standard deviation for catastrophic mortality
      cat.mort.alpha <- estBetaParams(cat.mort, cat.mort.sd^2)$alpha  # Alpha for beta distribution
      cat.mort.beta <- estBetaParams(cat.mort, cat.mort.sd^2)$beta  # Beta for beta distribution
      cat.mort.stoch <- rbeta(1, cat.mort.alpha, cat.mort.beta)  # Sample catastrophic mortality
      
      if (catastrophe.used == 1) {
        catastrophe <- rbinom(1, 1, cat.pr)# Randomly determine if catastrophe occurs
        if (catastrophe == 1) {
          diag(popmat[2:(stages), ]) <- (cat.mort.stoch*Sx.stoch[-stages])} # Apply catastrophe
        if (catastrophe == 0) {
          diag(popmat[2:(stages), ]) <- Sx.stoch[-stages]}
        #popmat[stages,stages] <- 0
      }
      
      # Project population to the next year
      n.mat[,i+1] <- popmat %*% n.mat[,i] # Matrix multiplication to project population
    } # End year loop (i)
    
    # Store population sums for each iteration
    n.sums.mat[e,] <- as.vector(colSums(n.mat))
    nad.sums.mat[e,] <- as.vector(colSums(n.mat[2:stages,]))# Adult females
    njv.sums.mat[e,] <- as.vector(sum(n.mat[1,], na.rm=T))# Juveniles

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Calculate the number of females sterilized (adults and juveniles)
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (ster.which == 2) {# Sterilize both adults and juveniles
      nad.ster.mat[e,] <- nad.sums.mat[e,] * ifelse((ster.stoch) > 1, 1, (ster.stoch)) # adult females
      njv.ster.mat[e,] <- njv.sums.mat[e,] * ifelse((ster.stoch*sex.ratio) > 1, 1, (ster.stoch*sex.ratio)) # juveniles
    }
    if (ster.which == 1) {# Sterilize only adults
      nad.ster.mat[e,] <- nad.sums.mat[e,] * ifelse((ster.stoch) > 1, 1, (ster.stoch)) # young of year excluded
      njv.ster.mat[e,] <- 0# No juvenile sterilization
    }
    
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Calculate Search Time for Koala Capture (Stochastic)
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    D.it <- nad.sums.mat[e,] / area# Population density
    srch.sd <- 0.05 # # Standard deviation for search time arbitrarily set to 5%
    srch.det <- as.numeric(coef(fit.sexpd)[1]) * (D.it - as.numeric(coef(fit.sexpd)[2]))^as.numeric(coef(fit.sexpd)[3]) # Search time based on population density
    srch.stoch <- rnorm(t+1, srch.det, srch.sd*srch.det) # Stochastic search time
    
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Calculate Sterilization Cost
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    cost.mat[e, ] <- (nad.ster.mat[e, ] * (((srch.stoch / 60) * cph) + (hpkc * cph) + spk)) + (njv.ster.mat[e, ] * spk)  # Total cost for adults + juveniles (if any)
  }  # End iteration loop (e)
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Calculate Median and Confidence Intervals for Population Size
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  n.up <- n.lo <- n.mn <- rep(0,(t+1))
  for (q in 1:(t+1)) {
    n.mn[q] <- median(n.sums.mat[,q], na.rm=T) # Median population size
    n.up[q] <- as.vector(quantile(n.sums.mat[,q], probs=0.975, na.rm=T)) ## Upper bound
    n.lo[q] <- as.vector(quantile(n.sums.mat[,q], probs=0.025, na.rm=T)) # Lower bound
  }
  # Plot population over time with confidence intervals
  plot(yrs, n.mn, type="l", xlab="year",ylab="Nf",xlim=c(yr.now,yr.end),ylim=c(min(n.lo),max(n.up)))
  lines(yrs, n.up, lty=2, col="red") # Upper bound
  lines(yrs, n.lo, lty=2, col="red")# Lower bound
  
  # Store the results for this sterilization proportion
  n.mn.mn[s] <- median(n.mn, n.rm=T)# Median population
  n.mn.lo[s] <- quantile(n.mn, probs=0.025, n.rm=T)# Lower bound
  n.mn.up[s] <- quantile(n.mn, probs=0.975, n.rm=T)# Upper bound
  
  ## store minimum n for this sterilisation proportion
  n.min <- p.qext <- rep(0,(t+1))
  for (q in 1:(t+1)) {
    n.min[q] <- min(n.sums.mat[,q], na.rm=T)  # Minimum population size
    p.qext[q] <- length(which(n.sums.mat[,q] < mvp)) / iter# Probability of extinction
  }

  # Calculate proportional change from starting population
  prop.mn[s] <- (n.mn[t+1] * 1/sex.ratio) / (start.pop * 1/sex.ratio) # Proportion of starting population at end
  prop.lo[s] <- (n.lo[t+1] * 1/sex.ratio) / (start.pop * 1/sex.ratio) # Lower bound of proportion
  prop.up[s] <- (n.up[t+1] * 1/sex.ratio) / (start.pop * 1/sex.ratio) # Upper bound of proportion
  
  # calculate minimum N over all stochastic time series
  min.mn[s] <- mean(n.min, na.rm=T)
  min.lo[s] <- quantile(n.min, probs=0.025, na.rm=T)
  min.up[s] <- quantile(n.min, probs=0.975, na.rm=T)
  
  # calculate probability of extinction over all e simuations
  p.qext.max[s] <- max(p.qext)

  # Calculate number of sterilized females (mean and confidence intervals)
  nster.mn[s] <- mean(rowSums(nad.ster.mat, na.rm=T) + rowSums(njv.ster.mat, na.rm=T), na.rm=T)
  nster.lo[s] <- quantile(rowSums(nad.ster.mat, na.rm=T) + rowSums(njv.ster.mat, na.rm=T), probs=0.025, na.rm=T)
  nster.up[s] <- quantile(rowSums(nad.ster.mat, na.rm=T) + rowSums(njv.ster.mat, na.rm=T), probs=0.975, na.rm=T)
  
  ##   # Calculate total sterilization cost + CI
  cost.mn[s] <- mean(rowSums(cost.mat, na.rm=T), na.rm=T) / 10^6 # (in $millions)
  cost.lo[s] <- as.numeric(quantile(rowSums(cost.mat, na.rm=T), probs=0.025, na.rm=T) / 10^6) # (in $millions)
  cost.up[s] <- as.numeric(quantile(rowSums(cost.mat, na.rm=T), probs=0.975, na.rm=T) / 10^6) # (in $millions)
  
  # Calculate cost per year
  cost.yr.mn.mat[s,] <- apply(cost.mat, MARGIN=2, mean, na.rm=T)
  cost.yr.lo.mat[s,] <- apply(cost.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
  cost.yr.up.mat[s,] <- apply(cost.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
  
  # Print progress
  print("################################################")
  print(paste("proportion females sterilised = ", ster.vec[s], sep=""))
  print("################################################")
  
} # End sterilization proportion loop (s)


##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## VISUALISATION & SAVING MODEL OUTPUTS
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#visualizing the results of the stochastic population projection with sterilization across various metrics, 
#including population size, sterilization costs, and the number of females sterilized

##++++++++++++
# 1) Mean Population vs. Sterilization Proportion Plot
#.   Plot: This plot shows how the mean population size changes as the proportion of sterilized females increases.
#.   Red Lines: Represent the 95% confidence interval for the population size.
#.   NF Data Frame: Stores the sterilization proportion and corresponding population statistics.
par(mfrow=c(1,1))
plot(ster.vec, n.mn.mn, type="l", xlab="proportion females sterilised", ylab="mean population (total)", ylim=c(min(n.mn.lo), max(n.mn.up)))
lines(ster.vec, n.mn.lo, lty=2, col="red")
lines(ster.vec, n.mn.up, lty=2, col="red")

NF<-data.frame(ster.vec,n.mn.mn,n.mn.lo,n.mn.up)
colnames(NF)<-c("prop steril","median pop", "CImin pop", "CImax pop")
file_nameNF <- paste("Total_popFemale(Median)_scenario#", ster.which, ".csv", sep="")# Dynamically create the file name based on the value of ster.which
write.table(NF, file=file_nameNF, sep=",", dec=".", row.names=FALSE, col.names=TRUE)# Save the table with the dynamic file name

##++++++++++++
# 2) Minimum Population vs. Sterilization Proportion
#.   Plot: Shows the minimum population size over the projection, depending on sterilization proportion.
#.   Red Lines: Represent the confidence interval for the minimum population.
#.   mean.min.NF Data Frame: Stores the results for the minimum population size.
par(mfrow=c(2,2))
plot(ster.vec, min.mn, type="l", xlab="proportion females sterilised", ylab="mean min Nf (total)", ylim=c(min(min.lo), max(min.up)))
lines(ster.vec, min.lo, lty=2, col="red")
lines(ster.vec, min.up, lty=2, col="red")

mean.min.NF <- data.frame(ster.vec,min.mn,min.lo,min.up)
colnames(mean.min.NF)<-c("prop steril","median min NF", "CImin minNF", "CImax minNF")
file_nameminNF <- paste("Min_initial_popFemale(Median)_scenario#", ster.which, ".csv", sep="")# Dynamically create the file name based on the value of ster.which
write.table(mean.min.NF, file=file_nameminNF, sep=",", dec=".", row.names=FALSE, col.names=TRUE)# Save the table with the dynamic file name

##++++++++++++
# 3) Proportion of Starting Population
#.   Plot: This plot shows the proportion of the starting population remaining by the end of the projection, as a function of the sterilization proportion.
#.   prop.min.NF Data Frame: Stores the proportion results.
plot(ster.vec, prop.mn, type="l", xlab="proportion females sterilised", ylab="mean proportion of start N (f + m)", ylim=c(min(prop.lo), max(prop.up)))
lines(ster.vec, prop.lo, lty=2, col="red")
lines(ster.vec, prop.up, lty=2, col="red")
prop.min.NF<-data.frame(ster.vec,prop.mn,prop.lo,prop.up)
colnames(prop.min.NF)<-c("prop steril","median prop minNF", "CImin prop minNF", "CImax prop minNF")
file_namePropminNF <- paste("PMin_initial_popFemale&Male(Median)_scenario#", ster.which, ".csv", sep="")# Dynamically create the file name based on the value of ster.which
write.table(prop.min.NF, file=file_namePropminNF, sep=",", dec=".", row.names=FALSE, col.names=TRUE)# Save the table with the dynamic file name

##++++++++++++
# 4) Number of Females Sterilized
#.   Plot: This shows the number of females sterilized based on the proportion of sterilization.
#.   nster.min.NF Data Frame: Stores the number of females sterilized and their confidence intervals.
plot(ster.vec, nster.mn, type="l", xlab="proportion females sterilised", ylab="# actually sterilised", ylim=c(min(nster.lo), max(nster.up)))
lines(ster.vec, nster.lo, lty=2, col="red")
lines(ster.vec, nster.up, lty=2, col="red")

nster.min.NF<-data.frame(ster.vec,nster.mn,nster.lo,nster.up)
colnames(nster.min.NF)<-c("prop steril","median nster minNF", "CImin nster minNF", "CImax nster minNF")
file_nameNsterminNF <- paste("Proportion_Actual_Femalsteril(Median)_scenario#", ster.which, ".csv", sep="")# Dynamically create the file name based on the value of ster.which
write.table(nster.min.NF, file=file_nameNsterminNF, sep=",", dec=".", row.names=FALSE, col.names=TRUE)# Save the table with the dynamic file name

##++++++++++++
# 5) Total Sterilization Costs
#.  Plot: This plot displays the total costs of sterilizing females based on the proportion sterilized.
#.  cost.min.NF Data Frame: Stores the sterilization costs and their confidence intervals.
plot(ster.vec, cost.mn, type="l", xlab="proportion females sterilised", ylab="total costs ($m)", ylim=c(min(cost.lo), max(cost.up)))
lines(ster.vec, cost.lo, lty=2, col="red")
lines(ster.vec, cost.up, lty=2, col="red")

cost.min.NF<-data.frame(ster.vec,cost.mn,cost.lo,cost.up)
colnames(cost.min.NF)<-c("prop steril","median cost minNF", "CImin cost minNF", "CImax cost minNF")
#write.table(cost.min.NF, file="Cost_sterilisation(Median)_scenario1.csv",sep=",", dec = ".", row.names = F, col.names = T)
file_namecost.min.NF <- paste("Cost_sterilisation(Median)_scenario#", ster.which, ".csv", sep="")# Dynamically create the file name based on the value of ster.which
write.table(cost.min.NF, file=file_namecost.min.NF, sep=",", dec=".", row.names=FALSE, col.names=TRUE)# Save the table with the dynamic file name


##++++++++++++
# 6) Yearly Sterilization Costs (Contour Plot)
rownames(cost.yr.mn.mat) <- seq(0, 0.9, 0.1)  # Set row names for sterilization proportions
colnames(cost.yr.mn.mat) <- seq(0, t, 1)  # Set column names for years

matplot <- cost.yr.mn.mat  # Median cost matrix
matplotlo <- cost.yr.lo.mat  # Lower bound matrix
matplotup <- cost.yr.up.mat  # Upper bound matrix

# Contour plot of costs
cpycont2 <- plot_ly(z = ~matplot, autocontour = TRUE, type = "contour", line = list(smoothing = 0.85), 
                    contours = list(showlabels = TRUE, labelfont = list(size = 14, color = "white"))) %>%
  colorbar(title = "cost ($)") %>%
  layout(xaxis = list(title = "year"), yaxis = list(title = "proportion sterilised"))
cpycont2

file_nameYrCost <- paste("MedianSterilisationCost_Years_(Median)_scenario#", ster.which, ".csv", sep="")# Dynamically create the file name based on the value of ster.which
write.table(matplot, file=file_nameYrCost, sep=",", dec=".", row.names=FALSE, col.names=TRUE)# Save the table with the dynamic file name

