# Apply a Wiener diffusion model to test for SAT variations in a choice RT task

rm( list=ls() )

# Data are modelled separately for each subject
# the following script demonstrates the analysis procedure for data from a single subject
# participant data should be saved in individual files and the procedure can be repeated

# DATA SIMULATION

# Assume that the data in this script is from a color discrimination task  
# The task is chosen since it is a typical 2CRT paradigm
# In the task, stimuli are classified according to their color
# For the simulation, assume that the stimuli are recoded into 2 responses: correct, incorrect

# To investigate a speed-accuracy trade-off, assume that there are 2 conditions
# one where subjects are instructed to emphasize speed, and one where accuracy is emphasized
# An extension to the paradigm could be having a range of stimuli difficulty
# as it would increase robustness of the estimation (Voss et al., 2013)

# Model assumptions
# The task has previously been validated with the model 
# (Voss, Rothermund, & Voss, 2004; Spaniol et al, 2011)
# and similar perceptual discrimination tasks have been successfully used
# with the diffusion model (Ratcliff and Rouder, 1998; Ratcliff & McKoon, 2008)
# If a new paradigm was used, the assumptions and validity of the model would have to be tested first

# Responses are classified on a binary scale (correct, incorrect) 
# so it meets the model assumptions
# Note that an alternative to the 2 response boundaries would be the actual responses
# for example, 'choice 1' and 'choice 2' 
# If actual responses are used, then a priori bias toward each boundary
# and separate drift rates for each boundary have to be used
# It is useful to use actual boundaries for data with low number of incorrect trials

# Predict that boundary separation varies across conditions 
# as a result of the experimental manipulation
# The variation will therefore produce changes in the the speed-accuracy tradeoff
# It is expected that boundary separation will be higher for the accuracy emphasis condition
# than the speed emphasis condition

# Decide on which parameters are allowed to vary
# Assumptions
# 1. Constant drift rate: stimulus difficulty is the same for both response options 
# Most two choice task paradigms meet this assumption (Ratcliff & McKoon, 2008) 
# If the task uses varying stimulus difficulty, then this assumption is not met 
# ex. the flanker task, where congruent stimuli are easier than incongruent (White et al., 2011)
# 2. Constant starting point: no initial bias towards either of the two boundaries

# Choose number of trials
# Large trial numbers help with estimation
# but they also produce learning effects (Dutilh, Kryptos, & Wagenmakers, 2011)
# A smaller number of trials is acceptable in this situation 
# since the assumed paradigm has been previously tested
# and it avoids biasing effects of practice on cognitive processes
ntrials <- 200

# Collect data
# Here, random data are generated 
# to model a situation with a notable shift in the speed-accuracy tradeoff

# The script uses a package by Wabersich and Vanderckhove (2014) 
# which provides options for calculating diffusion functions in R software
install.packages ("RWiener")    # install and load package for diffusion functions
library("RWiener")              # load package

# Use rwiener function to create samples of data
# rwiener uses the rejection-based algorithm outlined in Tuerlinckx et al (2001)
speed.condition <- rwiener (ntrials, alpha=1, tau=.3, beta=0.5, delta=0.5)
accuracy.condition <- rwiener (ntrials, alpha=2, tau=.3, beta=0.5, delta=0.5)
# The chosen parameter values are based on those typically obtained 
# when fitting the model to real data (Ratcliff & Tuerlinckx, 2002; Ratcliff & McKoon, 2008)
# drift rate: delta = 0.5
# starting point: beta = 0.5, since there shouldn't be any bias towards the correct or incorrect response 
# boundary separation: alpha = 2 and 1
# alpha chosen to be larger in the condition with accuracy emphasis (boundaries further apart) than speed
# non-decision time: tau = 0.3, all RTs > 0.3

# View start of data
# Samples contain two columns, one for RTs (q) and one for responses (resp)
head (speed.condition)
head (accuracy.condition)

# Use filter function to remove RT outliers
# Responses removed if RTs > 2SD away from the mean
# Cutoffs are good for RT outliers in some cases but not others (Whelan, 2008)
# Refer Ratcliff and Tuerlinckx (2002) for a more detailed iterative procedure
filter(log(speed.condition$q) > mean(log(speed.condition$q)) - 2 * sd(log(speed.condition$q)),
       log(speed.condition$q) < mean(log(speed.condition$q)) + 2 * sd(log(speed.condition$q))) 
# Note that filter function may not work for RT data. Instead, it may be easier to model function as rt <- rt[ rt < [ put your criterion here ] ]
# For more detailed RT trimming, the package trimr by James Grange offers multiple options for researchers

# Plot simulated data using the wiener_plot function
wiener_plot(speed.condition)
legend( "bottomright",  legend='speed condition' )
wiener_plot(accuracy.condition)
legend( "bottomright",  legend='accuracy condition' )
# The plots demonstrate that responses are fasted in the speed condition, but with more errors
# Note that the scales are different for the graphs as the depends on the range of responses

# Calculate performance measures for each condition
# create a function to compute speed (MRT, VRT) and accuracy (Pc)
performance <- function (dat) {
Pc  <- (sum(dat$resp=='upper'))/length(dat$resp) # proportion correct
VRT <- var(dat$q)                                # variance of reaction time
MRT <- mean(dat$q)                               # mean reaction time
return (c(Pc, VRT, MRT))
}
speed.performance <- performance (speed.condition)
accuracy.performance <- performance (accuracy.condition)

# Display obtained performance for each condition
measures <- as.data.frame (rbind (speed.performance, accuracy.performance))
performance.measures <- round  (measures, digits=2)      # round parameter estimates to 3 digits
names (performance.measures) <- c('percent correct', 'RT variance', 'RT mean')
print (performance.measures)
# Values shown confirm much faster mean speeds in the speed condition 
# and a slightly higher percentage of correct responses in the accuracy condition
# As expected with speed-accuracy tradeoffs, 
# small increases in accuracy cause large decreases in speed

# Check the similarity in accuracy values for each condition
# to see the proportion of correct to incorrect trials (since correct and incorrect are boundaries)
# If values are extremely different then consider using actual responses as boundaries
# so that there are enough trials for the distributions at each threshold 

# STEP 2: Fit model to data

# Use the EZ-diffusion computation to get starting values
# The script contains a function 'get.vaTer'
# Which uses data inputs (Pc, VRT, MRT) to estimate model parameters 
source('EZdiffusion.R')                       # loads EZdiffusion file
# EZ.diffusion script obtained from Wagenmakers et al (2007)
# Copy-pasted code for EZ function in Appendix C
# In case the file doesn't open with command, open it manually and run the script

# EZ estimates are used as initial values for maximum likelihood optimization 
# in other existing diffusion model software (Vanderchkove et al., 2007; Voss & Voss, 2008)

# Use the EZ-diffusion script to calculate parameters for each condition
EZpars.speed <- get.vaTer (speed.performance[1], speed.performance[2], speed.performance[3] )
EZpars.accuracy <- get.vaTer (accuracy.performance[1], accuracy.performance[2], accuracy.performance[3] )

# Convert obtained EZ parameter estimates to values for diffusion model
# Create a function to convert EZ estimates into a vector 
# that can be used a starting values for the diffusion model
# Note that EZ assumes that beta=alpha/2
EZ.estimates <- function (EZpars){
delta <- EZpars[[1]]                      # drift rate ('v' in EZ, 'delta' in DDM)
alpha <- EZpars[[2]]                      # boundary separation ('alpha', i.e. 'a' in EZ and DDM)
beta <- (EZpars[[2]])/2                   # starting point (alpha/2 in EZ, beta in DDM)
tau <- EZpars[[3]]                        # non-decision time ('Ter' in EZ, 'tau' in DDM)
EZ.values <- c(alpha, tau, beta, delta)
}

# starting parameter values for each condition
speed.init <- EZ.estimates (EZpars.speed)
accuracy.init <- EZ.estimates (EZpars.accuracy)

# Maximum likelihood is used here as the optimization criterion 
# Suitable since it is efficient with small datasets 
# after outliers have been removed (Ratcliff & Tuerlinckx, 2001)
# Use in-built likelihood function from Rwiener package with optim () to estimate parameters

# Define a two step optimization process 
# based on recommendations by Vanderchkove (2008) to avoid local minima
# First step uses default Nelder-Mead method
# Second step uses a quasi-Newton method and approximates the Hessian matrix
# optimization method guided by the script example in Appendix C (Wabersich & Vanderchkove, 2014 )
twostep.optim <- function (init.estimates, data) 
  {
step1 <- optim (init.estimates, wiener_deviance, dat=data, hessian=FALSE)
step2 <- optim (step1$par, wiener_deviance, dat=data, method="BFGS", hessian=TRUE)
}

# Apply optimization function to each condition
# In most cases, using EZ estimates as starting values works
# but when it doesn't, the try command handles error-recovery 
try(speed <- twostep.optim (init.estimates=speed.init, data=speed.condition), silent=TRUE)                                               
try(accuracy <- twostep.optim (init.estimates=accuracy.init, data=accuracy.condition), silent=TRUE)

# In the few cases where the EZ initial estimates don't work, use the estimates below 
# Reliability of the estimate below is tested with 10000 simulations
speed <- twostep.optim (init.estimates=c(1,.1,.1,1), data=speed.condition)
accuracy <- twostep.optim (init.estimates=c(1,.1,.1,1), data=accuracy.condition)

# Print obtained parameters in a list
predicted.values <- as.data.frame (rbind (speed$par, accuracy$par))
parameters <- round  (predicted.values[1:4], digits=2)      # round parameter estimates 
names (parameters) <- c('Boundary separation', 'Non-decision time', 'Starting point', 'Drift')
parameters$condition <- c ('speed', 'accuracy')
print (parameters)

# Check observed parameters to see if they fit initial predictions
# Below is a table showing a sample of results that was obtained for each condition

#  Boundary separation Non-decision time Starting point Drift Condition
#                0.95              0.31           0.49  0.58     speed
#                2.04              0.30           0.50  0.56  accuracy

# As the table demonstrates, boundary separation varied by almost 1
# while the other parameters were mostly constant
# In some repetitions of the simulation, drift rate also varies up to a difference of 0.5 

# Once fit has been assessed for parameter estimates, 
# they can then be used as dependent variables for statistics analyses 
# For analyses, standard errors and confidence intervals on parameter estimates 
# should be calculated to assess uncertainty in the estimates

# Use wiener plot function to view fitted model
fitted.speed <- rwiener (ntrials, speed$par [1], speed$par [2], speed$par [3], speed$par [4])
wiener_plot(fitted.speed)
legend( "bottomright",  legend='speed condition' )

fitted.accuracy <- rwiener (ntrials, accuracy$par [1], 
    accuracy$par [2], accuracy$par [3], accuracy$par [4])
wiener_plot(fitted.accuracy)
legend( "bottomright",  legend='accuracy condition' )

# STEP 3: MODEL TESTING AND SELECTION

# The script assumes that the actual model parameters were unknown
# as would be the case in a real experiment
# so model selection will be done by comparing the fit of competing models for the data

# Model selection test 1: Likelihood ratio test
# Null model (H0): no change in boundary separation (alpha=0)
# Alternative model (H1): boundary separation varies across conditions (alpha=1)
# Alternative value (alpha=1) based on obatined parameters in which alpha varied by 1
# Both hypotheses assume that all other parameters have no change

# Calculate likelihood for null and alternate models
# H0 assumes that alpha is the same across conditions, and H1 that alpha varies by 1
# Values based on parameters obtained from fitted model, in which there was a difference of 1
H0.likelihood <- wiener_likelihood(speed$par, speed.condition)+wiener_likelihood(
  c(speed$par[1]+0, speed$par[2:4]), accuracy.condition)
H1.likelihood <- wiener_likelihood(speed$par, speed.condition)+wiener_likelihood(
  c(speed$par[1]+1, speed$par[2:4]), accuracy.condition)

# Examine goodness of fit for each model using the likelihood ratio test
likelihood.ratio.test <- -2*H0.likelihood + 2*H1.likelihood

# Use chi-square distribution to assess goodness of fit
# df calculated using (K*(N-1)) - M (based on White et al, 2010)
# where K is the number of conditions (2), M is the number of free parameters (1)
# and N is the number of bins for free parameters (typically 12)
df = (2*(12-1)-1)
pchisq (likelihood.ratio.test, df, lower.tail=FALSE)          # p-value of obtained likelihood ratio
likelihood.ratio.test > qchisq(0.95, df=19)                      # test obtained value at alpha of 0.05
# provides a value of TRUE, so reject the H0

# Note that the likelihood ratio test is very conservative 
# since it uses predictions obtained from fitting observed distributions
# It is also based on the number of trials and conditions
# To avoid this bias, Monte-Carlo simulations can be used 
# to get a distribution of the fit (Voss et al., 2013) 

# Choose number of simulations
# Note: Around 10000 is ideal, but based on the slow processing time 
# of >30 minutes minutes for 1000 trials, a smaller value of 100 is shown here
nsim <- 100                                  
likelihood.ratio.test.sim <- rep( NaN, nsim )       # empty vector to place results of simulations

# Repeat the likelihood ratio test to obtain a distribution of the test results
for (i in 1:nsim) {
  
  # generate new datasets using the estimated parameters
  speed.condition.sim <- rwiener (ntrials, alpha=speed$par[1], tau=speed$par[2], beta=speed$par[3], delta=speed$par[4])
  accuracy.condition.sim <- rwiener (ntrials, alpha=accuracy$par[1], tau=accuracy$par[2], beta=accuracy$par[3], delta=accuracy$par[4])
  
  # filter outliers
  filter(log(speed.condition.sim$q) < mean(log(speed.condition.sim$q)) + 2 * sd(log(speed.condition.sim$q)),
         log(speed.condition.sim$q) > mean(log(speed.condition.sim$q)) - 2 * sd(log(speed.condition.sim$q))) 
  
  # obtain performance measures for each condition
  speed.performance.sim <- performance (speed.condition.sim)
  accuracy.performance.sim <- performance (accuracy.condition.sim)
  
  # use EZ to convert DVs to model parameters
  EZpars.speed.sim <- get.vaTer (speed.performance.sim[1], speed.performance.sim[2], speed.performance.sim[3] )
  EZpars.accuracy.sim <- get.vaTer (accuracy.performance.sim[1], accuracy.performance.sim[2], accuracy.performance.sim[3] )
  
  # get initial estimates for fitting
  speed.estimates.sim <- EZ.estimates (EZpars.speed.sim)
  accuracy.estimates.sim <- EZ.estimates (EZpars.accuracy.sim)
  
  # fit the DDM 
  # try uses the EZ values as initial estimates, but deals with any errors
  # that arise when the estimates are not suitable as starting values
  try (speed.sim <- twostep.optim (init.estimates=speed.estimates.sim, data=speed.condition.sim),
       silent=TRUE)                
  try (accuracy.sim <- twostep.optim (init.estimates=accuracy.estimates.sim, data=accuracy.condition.sim),
       silent=TRUE) 
  
  # calculate likelihood ratio for simulated data distribution
  try (H0.likelihood <- wiener_likelihood(speed.sim$par, speed.condition.sim)+wiener_likelihood(c
      (speed.sim$par[1]+0, speed.sim$par[2:4]), accuracy.condition.sim), silent=TRUE)
  try (H1.likelihood <- wiener_likelihood(speed.sim$par, speed.condition.sim)+wiener_likelihood(c
      (speed.sim$par[1]+1, speed.sim$par[2:4]), accuracy.condition.sim), silent=TRUE)
  
  likelihood.ratio.test.sim [i] <- -2*H0.likelihood + 2*H1.likelihood
}

# Use distribution of simulated likelihood ratios to assess fit 
# remove values that are NaN (when the try function suppressed errors)
likelihood.ratio.dist <- na.omit(likelihood.ratio.test.sim)    
mean (likelihood.ratio.dist) > qchisq(0.95, df)
pchisq (mean (likelihood.ratio.dist) , df, lower.tail=FALSE)

# show distribution of deviance
# This step was based on the method taught in class
hist( likelihood.ratio.dist)
# show the original likelihood ratio test result we measured on the simulated distribution
abline( v=likelihood.ratio.test, col='blue' )
# Compared the obtained likelihood to the distribution of model likelihoods
pLRT <- mean( likelihood.ratio.dist > likelihood.ratio.test )
print(pLRT)
# Provides a value of TRUE, which suggests that the obtained likelihood ratio falls within the range of possible ratios

# Model selection test 2
# The previous test assumes that only alpha varies and other parameters were fixed
# However, we can also run a test to see if drift also varies across conditions
# The fitted parameters indicate that the information accumulation (drift) decreases in the accuracy condition
# Thus see if a better fit is obtained by a model in which drift is also allowed to vary

# Although it is parsimonious not to vary too many parameters,
# it is a good idea to also test a model in which other parameters are allowed to vary
# as varying only one parameter may bias results by mapping any effect onto that parameter

# Try also varying drift
# in addition to the parameter of interest (boundary separation for the SAT),
# Null model (H0): no change in any parameters 
# Alternative model (H1): boundary separation and drift vary across conditions (alpha=1, drift=0.5)
# Alternative values based on obtained parameters 
# in which alpha varies by around 1 and drift by around 0.5
# Both hypotheses assume that all other parameters have no change

H0.likelihood.2 <- wiener_likelihood(speed$par, speed.condition)+wiener_likelihood(c(speed$par), accuracy.condition)
H1.likelihood.2 <- wiener_likelihood(speed$par, speed.condition)+wiener_likelihood(c(speed$par[1]+1, speed$par[2:3], speed$par[4]-0.5), accuracy.condition)
likelihood.ratio.test.2 <- -2*H0.likelihood.2 + 2*H1.likelihood.2

# Test probability of obtaining fitted parameters given the original parameters
df.2 = (2*(12-1)-2)
pchisq (likelihood.ratio.test.2, df.2, lower.tail=FALSE)
likelihood.ratio.test.2 > qchisq(0.95, df.2) 
# Gives a value of TRUE, suggesting that a better fit might be obtained by also allowing drift to vary

# Repeat using Monte-Carlo simulations to see if this is the case for multiple datasets

likelihood.ratio.test.sim.2 <- rep( NaN, nsim )       # empty vector to place results of simulations

# Repeat the likelihood ratio test to obtain a distribution of the test results
for (i in 1:nsim) {
  
  # generate new datasets using the estimated parameters
  speed.condition.sim <- rwiener (ntrials, alpha=speed$par[1], tau=speed$par[2], beta=speed$par[3], delta=speed$par[4])
  accuracy.condition.sim <- rwiener (ntrials, alpha=accuracy$par[1], tau=accuracy$par[2], beta=accuracy$par[3], delta=accuracy$par[4])
  
  # filter outliers
  filter(log(speed.condition.sim$q) < mean(log(speed.condition.sim$q)) + 2 * sd(log(speed.condition.sim$q)),
         log(speed.condition.sim$q) > mean(log(speed.condition.sim$q)) - 2 * sd(log(speed.condition.sim$q))) 
  
  # obtain performance measures for each condition
  speed.performance.sim <- performance (speed.condition.sim)
  accuracy.performance.sim <- performance (accuracy.condition.sim)
  
  # use EZ to convert DVs to model parameters
  EZpars.speed.sim <- get.vaTer (speed.performance.sim[1], speed.performance.sim[2], speed.performance.sim[3] )
  EZpars.accuracy.sim <- get.vaTer (accuracy.performance.sim[1], accuracy.performance.sim[2], accuracy.performance.sim[3] )
  
  # get initial estimates for fitting
  speed.estimates.sim <- EZ.estimates (EZpars.speed.sim)
  accuracy.estimates.sim <- EZ.estimates (EZpars.accuracy.sim)
  
  # fit the DDM 
  # try uses the EZ values as initial estimates, but deals with any errors
  # that arise when the estimates are not suitable as starting values
  try (speed.sim <- twostep.optim (init.estimates=speed.estimates.sim, data=speed.condition.sim),
       silent=TRUE)                
  try (accuracy.sim <- twostep.optim (init.estimates=accuracy.estimates.sim, data=accuracy.condition.sim),
       silent=TRUE) 
  
  # calculate likelihood ratio for simulated data distribution
  try (H0.likelihood <- wiener_likelihood(speed.sim$par, speed.condition.sim)+wiener_likelihood(c
    (speed.sim$par), accuracy.condition.sim), silent=TRUE)
  try (H1.likelihood <- wiener_likelihood(speed.sim$par, speed.condition.sim)+wiener_likelihood(c
    (speed$par[1]+1, speed$par[2:3], speed$par[4]-0.5), accuracy.condition.sim), silent=TRUE)
  likelihood.ratio.test.sim.2 [i] <- -2*H0.likelihood + 2*H1.likelihood  
}

likelihood.ratio.dist.2 <- na.omit(likelihood.ratio.test.sim.2)    
hist( likelihood.ratio.dist.2)
# show the original likelihood ratio test result we measured on the simulated distribution
abline( v=likelihood.ratio.test.2, col='blue' )
# Compared the obtained likelihood to the distribution of model likelihoods
pLRT.2 <- mean( likelihood.ratio.dist.2 > likelihood.ratio.test.2 )
print(pLRT.2)

# This simulation shows a similar result to the previous model with only one parameter
# which shows that varying both 1 and 2 parameters offer a good fit to the data

# Model selection: avoiding complexity
# Both the models in tests 1 and 2 provided a good fit to the data
# However, the second model varied an extra parameter
# To avoid over-fitting the model (Pitt & Myung, 2002) the AIC criterion is used to test models
AIC <- - 2*log(MLk) + 2*(pk)
# Formula obtained from Myung (2000)
# where log(MLk) is the log of the maximum likelihood (ML) function for model k
# and pk is the number of free parameters in the model
# The AIC calculations use the log likelihoods calculated above in the likelihood ratio test

# Model 1: varying no parameters
AIC.0k <- - 2*(H0.likelihood) + 2*(0)
# Model 2: varying only alpha
AIC.1k <- - 2*(H1.likelihood) + 2*(1)
# Model 3: varying alpha and delta
AIC.2k <- - 2*(H0.likelihood.2) + 2*(2)
print (c(AIC.0k, AIC.1k, AIC.2k))

# The ideal model is the one with the smallest AIC value
# In this case it is the second model (which varies the alpha parameter, but not delta)
# An example output of the AIC values is pasted below
# 1277.1043  630.9008 1281.1043
# The output shows that varying one parameter (AIC=630.9008)
# is much better than varying none (1277.10) or two (1281.10)
# This makes sense since the data generated for this experiment varied only 1 parameter (alpha)

# Note there is also a function in the RWiener package for calculating AIC
# But I wanted to try calculating it myself so I didn't use their option

###
# Model assessment
# Since this is a simulation, can test if model differs from values used to simulate parameters
# see how close fitted model is to the actual parameters
# only works in this case since original parameters are known

original.speed.parameters <- c(1, .3, 0.5, 0.5)
original.accuracy.parameters <- c(2, .3, 0.5, 0.5)

# H0-model: original parameters
# H1-model: estimated parameters after fitting model
H0.likelihood.3 <- wiener_likelihood(original.speed.parameters, speed.condition)+wiener_likelihood(original.accuracy.parameters , accuracy.condition)
H1.likelihood.3  <- wiener_likelihood(speed$par, speed.condition)+wiener_likelihood(accuracy$par, accuracy.condition)
likelihood.ratio.test.3 <- -2*H0.likelihood.3 + 2*H1.likelihood.3

# Test probability of obtaining fitted parameters given the original parameters
df = (2*(12-1)-3)
pchisq (likelihood.ratio.test.3, df, lower.tail=FALSE)
likelihood.ratio.test.3 > qchisq(0.95, df) 
# Gives an output of false, so fail to reject H0
# Which shows that the fitted model does not differ significantly from the model used to fit it
# Suggest that the fitting procedure is able to recover original parameters 

# Can repeat the simulation shown above to test ability to recover parameters in multiple cases

# Another option is to use the Rwiener in-built functions  
# Compare the deviance from the dataset to the original parameters used to generate the data
wiener_deviance (original.speed.parameters, speed.condition)
wiener_deviance(speed$par, speed.condition)
# deviance value for obtained and estimated parameters is similar

# plot obtained data and data that would be predicted based on the obtained parameters
# first plot correct responses
# speed condition - points for obtained data and a line for predicted data
curve(dwiener(x, speed$par [1], speed$par [2], speed$par [3], speed$par [4], rep("upper", length(x))),
      main="Density of RT responses", type="p", ylab="density", xlab="RT quantile", 
     xlim=c(0,3), col="green", add=FALSE)
curve(dwiener(x, alpha=1, tau=.3, beta=0.5, delta=0.5, rep("upper", length(x))),
      col="darkgreen", add=TRUE)
# accuracy condition
curve(dwiener(x, accuracy$par [1], accuracy$par  [2], accuracy$par  [3], accuracy$par  [4], rep("upper", length(x))),
    type="p", col="grey", add=TRUE)
curve(dwiener(x, alpha=2, tau=.3, beta=0.5, delta=0.5, rep("upper", length(x))),
      col="black", add=TRUE)
# repeat for error responses
# speed condition
curve(dwiener(x, speed$par [1], speed$par [2], speed$par [3], speed$par [4], rep("lower", length(x))),
      type="p",col="pink", add=TRUE)
curve(dwiener(x, alpha=1, tau=.3, beta=0.5, delta=0.5,  rep("lower", length(x))),
      col="red", xlab="RT quantile", add=TRUE)
# accuracy condition
curve(dwiener(x, accuracy$par [1], accuracy$par  [2], accuracy$par  [3], accuracy$par  [4], rep("lower", length(x))),
      type="p", col="blue", add=TRUE)
curve(dwiener(x, alpha=2, tau=.3, beta=0.5, delta=0.5, rep("lower", length(x))),
      col="darkblue", add=TRUE)

legend( "topright", legend=c('Speed (correct)', 'Accuracy (correct)', 
    'Speed (incorrect)','Accuracy (incorrect)'), title='Condition (response type)',
    col=c('green','black','red', 'darkblue'), pch=c(1, 1, 1, 1), lty=c(1, 1, 1, 1))
