#_________________ Applied Statistics 2024/2025 _________________________________

#### 7.3.Markov models [markovchain package] ####
#_______________________________________________#

library(markovchain)
library(igraph)


### Summary:
### - EXAMPLE 1: Weather example
### - EXAMPLE 2: Irreducible Markov Chain
### - EXAMPLE 3: Gambler Ruin example
### - ...

# Source: https://cran.r-project.org/web/packages/markovchain/vignettes/an_introduction_to_markovchain_package.pdf

#_______________________________________________________________________________
##### EXAMPLE 1: Weather example #####

# Let use the following transition matrix 
#          sunny  cloudy   rain
#   sunny  |0.70    0.2     0.1 |
#   cloudy |0.3     0.4     0.3 |
#   rain   |0.2     0.45    0.35|

# Create a markovchain object 
mcWeather <- new("markovchain", 
                 # states: a character vector, listing the states for which transition probabilities are defined.
                 states = c("sunny", "cloudy", "rain"),
                 # transitionMatrix: the probabilities of the transition matrix.
                 transitionMatrix = matrix(data = c(0.70, 0.2, 0.1, 
                                                    0.3, 0.4, 0.3, 
                                                    0.2, 0.45, 0.35), byrow = TRUE, nrow = 3), 
                 # name: optional character element to name the DTMC (Discrete Time Markov Chain)
                 name = "Weather")  

# The examples that follow shows how operations on markovchain objects can be easily performed.
# For example, using the previously defined matrix we can find what is the probability
# distribution of expected weather states in two and seven days, given the actual state to be
# cloudy.

initialState <- c(0, 1, 0)
after2Days <- initialState * (mcWeather * mcWeather)
after2Days

after7Days <- initialState * (mcWeather ^ 7)

round(after7Days, 3)

# The transition matrix of a markovchain object can be displayed using print or show methods
# Similarly, the underlying transition probability diagram can be plotted 
# by the use of plot method which is based on igraph package.

print(mcWeather)
x11()
plot(mcWeather)


# The conditional distribution of weather states, given that current day’s weather is sunny,
# is given by following code.

conditionalDistribution(mcWeather, "sunny")

steadyStates(mcWeather)


## Classification of states 

# Absorbing states are determined by means of absorbingStates method.
absorbingStates(mcWeather)  # no absorbing states
 
summary(mcWeather)  ##classification of states

# The function 'is.accessible' permits to investigate whether a state sj is accessible from state
# si, that is whether the probability to eventually reach sj starting from si is greater than zero.

is.accessible(object = mcWeather, from = "sunny", to = "rain")

# If a DTMC is irreducible, all its states share the same periodicity. 
# Then, the period function returns the periodicity of the DTMC, provided that it is irreducible.
# The example that follows shows how to find if a DTMC is reducible or irreducible by means 
# of the function is.irreducible and, in the latter case, the method
# period is used to compute the periodicity of the chain.

is.irreducible(mcWeather)

period(mcWeather)

#_______________________________________________________________________________
##### EXAMPLE 2: Irreducible Markov Chain #####

IrrMatr <- matrix(0,nrow = 5, ncol = 5)
IrrMatr[1,] <- c(0, 1/3, 0, 2/3, 0)
IrrMatr[2,] <- c(1/2, 0, 0, 0, 1/2)
IrrMatr[3,] <- c(0, 0, 1/2, 1/2, 0)
IrrMatr[4,] <- c(0, 0, 1/2, 1/2, 0)
IrrMatr[5,] <- c(0, 0, 0, 0, 1)
statesNames <- letters[1:5]
IrrMatr
IrrMc <- new("markovchain", 
             transitionMatrix = IrrMatr, 
             name = "Irreducible MC", 
             states = statesNames)

summary(IrrMc)

x11()
plot(IrrMc)

#_______________________________________________________________________________
##### EXAMPLE 3: Gambler Ruin example #####

# It is possible for a Markov chain to have more than one stationary distribution, 
# as the gambler ruin example shows.
gamblerRuinMarkovChain <- function(moneyMax, prob = 0.5) {
                              m <- matrix(0,nrow = moneyMax + 1, ncol = moneyMax + 1)
                              m[1,1] <- m[moneyMax + 1,moneyMax + 1] <- 1
                              states <- as.character(0:moneyMax)
                              rownames(m) <- colnames(m) <- states
                              for(i in 2:moneyMax){
                                    m[i,i-1] <- 1 - prob
                                    m[i, i + 1] <- prob
                              } 
                              new("markovchain", 
                                  transitionMatrix = m, 
                                  name = paste("Gambler ruin", moneyMax, "dim", sep = " ")
                                  )
                          }

mcGR4 <- gamblerRuinMarkovChain(moneyMax = 4, prob = 0.5)
steadyStates(mcGR4)
plot(mcGR4)
summary(mcGR4)


#_______________________________________________________________________________
# 1) Simulation of a DTMC (Discrete Time Markov Chain)
# Simulating a random sequence from an underlying DTMC is quite easy thanks to the function
# 'rmarkovchain'. The following code generates a year of weather states according to mcWeather
# underlying stochastic process.

weathersOfDays <- rmarkovchain(n = 365, object = mcWeather, t0 = "sunny")
weathersOfDays[1:30]

# 2) Fitting a DTMC 
# A time homogeneous Markov chain can be fit from given data. Four methods have been
# implemented within current version of markovchain package: maximum likelihood, maximum
# likelihood with Laplace smoothing, Bootstrap approach, maximum a posteriori.

weatherFittedMLE <- markovchainFit(data = weathersOfDays, method = "mle", name = "Weather")
weatherFittedMLE$estimate
weatherFittedMLE$standardError

# 3) Predict from a DTMC
# The n-step forward predictions can be obtained using the predict methods explicitly written
# for markovchain and markovchainList objects.
# The prediction is the mode of the conditional distribution of Xt+1 given Xt = sj, 
# being sj the last realization of the DTMC.

# The 3-days forward predictions from markovchain object can be generated as follows, 
# assuming that the last two days were respectively “cloudy” and “sunny”.

predict(object = weatherFittedMLE$estimate, newdata = c("cloudy", "sunny"), n.ahead = 3)

# Statistical tests: assessing the Markov property (verifyMarkovProperty) 
# and the stationary (assessStationarity) of a Markov chain sequence

# The verifyMarkovProperty function verifies whether the Markov property holds for the given chain.
verifyMarkovProperty(weathersOfDays)


# Another artificial example 
sample_sequence <- c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", "b", "b", "b", "a")
verifyMarkovProperty(sample_sequence)


# ... Other Applications (pag 56)
# from https://cran.r-project.org/web/packages/markovchain/vignettes/an_introduction_to_markovchain_package.pdf

