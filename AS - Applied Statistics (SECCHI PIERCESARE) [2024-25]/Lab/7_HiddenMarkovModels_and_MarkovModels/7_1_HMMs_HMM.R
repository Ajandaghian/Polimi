#_________________ Applied Statistics 2024/2025 _________________________________

#### 7.1.Hidden Markov models [HMM package]  ####
#_______________________________________________#

library(HMM)

### Summary:
### - EXAMPLE 0: A simulated example for understanding how to use the HMM package
### - EXAMPLE 1: Dishonest casino example
### - EXAMPLE 2: DNA sequence


#_______________________________________________________________________________
##### EXAMPLE 0: A simulated example for understanding how to use the HMM package #####

# How to initialize our HMM
?initHMM

# Initialise HMM nr.1
hmm1 = initHMM(c("X","Y"),               # States: Vector with the names of the states
               c("a","b","c")            # Symbols: Vector with the names of the symbols.
               )
hmm1
# Initialise HMM nr.2
hmm2 = initHMM(c("X","Y"),               # States: Vector with the names of the states
               c("a","b"),               # Symbols: Vector with the names of the symbols.
               c(.3,.7),                 # startProbs: Vector with the starting probabilities of the states.
               matrix(c(.9,.1,.1,.9),2), # transProbs: Stochastic matrix containing the transition probabilities between the states.
               matrix(c(.3,.7,.7,.3),2)  # emissionProbs: Stochastic matrix containing the emission probabilities of the states.
               )
hmm2

### A simulated example

# Initial HMM
hmm3 = initHMM(c("A","B"),
               c("L","R"),
               transProbs = matrix(c(.9,.1,.1,.9),2),
               emissionProbs = matrix(c(.5,.51,.5,.49),2)
               )
print(hmm3)

# Simulate from the HMM
# Simulates a path of states and observations for a given Hidden Markov Model.
simHMM(hmm3, 100)

# Let's define a sequence of observation
a = sample(c(rep("L",100), rep("R",300)))
b = sample(c(rep("L",300), rep("R",100)))
observation = c(a,b)
observation # vector of symbols

# Viterbi
?viterbi
viterbi(hmm3, observation) # most probable path of hidden states

# Baum-Welch
?baumWelch
bw = baumWelch(hmm3, observation, maxIterations = 100) 
print(bw$hmm) # optimal parameters to the HMM

# Posterior probabilities for the states
# This function computes the posterior probabilities of being in state X at time k 
# for a given sequence of observations and a given Hidden Markov Model.
?posterior
posterior(hmm3, observation) 


#_______________________________________________________________________________
##### EXAMPLE 1: Dishonest casino example (from Durbin et. al. 1999) #####
# A dishonest Casino uses a fair dice most of the time, but switches to the loaded dice once in a while.
# The probabilities of the fair die are (1/6, ..., 1/6) for throwing ("1",...,"6").
# The probabilities of the loaded die are (1/10, ..., 1/10, 1/2) for throwing ("1",...,"5","6"). 
# The observer doesn’t know which die is actually taken (the state is hidden), 
# but the sequence of throws (observations) can be used to infer which die (state) was used.
# Can we detect which die is in use at any given time, just by observing the sequence of rolls?
# States = {Fair, Loaded}
# P(Fair --> Fair) = 0.95
# P(Fair --> Loaded) = 0.05
# P(Loaded --> Fair) = 0.1
# P(Loaded --> Loaded) = 0.9

States = c('Fair', 'Loaded')
Symbols = c(1, 2, 3, 4, 5, 6)
observations = c(3, 2, 2, 1, 2, 3, 6, 6, 6, 6)
transProbs = matrix( c(0.95, 0.1, 
                       0.05,  0.9), 2)
transProbs
emissProbs = matrix( c(1/6, 1/10, 
                       1/6, 1/10, 
                       1/6, 1/10,
                       1/6, 1/10,
                       1/6, 1/10,
                       1/6, 1/2 ), 2)
emissProbs
hmm = initHMM(States, Symbols, c(0.5, 0.5), transProbs, emissProbs)
hmm
viterbi(hmm, observations) # we get the most probable path of hidden states


# Dishonest casino example
dishonestCasino()


#_______________________________________________________________________________
##### EXAMPLE 2: DNA sequence #####
# In a HMM, the nucleotides (A, T, G, C) found at a particular position 
# in a sequence depends on the "state" at the previous nucleotide position in the sequence. 
# The "state" at a sequence position is a property of that position of the sequence, 
# for example, a particular HMM may model the positions along a sequence as belonging to 
# either one of two states, GC-rich or AT-rich. 
# It is of interest because GC-rich zones are less subject to mutations (the bond is stronger).

# Consider the following HMM model that consists of two states 'GC' and 'AT'-rich as hidden states. 
# The tentative guess about the transition and emission probabilities are given by T and E below.
#           GC    AT
# T =  GC [0.8   0.2
#      AT  0.75  0.25]

#         GC   AT
# E = A [0.1  0.4   
#     C  0.4  0.1   
#     G  0.4  0.1    
#     T  0.1  0.4]  

# Given the observed sequence of nucleotides
# X = CAGCCCTAGTTGCCCCCAGAGGCAGGTAAATAGCCA
# (a) Estimate E (Emission matrix) and T (Transition matrix) from the observed sequence
# (b) What is the most probable path for generating the observed sequence, 
#     given the estimated E and T ?
# (c) For each position i in the sequence, what are the posterior probabilities of P(Yi = GC|X) and P(Yi =AT|X)

# Solution: 

States = c('GC','AT')
Symbols = c('A','C','G','T')
startProbs = c(0.5,0.5)
transProbs = matrix(c(0.80, 0.20, 
                      0.75, 0.25), nrow=2, ncol=2, byrow=TRUE)
transProbs
emissionProbs = matrix(c(0.10,0.40,
                         0.40,0.10,
                         0.40,0.10,
                         0.10,0.40), nrow=4, ncol=2,byrow=TRUE)
emissionProbs
observation = c('C','A','G','C','C','C','T','A','G','T','T','G','C','C','C','C','C','A',
                'G','A','G','G','C','A','G','G','T','A','A','A','T','A','G','C','C','A')


## Point a
# Estimate E (Emission matrix) and T (Transition matrix) from the observed sequence 

# For an initial Hidden Markov Model and a given sequence of observations, 
# the Baum-Welch algorithm infers optimal parameters to the HMM. 

hmm = initHMM(States, Symbols, startProbs, transProbs, emissionProbs)
HMM = baumWelch(hmm, observation, maxIterations=100, delta=1E-9)

HMM
# estimation of E from the observed sequence
HMM$hmm$emissionProbs

# estimation of T from the observed sequence
HMM$hmm$transProbs

# The inferred HMM is used for the next tasks.


## Point b
# What is the most probable path for generating the observed sequence, given the estimated E and T ?
# The Viterbi-algorithm computes the most probable 
# path of states for a sequence of observations for a given Hidden Markov Model.
viterbi(HMM$hmm, observation)


## Point c
# For each position i in the sequence, 
# what are the posterior probabilities of P(Yi = GC|X) and P(Yi =AT|X)
posterior(HMM$hmm, observation)

# This function computes the posterior probabilities of being in state X at time k 
# for a given sequence of observations and a given Hidden Markov Model. 


detach("package:HMM", unload = TRUE) 
# BE CAREFUL! it create conflicts with the next library we will use




