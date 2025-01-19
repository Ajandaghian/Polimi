#_____________ Applied Statistics 2024/2025 _______________________
#### 3.3 Illustration of Correlation Structures (Chapter 11.4) ####
#_________________________________________________________________#

# DATA GENERATION
# hypothetical data frame df containing two subjects: 
# -The first subject has four consecutive observations (variable occ)
#  The observations are made at different locations in a two-dimensional space, 
#  with the coordinates given by the loc1 and loc2 position variables. 
# -The second subject has only three observations, with the
#  observation for the third occasion missing. 

#  Note that, for the two subjects, the
#  coordinates of the observations made at the same occasion differ.

subj <- rep(1:2, each = 4)               # Two subjects
occ  <- rep(1:4, 2)                      # Four observations each
loc1 <- rep(c(0, 0.2, 0.4, 0.8), 2)      # First coordinate
loc2 <-                                  # Second coordinate 
  c(0, 0.2, 0.4, 0.8, 0, 0.1, 0.2, 0.4) 
df0  <-                                  # Hypothetical data frame
  data.frame(subj, occ, loc1, loc2)
(df  <-                                  # Occ = 3 for subj.2 deleted
    subset(df0, subj != 2 | occ != 3))    



#### 1) Compound symmetry: corCompSymm class (SERIAL) ####

# h(k,rho)=rho  where k=|j-j'|
# i.e. constant correlation coefficient between any two observations
cs <- corCompSymm(value = 0.3, form = ~1|subj)  # Object defined with constant correlation = 0.3
cs <- Initialize(cs, df)                        # initialized for the dataset df
# same correlation structure (with different dimensions) is obtained for both subjects.

# NB. providing corCompSymm(value = 0.3, form = ~occ|subj) would lead to the same object,
# because in corCompSym any position variable (time t in this case) is ignored;
# REMEMBER: h(k,rho) = h(d(t_ij, t_ij'), rho) and in this case: h(k,rho)=rho

getCovariate(cs)                   # Positions in series
corMatrix(cs, corr=TRUE)           # Corr. matrix C_i

?corMatrix.corStruct



#### 2) Autoregressive structure of order 1: corAR1 class (SERIAL) #### 

# h(k,rho)=rho^k   k=|j-j'|
cs1 <- corAR1(0.3, form = ~tx)         # Uninitialized corAR1 struct

tx  <- 1:4                             # A covariate with values 1, 2, 3, 4
df2 <- data.frame(tx)                  # An auxiliary data frame
cs1i <- Initialize(cs1, data = df2)    # Initialized corAR1 object

corMatrix(cs1i, corr=TRUE)             # Corr. matrix C_i



# NB. [ISSUES regarding serial correlation classes other than corCompSymm]
# form = ~1|group means "use of the order of the observations in the group as the position index"
# When data are balanced (i.e., all subjects have all measurements, or have monotone
# missingness patterns), this works fine. 
# However, if, for some subjects, intermittent measurements are missing, 
# the use of the observation order can result in the wrong correlation structure.

# EXAMPLE
df
# Not-recommended syntax ...
car <- corAR1(value = 0.3, form = ~1|subj) 
carI <- Initialize(car, df)       # corAR1 class object initialized
getCovariate(carI)                # Position=order of observations for a subject     
corMatrix(carI)[[1]]              # Correct matrix for the 1st subject
corMatrix(carI)[[2]]              # Incorrect matrix for the 2nd subject
# for example, in position [1,3]
# we get rho^(3-1) = rho^2 = 0.09 cause we are assuming distance between first and third observation
# but we should have rho^(4-1) = rho^3 = 0.027

# Recommended syntax!!!
car1 <- corAR1(value = 0.3, form = ~occ|subj)   
car1 <- Initialize(car1, df)      # corAR1 classs object initialized  
getCovariate(car1)                # Correct positions based on the occ variable
corMatrix(car1)[[2]]              # Correct matrix for the 2nd subject

# For data with measurement timepoints common to all subjects, this
# caution is required only for nonmonotone missing data patterns. 
# Nevertheless, it is prudent to always use a position variable, which reflects the proper
# positions of the observations in a sequence for each group (subject), in the form argument.




#### 3) Exponential structure: corExp class (SPATIAL) #### 

# h(s,rho)=exp(-s/rho)   s=distance; rho=range

# case 1
ceE <- corExp(value = 1, form = ~ loc1 + loc2 | subj)            # Euclidean metric (default)
ceE <- Initialize(ceE, df)     
corMatrix(ceE)                    # List with corr matrices for both subjects
# Note that the matrices differ, because spatial coordinates of the measurements differ for the subjects.


# case 2
ceM <- corExp(1, ~ loc1 + loc2 | subj, metric = "man")          # Manhattan metric
ceM <- Initialize(ceM, df)
corMatrix(ceM)[[1]]              # Corr matrix for the 1st subject
# note that this matrix is different from the case 1 (for subject 1)


# case 3
ceEn <- corExp(c(1, 0.2), ~ loc1 + loc2 | subj, nugget = TRUE)  # nugget = 0.2 (euclidean metric)
# value = c(1, 0.2) (i.e., range=1, nugget=0.2)
ceEn <- Initialize(ceEn, df)
corMatrix(ceEn)[[1]]             # Corr matrix for the 1st subject
