#_______________________ Applied Statistics 2024/2025 _____________________

#### 1.6 Additional exercises: LMs with Dummies and Hypothesis Testing ####
#__________________________________________________________________________

# SET THE DIRECTORY!

library(MASS)
library(car)


##### EXERCISE 1 #####

# The file vehicles.txt reports the number Y (expressed in thousands of units)
# of vehicles registered annually in three countries of the European Union
# (France, Germany and Italy) during a reference period of 10 years.
# Recent economic models describe the behavior of this variable according
# the model:
# Y | (X = x, G = g) = beta0.g + beta1.g * x^2 + eps
# with eps ~ N (0, sigma ^ 2), x = 1, 2,. . . , 10 (year) and
# g = France, Germany, Italy (EU country).
# (a) With the method of least squares, estimate the 7 parameters of the model.
# (b) using appropriate statistical tests, state if you deem necessary to
#     include into the model:
#     1. the variable x^2;
#     2. the variable G;
#     3. the effect of the variable G onto the coefficient that multiplies the
#        regressor x^2;
#     4. the effect of the variable G on the intercept.
# (c) Once identified the "best model", build three prediction intervals
#     for the number of vehicles registered in the three countries 
#     during the eleventh year, so that the three new observations
#     will fall simultaneously within the respective ranges with 95%
#     of probability.

pb4  <- read.table('vehicles.txt')
pb4

matplot(pb4, type='l',lwd=2, xlab='year', ylab='(thousands of) vehicles')


### question (a)

# We first build the design matrix and the vector of the responses
Anno <- rep(1:10,3)
Anno

Imm <- c(pb4[,1], pb4[,2], pb4[,3])
Imm

# Model: Imm = beta0.g + beta1.g*Anno^2 + eps (Anno=Year)
# con g=0,1,2 [Italy,France,Germany]; E[eps]=0, Var(eps)=sigma^2

# We need to build appropriate dummies to account for the levels of
# the categorical variable G=Italy/France/Germany (3 levels)
# g groups => g-1 dummies (3 groups => 2 dummies)

dFr <- rep(c(1,0), c(10,20))      # dFr = 1 if France,  0 otherwise
dGer<- rep(c(0,1,0), c(10,10,10)) # dGer= 1 if Germany, 0 otherwise
dFr
dGer

# Equivalent Model:
# Imm = b.0 + b.1*dFr + b.2*dGer + b.3*Anno^2 + b.4*dFr*Anno^2 + b.5*dGer*Anno^2

# Indeed:
# beta0.It=b.0;      beta1.It=b.3;
# beta0.Fr=b.0+b.1;  beta1.Fr=b.3+b.4;
# beta0.Ger=b.0+b.2; beta1.Ger=b.3+b.5

dati <- data.frame(Imm   = Imm,
                   Anno2 = Anno^2,
                   dFr   = rep(c(1,0), c(10,20)),      # dummy that selects France
                   dGer  = rep(c(0,1,0), c(10,10,10))) # dummy that selects Germany
dati

fit <- lm(Imm ~ dFr + dGer + Anno2 + Anno2:dFr + Anno2:dGer, data=dati)


summary(fit)
summary(fit)$sigma

### question (b)
shapiro.test(residuals(fit))
x11()
par(mfrow=c(2,2))
plot(fit)

dev.off()

# 1. the variable x^2;
linearHypothesis(fit,
                 rbind(c(0,0,0,1,0,0),
                       c(0,0,0,0,1,0),
                       c(0,0,0,0,0,1)),
                 c(0,0,0))

# 2. the variable G;
linearHypothesis(fit,
                 rbind(c(0,1,0,0,0,0),
                       c(0,0,1,0,0,0),
                       c(0,0,0,0,1,0),
                       c(0,0,0,0,0,1)),
                 c(0,0,0,0))

# 3. the effect of the variable G onto the coefficient that multiplies the regressor x^2;
linearHypothesis(fit,
                 rbind(c(0,0,0,0,1,0),
                       c(0,0,0,0,0,1)),
                 c(0,0))

# 4. the effect of the variable G on the intercept.
linearHypothesis(fit,
                 rbind(c(0,1,0,0,0,0),
                       c(0,0,1,0,0,0)),
                 c(0,0))

#- Interaction terms with `Anno2` and `G` are statistically significant.
#- The effect of `G` on the coefficient of `x²` is significant.
#- `G` does not significantly affect the intercept.
#- The hypothesis tests for `x²` and `G` interactions are strongly rejected.

### question (c)
fit2 <- lm(Imm ~ Anno2 + Anno2:dFr + Anno2:dGer, data=dati)
summary(fit2)

nuovi <- data.frame(Anno2 = c(11,11,11)^2, dFr=c(1,0,0), dGer=c(0,1,0))
IP <- predict(fit2, newdata=nuovi, interval='prediction', level=1-0.05/3)
rownames(IP) <- c('Fr','Ger','It')
IP



##### EXERCISE 2 #####

# At the Tenaris steel mills, the relationship between length [m] and
# Temperature [° C] of some steel bars that will be sold to Pirelli
# is under study (the data are contained in tenaris.txt file). The relation
# is hypothesized of the kind:
#   L = L0 + C* T + D  * T ^ 2 + eps
# with L the length of the bar, T the temperature of the bar, L0 the length 
# of the bar at 0 °C, C the coefficientof of linear thermal expansion, D
# the coefficient of quadratic thermal expansion and eps a measurement error
# of zero mean.
# Answer the following questions using appropriate statistical arguments:
# a) Estimate the parameters L0, C, D and the variance of error eps.
# b) Based on the analysis of residuals, do you think that there are the
#    conditions to make inference on the coefficients based on a Gaussian
#    model? (In case of Yes proceed to step (c); in case of negative answer
#    identify the problem, remove it and return to point (a))
# c) Do you think that the model explains the possible dependence between 
#    the temperature T and the length L?
# d) do you deem plausible to consider that the length of the bars at 0 °C
#    is equal to 2?
# E) do you think that you can eliminate from the model the quadratic term?

ten <- read.table('tenaris.txt', header=TRUE)
ten

attach(ten)

### question a)

fit <- lm(L ~ T + I(T^2))

summary(fit)

summary(fit)$sigma^2
# oppure
e <- residuals(fit)
S2 <- t(e)%*%e / (df.residual(fit))
S2

### question b)

shapiro.test(e)

x11()
par(mfrow=c(2,2))
plot(fit)

x11()
plot(T,L)
points(ten[1,1],ten[1,2],pch=19)

graphics.off()

detach(ten)

# Remove the outlier
ten1 <- ten[-1,]

fit <- lm(L ~ T + I(T^2), data=ten1)
summary(fit)

e <- residuals(fit)
S2 <- t(e)%*%e / (df.residual(fit))
S2

shapiro.test(e)

x11()
par(mfrow=c(2,2))
plot(fit)

dev.off()

### question c)

attach(ten1)
x11()
plot(T,L)
points(T,fitted(fit),col='blue', pch=19)

dev.off()

summary(fit)$r.squared

### question d)

linearHypothesis(fit, c(1,0,0), 2)

### question e)
summary(fit)

# or
linearHypothesis(fit, c(0,0,1), 0)

#_______________________________________________________________________________
