#______________________ Applied Statistics 2024/2025 _________________#

####  1.4 Anscombe Quartet ####
#_____________________________#

# This example show some cases for which examining data and residuals is crucial! 
# (don't look at R2 only!)

anscombe
# Four x-y datasets which have the same traditional statistical 
# properties (mean, variance, correlation, regression line, etc.), 
# yet are quite different.

attach(anscombe)

# dataset 1
lm1 <- lm(y1 ~ x1)
summary(lm1)

x11(width=14, height=7)
par(mfcol=c(2,4))
plot(x1,y1, main='Dataset 1')
abline(lm1)

plot(x1,residuals(lm1))
abline(h=0)

# dataset 2
lm2 <- lm(y2 ~ x2)
summary(lm2)

plot(x2,y2, main='Dataset 2')
abline(lm2)

plot(x2,residuals(lm2))
abline(h=0)

# same R^2, same coefficient estimate, same residual std error

# dataset 3
lm3 <- lm(y3 ~ x3)
summary(lm3)

plot(x3,y3, main='Dataset 3')
abline(lm3)

plot(x3,residuals(lm3))
abline(h=0)

# dataset 4
lm4 <- lm(y4 ~ x4)
summary(lm4)

plot(x4,y4, main='Dataset 4')
abline(lm4)

plot(x4,residuals(lm4))
abline(h=0)

graphics.off()
detach(anscombe)
