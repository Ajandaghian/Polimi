#______________________ Applied Statistics 2024/2025 _________________

####  1.3 Using lm() VS gls() to fit a Linear Model (Chapter 5.4) ####
#____________________________________________________________________#

?lm()
?gls()
#_Component_________________lm() (stats)_____________gls() (nlme)____________
# Formula                   formula =                model =
# Data                      data =                   data =
# Subset                    subset =                 subset =
# Missing values            na.action = na.fail      na.action = na.fail 
# Estimation method         method="qr"              method= "REML" (restricted likelihood)
# Offset                    offset =                 –

# NB. An offset is a known additional term x^0_i for all i.
#     y_i = x^(0)_i + beta_1 * x^(1)_i + ... + beta_p * x^(p)_i + eps_i
#     It is included as an additional (the first) column of the design matrix, and the
#     corresponding parameter beta0 is assumed to be known and equal to 1.
#     An offset can be accommodated in the classical LM defining \tilde{y_i} = y_i - x^(0)_i
#     without the need for an explicit modification of the design matrix.



#_____________________________________________________________________________________#
##### Extracting Information from a Model-Fit Object: lm() VS gls() (Chapter 5.5) #####
#_____________________________________________________________________________________#

#_Component_________________lm() (stats)____________________gls() (nlme)______________
# Summary                   (summ <- summary(lm.fit))       (summ <- summary(gls.fit))
# Est. method                                               gls.fit$method
# beta_hat                  coef(lm.fit)                    coef(gls.fit)
# ^beta, se(^beta), t-test  coef(summ) or summ              coef(summ) or summ
# ^Var(^beta)               vcov(lm.fit)                    vcov(gls.fit)
# 95% CI for ^beta          confint(lm.fit)                 confint(gls.fit)
#                                                           intervals(gls.fit, which="coef")
# ^sigma                    summ$sigma                      summ$sigma
# 95% CI for ^sigma                                         intervals(gls.fit, which="var-cov")

# ML value                  logLik(lm.fit)                  logLik(gls.fit, REML = FALSE)
# REML value                logLik(lm.fit, REML = TRUE)     logLik(gls.fit, REML = TRUE)
# AIC                       AIC(lm.fit)                     AIC(gls.fit)
# BIC                       BIC(lm.fit)                     BIC(gls.fit)

# Fitted values             fitted(lm.fit)                  fitted(gls.fit)
# Raw residuals             residuals(lm.fit,               residuals(gls.fit,
#                                    type="response")                type="response")
# Predicted                 predict(lm.fit, newdata)        predict(gls.fit, newdata)

# R-call                    (cl <- getCall(lm.fit))         (cl <- getCall(gls.fit))
# Formula for mean          (form <- formula(lm.fit))       (form <- formula(gls.fit))
# Data name                 (df.name <- cl$data)            (df.name <- cl$data)
# Data frame                eval(df.name)                   eval(df.name)
# Model frame               (mf <- model.frame(lm.fit))     mfDt <- getData(gls.fit)
#                                                           (mf <- model.frame(form, mfDt))
# Design matrix             model.matrix(lm.fit)
#                           model.matrix(form, mf)          model.matrix(form, mf)

