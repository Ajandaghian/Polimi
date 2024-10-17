#______________ Applied Statistics 2024/2025 _________________

#### 1.2 Recap on Linear Models formulation (Chapter 5.2) ####
#____________________________________________________________#

# We specify the Mean Structure using a Model Formula
# y (dependent variable) ~ term1 + term2 +... (mean structure of the model)
# (i.e. y = beta_0 + beta_1 * term1 + ... + eps )

# Operators used when specifying an R formula:
# + (Essential)                       Separates terms in the formula
# : (Essential)                       Separates predictors in interaction terms
# *, /, %in%, -, ^ (Non essential)    Used to keep the formula short

# Symbolic operations recap:
# let us consider:
# - dependent variable named y;
# - explanatory covariates include:
#     - three continuous variables (x1, x2, x3);
#     - three factors (f1, f2, f3).

# y ~ x1                       # Univariate linear regression 
# formula(y ~ x1)              # ... equivalent specification

## NB: an intercept is implicitly included by default.
## To explicitly specify the inclusion of an intercept in the model, we use 1 as a separate term
# y ~ 1 + x1                   # Explicit indication for intercept 

## NB: to indicate that there is no intercept in the model, 
## we can use 0 or -1 as a separate term
# y ~ 0 + x1                   # No intercept using term 0
# y ~ -1 + x1                  # No intercept using term -1

# y ~ f1 + x1                  # ANCOVA with main effects only 
# y ~ f1 + x1 + f1:x1          # Main effects and factor by numeric interaction   
# y ~ f1 + f2 + f1:f2          # Main effects and f1 by f2 two way interaction 
# y ~ f1 + f1:f3               # f3 nested within f1  
# y ~ x1 + f1 + f2 +           # Main effects and ... 
#   x1:f1+ x1:f2 + f1:f2       # ... two-way interactions 


# Formulae employing nonessential operators (*, /, %in%, -, ^)
# y ~ f1*f2             means     y ~ f1 + f2 + f1:f2       # ANOVA with two-way interaction
# y ~ f1 + f3 %in% f1   means     y ~ f1 + f1:f3            # f3 nested within f1 
# y ~ f1/f3             means     y ~ f1 + f1:f3            # ... equivalent specification
# y ~ (x1+f1+f2)^2      means     y ~ f1 + f2 + f3 + 
#                                   f1:f2 + f1:f3 + f2:f3   # Up to 2nd order interactions

# Composite terms:
# -> Formulae syntax can be extended through mathematical functions, e.g.
# y ~ sqrt(x1) + x2     # Square root transformation of x1
# log(y) ~ x1 + x2      # log transform for y

# -> the functions I() and update():
# form2 = y ~ I(x1 + x2/100)   # I() function for to explicitate the arithmetic (vs default) meaning 
# update(form2, . ~ . + x3)    # x3 predictor added to form2
# update(form2, . ~ . -1)      # Intercept omitted from form2

