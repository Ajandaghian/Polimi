#______________________________ Applied Statistics 2024/2025 ________________________________

#### 1.5 Supplementary Materials ####
#____________________________________
# Summary:
# 1.5.1 Some details about Contrast matrices
# 1.5.2 From a Formula to the Design Matrix
  

##### 1.5.1 Some details about Contrast matrices (Section 5.3.2.2) #####
#____________________________________________________________________________________________

?factor

# factor() or ordered() functions are used to create unordered and ordered factors, respectively. 

# To decode a given factor into the columns of a design matrix, 
# it is necessary to associate the factor with an appropriate matrix of contrasts.

options()$contrasts             # Default contrasts
# The advantage of using the default choice of contrast matrices is that, 
# in studies with a balanced design, columns of the design matrix X become orthogonal. 
# Consequently, the matrix X'X (used in the estimation of the parameters beta), 
# becomes diagonal and the corresponding estimates of the elements of beta are uncorrelated. 

# Other choices of contrasts in LMs may introduce artificial correlations 
# between the estimates, even for balanced designs. 

# The choice of the contrasts involves a trade-off between simplicity and 
# the interpretability of the beta estimates.

# Some examples:
# 1) Dummy (or Treatment) coding 
#    It's the coding scheme we are most familiar with.
#    It compares each level of a categorical variable to a reference level. 
#    By default, the reference level is the 1st level of the categorical variable (alphabetical order).
#    Creates a contrast column for each level except the base level. 
#    The column for level k has a 1 at row k and 0 elsewhere.
contr.treatment(3)              # Default base level = 1 
contr.treatment(3, base = 3)    # base level = 3. Same as contr.SAS(3)

# 2) Helmert contrasts
#    produces orthogonal contrasts
#    contrast i is the difference between level i+1 and the average of levels 1:i.
#    it contrasts the second level with the first, the third with the average of the first two, and so on 
contr.helmert(3)                # Helmert contrasts 

# 3) Sum-to-zero contrasts
#    produces non-orthogonal contrasts and returns coefficients s.t. their sum is zero
contr.sum(3)                    # Sum-to-zero

# 4) Polynomial contrast
#    produces orthogonal contrasts and creates orthogonal polynomials of degree 1, 2, etc., 
#    when evaluated at scores. Columns are scaled to have Euclidean norm 1
contr.poly(3)                   # Contrasts based on orthogonal polynomials
contr.poly(3, scores=c(1,5,7))  # Polynomial contrasts


# Assigning and extracting a set of contrasts for a factor:
lesion.f <- factor(armd.wide$lesion)   # Factor created
str(lesion.f)                          # Structure
levels(lesion.f)                       # Levels extracted
contrasts(lesion.f)                    # Contrasts extracted


# Assigning contrasts using the “contrasts() <- contrast function” syntax
contrasts(lesion.f) <- contr.sum(4)    # ... "contrasts <- " syntax
# See in 1_1_LMs_homosc_indep_errors_ARMD.R how it is performed

# NB. All the contrast matrices have k rows and k−1 columns, 
# where k is the number of levels of the corresponding factor. 
# In this way, we avoid collinearity in the design matrices containing an intercept. 
# example: Treatment coding: if the values are "Red", "White" and "Blue";
# we might have a column named "Red," containing 1's for Red items and 0 for non-Red items, 
# and a column named "White," containing 1's for White items and 0 for non-White items. 
# That's all we need -- if we see an item with 0's for both Red and White, it must be blue. 
# Red   1  0 (0)
# White 0  1 (0)
# Blue  0  0 (1)

# More generally, by assigning a contrast matrix with at most k−1 linearly independent columns, 
# we avoid collinearity in a design matrix for any model containing all terms 
# of lower order than a given factor or, more broadly, a term involving factor(s). 

# However, in some cases like, e.g., of a model without an intercept, 
# it is more appropriate to use a k×k identity matrix instead of a k×(k −1) contrast. 
# Such a choice is possible using the contrasts=FALSE argument of the contrasts() function.



#____________________________________________________________________________________________
##### 1.5.2 From a Formula to the Design Matrix (Chapter 5.3) #####

# model.frame = model.frame(formula, data.frame)
# design.matrix = model.matrix(formula, model.frame)

data(armd.wide, package = "nlmeU")

# *** Formula  ***
form1 <- formula(visual52 ~                     # Dependent variable 
                   sqrt(line0) +                # Continuous explanatory variable
                   factor(lesion) +             # Factor with 4 levels
                   treat.f*log(visual24) +      # Crossing of two variables
                   visual0)
form1

form2 = visual52 ~ sqrt(line0) + factor(lesion) +  treat.f*log(visual24) +  visual0
form2 # equivalent to form1


# *** Model frame ***
armd.mf1 <- model.frame(form1,                           # Formula
                        data = armd.wide,                # Data frame
                        subset = !(subject %in% c(1,2)), # Exclude two subjects   
                        na.action = na.exclude,          # Dealing with missing data (excluding them)
                        SubjectId = subject)             # Identifier of data records
class(armd.mf1)           
dim(armd.wide)                  # Data frame dimensions
dim(armd.mf1)                   # Model frame dimensions
names(armd.mf1)                 # Components of the model frame
head(armd.mf1, n = 4)           # First four records


# *** Model matrix ***
Xmtx <-  model.matrix(form1, armd.mf1)         # DESIGN MATRIX
dim(Xmtx)                                      # No rows and cols
print(head(Xmtx, n = 6), digits = 4)           # First 6 rows
length(unique(armd.wide$lesion))
names(attributes(Xmtx))                        # Attribute names
attr(Xmtx,"assign")                            # Cols to terms map
attr(Xmtx,"contrasts")                         # Contrasts attribute


# From the code in 1_1_LMs_homosc_indep_errors_ARMD.R 
# we can compute the design matrix and the contrasts for the (unordered) factors
# NB: model.frame   = model.frame (formula, data.frame)
#     design.matrix = model.matrix(formula, model.frame) 
data(armd)

lm1.form <- formula(visual ~ -1 + visual0 + time.f + treat.f:time.f ) 

vis.lm1.mf <- model.frame(lm1.form, data = armd)      # Model frame
vis.lm1.dm <- model.matrix(lm1.form, vis.lm1.mf)      # Design matrix X
vis.lm1.dm
dim(vis.lm1.dm)                                       # Dimensions
(nms <- colnames(vis.lm1.dm))                         # Long column names ...
nms <- abbreviate(nms)                                # ... abbreviated
colnames(vis.lm1.dm) <- nms                           # ... assigned.
head(vis.lm1.dm, n = 6)                               # X matrix. Six rows.
attr(vis.lm1.dm, "contrasts")                         # Contrasts attribute (commented above)
contrasts(armd$treat.f)                               # Contrasts for treat.f
