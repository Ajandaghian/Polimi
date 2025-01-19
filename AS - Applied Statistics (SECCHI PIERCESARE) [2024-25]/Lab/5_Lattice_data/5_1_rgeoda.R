#_________________ Applied Statistics 2023/2024 _________________________________

#### 5.1.Exploratory Spatial Data Analysis on Lattice Data with rgeoda ####
#_______________________________________________________________________________#

#### 1. Install and load `rgeoda` ####

# install.packages("rgeoda")

library(rgeoda)

# In addition, the package sf needs to be loaded, since it is a dependency:
  
library(sf)


#### 2. Load Spatial Data ####

# The rgeoda package for R relies on the sf package for basic spatial data
# handling functions. 
# In a typical R workflow, one first reads a shape file (or other GIS format file) 
# with the data using the sf st_read(file path) command. 
# For example, to load the Shapefile `Guerry.shp` comes with the package:
  
guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
print(guerry_path)
guerry <- st_read(guerry_path)

# Once the spatial object has been created, 
# it can be used to compute a spatial weights matrix
# using one of the several weights functions in rgeoda.


#### 3. Spatial Weights ####

# Spatial weights are central components in spatial data analysis: 
# they represent the possible spatial interactions between observations in space. 
# rgeoda provides 4 functions to create 3 different types of spatial weights:
  
# - Contiguity Based Weights: `queen_weights()`, `rook_weights()`
# - Distance Based Weights: `distance_weights()`
# - K-Nearest Neighbor Weights: `knn_weights()`


##### 3.1 Queen Contiguity Weights #####

# Contiguity means that two spatial units share a common border of non-zero length. 
# Operationally, we can further distinguish between a rook and a queen criterion of contiguity.
# The queen criterion is somewhat more encompassing and defines neighbors 
# as spatial units sharing a common edge or a common vertex.

# To create a Queen contiguity weights, one can call the function 
# queen_weights(sf_obj, order=1, include_lower_order = False, precision_threshold = 0)

# For example, to create a Queen contiguity weights using the sf object guerry:
queen_w <- queen_weights(guerry)
summary(queen_w)

# The function `queen_weights()` returns an instance of  `Weight` object. 


#### Attributes of `Weight` object

# To check if weights matrix is symmetric
is_symmetric(queen_w)

# To check if weights matrix has isolates, or if any observation has no neighbors
has_isolates(queen_w)

weights_sparsity(queen_w) # A numeric value of spatial weights sparsity

# To access the details of the weights: e.g. list the neighbors of a specified observation:
nbrs <- get_neighbors(queen_w, idx = 1)
cat("\nNeighbors of the 1-st observation are:", nbrs)

# To compute the spatial lag of a specified observation by passing the values of the selected variable:
lag <- spatial_lag(queen_w, guerry['Crm_prs']) # "Crm_prs" is "Population per Crime against persons"
head(lag)


##### 3.2 Rook Contiguity Weights #####

# The rook criterion defines neighbors by the existence of a common edge between two spatial units.
# To create a Rook contiguity weights, one can call function: 
#rook_weights(sf_obj, order=1,include_lower_order=False, precision_threshold = 0)

# For example, to create a Rook contiguity weights using the sf object guerry:
rook_w <- rook_weights(guerry)
summary(rook_w)

# The weights we created are in memory. To save the weights to a file, one can call the function:
# save_weights(gda_w, id_variable, out_path, layer_name = "")

# The id_variable defines the unique value of each observation when saving a weights file

# The layer_name is the layer name of loaded dataset. 
# For a shapefile, the layer name is the file name without the suffix (e.g. Guerry). 

# For example, using Guerry dataset, the column "CODE_DE" can be used as a key to save a weights file:
  
save_weights(rook_w, guerry['CODE_DE'], out_path = '/Users/alessandragni/Downloads/Guerry_r.gal', 
             layer_name = 'Guerry')



##### 3.3 Distance Based Weights #####

# The most straightforward spatial weights matrix constructed from a distance measure 
# is obtained when i and j are considered neighbors whenever j falls within 
# a critical distance band from i. In order to start the distance based neighbors, 
# we first need to compute a threshold value. 
# rgeoda provides a function min_distthreshold to help you find a optimized 
# distance threshold that guarantees that every observation has at least one neighbor:
# min_distthreshold(GeoDa gda, bool is_arc = False, is_mile = True)

#To create a Distance based weights, one can call the function distance_weights:
# Then, with this distance threshold, we can create a distance-band weights using the function:
#distance_weights(geoda_obj, dist_thres, power= 1.0,  is_inverse= False, is_arc= False, is_mile= True)

# For example:   
dist_thres <- min_distthreshold(guerry)
dist_thres
dist_w <- distance_weights(guerry, dist_thres)
summary(dist_w)



##### 3.4 K-Nearest Neighbor Weights #####

# A special case of distance based weights is K-Nearest neighbor weights, 
# in which every observation will have exactly k neighbors. It can be used to avoid 
# the problem of isolate in distance-band weights when a smaller cut-off distance is used. 
# To create a KNN weights, we can call the function `knn_weights`:
# knn_weights(gda, k, power = 1.0,is_inverse = False, is_arc = False, is_mile = True)

# For example, to create a 6-nearest neighbor weights using Guerry:

knn6_w <- knn_weights(guerry, 6)
summary(knn6_w)





#### 4 Local Indicators of Spatial Association–LISA ####

# rgeoda provides following methods for local spatial autocorrelation statistics:

# - Local Moran: local_moran()
# - Quantile LISA: local_quantilelisa()



##### 4.1 Local Moran #####

# The Local Moran statistic is a method to identify local clusters and local spatial outliers. 
# For example, we can call  the function `local_moran()` with the created Queen weights 
# and the data "crm_prp = guerry['Crm_prp']" as input parameters:

crm_prp = guerry["Crm_prp"]
lisa <- local_moran(queen_w, crm_prp)


# The `local_moran()` function will return a `lisa` object, and we can access its 
# values/results of lisa computation using the following functions:

# - lisa_clusters(): Get the local cluster indicators returned from LISA computation.
# - lisa_colors(): Get the cluster colors of LISA computation.
# - lisa_labels(): Get the cluster labels of LISA computation.
# - lisa_values(): Get the local spatial autocorrelation values returned from LISA computation.
# - lisa_num_nbrs(): Get the number of neighbors of every observations in LISA computation.
# - lisa_pvalues(): Get the local pseudo-p values of significance returned from LISA computation.


# For example, we can call the function `lisa_values()` to get the values of the local Moran:

lms <- lisa_values(gda_lisa = lisa)
lms


# To get the pseudo-p values of significance of local Moran computation:
  
pvals <- lisa_pvalues(lisa)
pvals


# To get the cluster indicators of local Moran computation:
  
cats <- lisa_clusters(lisa, cutoff = 0.05)
cats


# The predefined values of the indicators of LISA cluster are:

# 0 Not significant
# 1 High-High
# 2 Low-Low
# 3 High-Low
# 4 Low-High
# 5 Undefined
# 6 Isolated

# which can be accessed via the function lisa_labels():
lbls <- lisa_labels(lisa)
lbls

# By default, the local_moran() function will run with some default parameters, e.g.:
# significance_cutoff: 0.05
# permutation: 999
# permutation_method: 'complete'
# cpu_threads: 6
# seed (for random number generator): 123456789

# which are identical to GeoDa desktop software so to replicate the results in GeoDa software. 
# You can set different values when calling the lisa functions.

# For example, re-run the above local Moran example using 9999 permutations. 
lisa <- local_moran(queen_w, crm_prp, permutations = 9999)

# Then, we can use the same lisa object to get the new results after 9999 permutations:

pvals <- lisa_pvalues(lisa)
pvals


# rgeoda uses GeoDa C++ code, in which multi-threading is used to accelerate the computation of LISA. 
# We can use the argument `ncpu` to specify how many threads to run the computation:
  
lisa <- local_moran(queen_w, crm_prp, cpu_threads = 4)



##### 4.2 Quantile LISA #####

# The quantile local spatial autocorrelation converts the continuous variable 
# to a binary variable that takes the value of 1 for a specific quantile. 
# Then apply a local join count to the data converted. 
# Two input parameters, k and q, need to be specified in the function local_quantilelisa(): 
# k is the number of quantiles (k > 2), and the q is the index of selected quantile lisa ranging from 1 to k.

# For example, 

qsa <- local_quantilelisa(queen_w, crm_prp, 5, 5)

# To get the p-values and cluster indicators of the quantile LISA computation:
lisa_pvalues(qsa)
lisa_clusters(qsa)



#### 6 Exploratory Spatial Data Analysis ####

# For exploratory spatial data analysis, rgeoda provides some utility functions 
# to allow users to easily work with sf to visualize the results and do exploratory spatial data analysis.

##### 6.1 Start from sf package #####

# The sf package has been popular tool to handle geospatial data. 
# It is a good substitute of sp package which will be deprecated soon.

# For example, we can simply call plot() function to render the first 9 choropleth maps 
# using the first 9 variables in the dataset:
  
plot(guerry)


##### 6.2 Exploratory Spatial Data Analysis with rgeoda #####

# Now, with the sf object guerry, you can call rgeoda's spatial analysis functions. 
# For example, to examine the local Moran statistics of variable "crm_prs" 
# (Population per Crime against persons):

queen_w <- queen_weights(guerry)
lisa <- local_moran(queen_w,  guerry['Crm_prs'])


##### 6.3 Create Local Moran Map #####

# With the LISA results, we can make a local Moran cluster map:
lisa_colors <- lisa_colors(lisa)
lisa_labels <- lisa_labels(lisa)
lisa_clusters <- lisa_clusters(lisa)
plot(st_geometry(guerry), 
     col=sapply(lisa_clusters, function(x){return(lisa_colors[[x+1]])}), 
     border = "#333333", lwd=0.2)
title(main = "Local Moran Map of Crm_prs")
legend('bottomleft', legend = lisa_labels, fill = lisa_colors, border = "#eeeeee")

# In the above code, we use th values of cluster indicators from rgeoda's LISA object 
# are used to make the LISA map. We can save the clusters back to the original sf data.frame:

guerry['moran_cluster'] <- lisa_clusters


# Checking the values of the cluster indicators, we will see they are integer numbers 
# 0 (not significant), 1 (high-high cluster), 2 (low-low cluster), 3 (low-high cluster), 
# 4 (high-low cluster), 5 (neighborless/island), 6 (undefined):

lisa_clusters

# To create a significance map that is associated with the local Moran map, 
# we can do the same as making the local Moran cluster map using the results from lisa_pvalues():
  
lisa_p <- lisa_pvalues(lisa)
p_labels <- c("Not significant", "p <= 0.05", "p <= 0.01", "p <= 0.001")
p_colors <- c("#eeeeee", "#84f576", "#53c53c", "#348124")
plot(st_geometry(guerry), 
     col=sapply(lisa_p, function(x){
       if (x <= 0.001) return(p_colors[4])
       else if (x <= 0.01) return(p_colors[3])
       else if (x <= 0.05) return (p_colors[2])
       else return(p_colors[1])
     }), 
     border = "#333333", lwd=0.2)
title(main = "Local Moran Map of Crm_prs")
legend('bottomleft', legend = p_labels, fill = p_colors, border = "#eeeeee")


