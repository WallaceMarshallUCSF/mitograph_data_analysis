# script to test for complate spatial randomness in mitograph data using Moran's I test

# revised verison that uses brute force to form the matrix

# first you run extract_n_p_largest_component_V3 for WT 
# name the output csv file results_WT.csv

# then the program will:
# - set up an array indexed with a and b, that is the number of outgrowth and fusion events
# as per alternative labeling scheme 3 in the paper.



# before running make sure to have the packages installed and set the current working directory to where all the .gnet files are
# install.packages("terra", type = "mac.binary")
# library(terra)
# install.packages("spdep", type = "mac.binary")
# library(spdep)
# library(sf)
#setwd("/Users/wallacemarshall/papers/Mitochondria_graph_theory/data")

# to run the script just click the Source button!

datarecords_WT <- read.table("results_WT.csv", sep=",", header=TRUE)

# count how many datapoints there are because we will need to convert to a density
datapoints_WT = nrow(datarecords_WT)

# determine the size of the array needed to hold all the data
#maxn1 = max(datarecords_WT$N, na.rm = TRUE)
#maxp1 = max(datarecords_WT$P, na.rm = TRUE)

# just analyze data where p and n are 10 or less
# to avoid having huge earth mover distances for the sparsely sampled
# data with larger n or p
nmax = 10
pmax = 10

bmax = nmax
amax = (nmax + 2)/2

# make arrays to hold the counts
# row or column 1 corespond to n or p equal to zero
dataset_WT <- matrix(0, nrow=amax+1, ncol=bmax+1)

# iterate over all rows in the WT data file
for (data_index in 1:datapoints_WT)
{
  current_p =datarecords_WT[data_index, 2]
  current_n =datarecords_WT[data_index, 3]
  
  # convert from [p,n] to (a,b) labeling 
  # i.e. the number of branch outgrowths and fusion events for a given class
  # add 1 to each value because zero is not allowed as an index
  b = current_n + 1
  a = (current_n + 2 - current_p)/2 + 1
  
  if (a <= amax && b <= bmax)
  {
    # update count for that value of a,b
    dataset_WT[a,b] = dataset_WT[a,b] + 1
  }
}



counts_matrix <- dataset_WT


# 2. Extract dimensions automatically
r_len <- nrow(counts_matrix)
c_len <- ncol(counts_matrix)

# 3. Build grid neighbors
grid_coords <- cell2nb(nrow = r_len, ncol = c_len, type = "queen")

# 4. Create weights matrix
weights_matrix <- nb2listw(grid_coords, style = "W")

# 5. Flatten matrix row-by-row into a 1D vector
cell_values <- as.vector(t(counts_matrix))

# Sanity Check: Ensure lengths match before testing
stopifnot(length(cell_values) == length(grid_coords))

# 6. Run Moran's I test
moran_result <- moran.test(cell_values, weights_matrix)

# Print the test output
print(moran_result)



