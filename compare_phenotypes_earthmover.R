# script to compare the distribution of mitographs in 
# wt and double delete mutants using Earthmover distance

# first you run extract_n_p_largest_component_V3 separately for WT and for DD
# name the output csv files results_WT.csf and results_DD.csv

# then the program will:
# - set up an array indexed with a and b, that is the number of outgrowth and fusion events
# as per alternative labeling scheme 3 in the paper.

# - for each data record use n and p to find a and b and then add one to that array entry
# do that for each genotype

# we now have two arrays that we can compare using earthmover distance



# script to compute the [n,p] representation of the largest component of all mitographs in a directory
# detects graphs that aren't mitographs, i.e. they have some vertices with degree other than 1 or 3
# and throws them out since they have something wrong with them


# before running make sure to have the packages installed and set the current working directory to where all the .gnet files are
# install.packages("emdist")
# library(emdist)
#setwd("/Users/wallacemarshall/papers/Mitochondria_graph_theory/data")

# to run the script just click the Source button!

datarecords_WT <- read.table("results_WT.csv", sep=",", header=TRUE)

datarecords_DD <- read.table("results_DD.csv", sep=",", header=TRUE)

# count how many datapoints there are because we will need to convert to a density
datapoints_WT = nrow(datarecords_WT)
datapoints_DD = nrow(datarecords_DD)

# determine the size of the array needed to hold all the data
maxn1 = max(datarecords_WT$N, na.rm = TRUE)
maxn2 = max(datarecords_DD$N, na.rm = TRUE)
maxp1 = max(datarecords_WT$P, na.rm = TRUE)
maxp2 = max(datarecords_DD$P, na.rm = TRUE)

#nmax = max(maxn1,maxn2)
#pmax = max(maxp1,maxp2)

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
dataset_DD <- matrix(0, nrow=amax+1, ncol=bmax+1)

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

# iterate over all rows in the DD data file
for (data_index in 1:datapoints_DD)
{
  current_p =datarecords_DD[data_index, 2]
  current_n =datarecords_DD[data_index, 3]
  
  # convert from [p,n] to (a,b) labeling 
  # i.e. the number of branch outgrowths and fusion events for a given class
  b = current_n
  a = (current_n + 2 - current_p)/2
  
  if (a <= amax && b <= bmax)
  {
    # update count for that value of a,b
    dataset_DD[a,b] = dataset_DD[a,b] + 1
  }
}


# now calculate the density ie the fraction of counts in each state
dataset_WT_density = dataset_WT / sum(dataset_WT)
dataset_DD_density = dataset_DD / sum(dataset_DD)

# find earthmover distance using manhattan distance function
result <- emd2d(dataset_WT_density, dataset_DD_density, dist = "manhattan")

print(result)



# now do permutation test
# -------------------------------------------------------------------
# 1. Helper function: Convert a density matrix back to (X, Y) coordinates
# -------------------------------------------------------------------
matrix_to_points <- function(mat, n_points = 5000) {
  # Get row and col indices
  grid <- expand.grid(x = 1:nrow(mat), y = 1:ncol(mat))
  weights <- as.vector(mat)
  
  # Sample bin indices according to the matrix density weights
  sampled_indices <- sample(seq_len(nrow(grid)), size = n_points, replace = TRUE, prob = weights)
  
  return(grid[sampled_indices, ])
}

# -------------------------------------------------------------------
# 2. Helper function: Bin 2D points back into a normalized matrix
# -------------------------------------------------------------------
points_to_matrix <- function(pts, nr, nc) {
  # Tabulate x and y into a 2D matrix grid
  mat <- table(factor(pts$x, levels = 1:nr), factor(pts$y, levels = 1:nc))
  mat <- as.matrix(mat)
  
  # Normalize to ensure it sums to 1
  return(mat / sum(mat))
}


# -------------------------------------------------------------------
# 3. Main Permutation Test Function
# -------------------------------------------------------------------
perm_test_emd <- function(mat_A, mat_B, n_permutations = 1000, n_simulated_points = 5000) {
  nr <- nrow(mat_A)
  nc <- ncol(mat_A)
  
  # Ensure inputs are normalized
  mat_A <- mat_A / sum(mat_A)
  mat_B <- mat_B / sum(mat_B)
  
  # Step A: Compute Observed EMD
  obs_emd <- emd2d(mat_A, mat_B, dist = "manhattan")
  
  # Step B: Reconstruct point clouds representing each density grid
  pts_A <- matrix_to_points(mat_A, n_points = n_simulated_points)
  pts_B <- matrix_to_points(mat_B, n_points = n_simulated_points)
  
  # Combine points into a single pooled pool
  pooled_pts <- rbind(pts_A, pts_B)
  n_total <- nrow(pooled_pts)
  n_A <- nrow(pts_A)
  
  # Step C: Permutation Loop
  perm_emds <- numeric(n_permutations)
  
  for (i in 1:n_permutations) {
    # Randomly assign indices to group A and group B
    shuffled_idx <- sample.int(n_total)
    
    perm_pts_A <- pooled_pts[shuffled_idx[1:n_A], ]
    perm_pts_B <- pooled_pts[shuffled_idx[(n_A + 1):n_total], ]
    
    # Bin points back into matrices
    perm_mat_A <- points_to_matrix(perm_pts_A, nr, nc)
    perm_mat_B <- points_to_matrix(perm_pts_B, nr, nc)
    
    # Calculate EMD for the null distribution
    perm_emds[i] <- emd2d(perm_mat_A, perm_mat_B, dist = "manhattan")
  }
  
  # Step D: Calculate p-value (adding +1 avoids p = 0)
  p_value <- (sum(perm_emds >= obs_emd) + 1) / (n_permutations + 1)
  
  return(list(
    observed_emd = obs_emd,
    p_value = p_value,
    null_distribution = perm_emds
  ))
}


# Run permutation test
results <- perm_test_emd(dataset_WT_density, dataset_DD_density, n_permutations = 1000, n_simulated_points = 5000)

cat("Observed EMD:", results$observed_emd, "\n")
cat("P-value:", results$p_value, "\n")

# Plot null distribution vs observed value
hist(results$null_distribution, main = "Null Distribution of Permuted EMDs", xlab = "EMD")
abline(v = results$observed_emd, col = "red", lwd = 2, lty = 2)




