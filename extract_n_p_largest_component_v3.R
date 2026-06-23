# script to compute the [n,p] representation of the largest component of all mitographs in a directory
# detects graphs that aren't mitographs, i.e. they have some vertices with degree other than 1 or 3
# and throws them out since they have something wrong with them

# modified from v2 to also read in a csv file of metadata and separate the files into wt versus DD (double deletion)

#NOTE ***********
# you have to give it a value for a parameter called "target_genotype" which is either WT or DD
# and then it will only analyze the one you give it.  if target_genotype is ALL then it will
# analyze everything

# before running make sure to have the packages installed and set the current working directory to where all the .gnet files are
#install.packages("igraph")
#library(igraph)
#library('dplyr')
#library(ggplot2)
#setwd("/Users/wallacemarshall/papers/Mitochondria_graph_theory/data")

# before running load the data in like this:
# datarecords <- read.table("MitoTable.csv", sep=",", header=TRUE)
# and then tell it what genotype to select like this:
# target_genotype <- "DD"
# to run the script just click the Source button!

if (target_genotype == "WT")
{genotype_list = list("WT", "WTFIX")
  }  else if (target_genotype == "DD")
  {genotype_list = list("DD", "DDFIX")
  }  else
{genotype_list = list("WT", "DD", "WTFIX", "DDFIX")}

# find all the .gnet files
filename_list = list.files(pattern = "\\.gnet$")
datarecords <- read.table("MitoTable.csv", sep=",", header=TRUE)

num_graphs = length(list.files(pattern = "\\.gnet$"))

print(c("total gnet files = ", num_graphs))

# declare vectors to hold the different result columns that will eventually be merged into a dataframe
filename_vec <- vector(mode = "character", length = num_graphs)
ncount_vec <- vector(mode = "numeric", length = num_graphs)
pcount_vec <- vector(mode = "numeric", length = num_graphs)
totalcomponents_vec <- vector(mode = "numeric", length = num_graphs)
next_biggest_vec <- vector(mode = "numeric", length = num_graphs)



graphcounter <- 0

# now iterate through all the .gnet files
for (graphindex in 1:num_graphs)
  # iterate over all .gnet files
{
  filename = filename_list[graphindex]
  
  # search through the metadata file to find the genotype of that file
  genotype = as.character(datarecords[datarecords$Name == sub(".gnet", "", filename), "Type"])
  
  if (any(genotype %in% genotype_list) )
  { 
  # load in the .gnet file removint the first line which igraph doesn't expect
  graph_data <- read.table(filename, skip="1")
  
  # convert from gnet into a graph representation
  # need to add 1 to every node label because igraph doesn't allow a label of 0
  g3 <- graph_from_edgelist(as.matrix(graph_data[,c(1,2)]) + 1, directed=FALSE)
  
  # decompose into separate components and then find the largest one
  comps <- decompose(g3)
  giant_component <- comps[[which.max(lengths(comps))]]
  
  # find the size of the giant component
  graph_size = vcount(giant_component)
  
  # get the number of degree 1 and 3 vertices
  degree1 = graph_size*degree_distribution(giant_component, cumulative = FALSE, loops = TRUE, normalized = FALSE)[2]
  
  degree3 = graph_size*degree_distribution(giant_component, cumulative = FALSE, loops = TRUE, normalized = FALSE)[4]
  
  # in the case that a missing value is returned set the count to zero
  degree1[is.na(degree1)] <- 0
  degree3[is.na(degree3)] <- 0
  
  
  
  # test for a valid mitograph
  # defined as those graphs that only have ivertices of degree 1 or 3
  countsum <- degree1 + degree3
  if (graph_size == countsum)
  {
    # store the values in the results vecors
    graphcounter <- graphcounter + 1
    filename_vec[graphcounter] = filename
    pcount_vec[graphcounter] = degree1
    ncount_vec[graphcounter] = degree3
  } # testing for valid mitograph
  

  } # test if genotype is in the list
  
} # iterate over all .gnet files


print(c("total gnet files = ", num_graphs))
print(c("total valid mitographs = ", graphcounter))

# combine all the results into one dataframe 
results_frame = data.frame(File = filename_vec, P = pcount_vec, N = ncount_vec)

# remove entries that did not get values stored because they were invalid mitographs
corrected_results_frame <- results_frame %>% filter(P != 0 & N != 0)

# write out the corrected dataframe
write.csv(corrected_results_frame, "results.csv", row.names = FALSE)

# now make a plot
max_n = max(ncount_vec)
max_p = max(pcount_vec)

# make a complete grid that will contain all the possible points plus impossible ones
# make a dataframe storing all the grid locations
# have a column that is the count
# go through and increment the count everytime a [p,n] combination comes up
# and then at the end eliminate any rows that have a count of zero.
pvals = seq(from = 0, to = max_p, by = 1)
nvals = seq(from = 0, to = max_n, by = 1)


# make a dataframe to hold a grid of n and p values for plotting
grid_data <- expand.grid(n = 0:(max_n), p = 0:(max_p))

# append a column for values that will store the counts
grid_data$value <- rep(0, times = (max_n + 1)*(max_p + 1))

for (validindex in 1:graphcounter)
  # iterate over all valid graphs in the output dataframe
{
  current_row = results_frame[validindex,,]
  p_coord = as.numeric(current_row[2]) 
  n_coord = as.numeric(current_row[3])  
  current_value = as.numeric(grid_data[grid_data$n == n_coord & grid_data$p == p_coord, "value"])
  update_value <-current_value + 1
  grid_data$value[grid_data$n == n_coord & grid_data$p == p_coord] <- update_value
}



# erase any entries that never occur in the dataset 
# and also require any value to occur more than once
corrected_grid_data <- grid_data %>% filter(value > 1)

# expand the scale on the x and y axis so that the markers are separated
corrected_grid_data$n <- corrected_grid_data$n 
corrected_grid_data$p <- corrected_grid_data$p 

plothandle <- ggplot(corrected_grid_data, aes(x = p, y = n, fill = value)) +
  # shape = 21 allows us to use an outline (color) and an inner fill (colormap)
  geom_point(shape = 21, size = 2, color = "black", stroke = 1) + 
  # Apply a continuous colormap (e.g., viridis)
  scale_fill_viridis_c(option = "plasma") + 
  # Fix aspect ratio to make sure markers are perfectly circular
  coord_fixed(ratio = 1) + 
  theme(legend.position = "right") +
  xlim(0,20) +
  ylim(0,30) +
  theme(
    # Inner plot area background
    panel.background = element_rect(fill = "white", color = "black"),
    
    # Outer plot area background (margins)
    plot.background = element_rect(fill = "white"),
    
    # Major grid lines (aligned with axis ticks)
    panel.grid.major = element_line(color = "black")
    
    
  )

print(plothandle)

print("done")



