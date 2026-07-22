# script to generate random graphs using the configuration model for a range of [p,n] values, and
# calculate the proportion of them that are planar

# the idea is to see whether any [p,n] values might become more or less favorable then the graphs become p lanar


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

#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("RBGL")
#library(RBGL)



#setwd("/Users/wallacemarshall/papers/Mitochondria_graph_theory/data")

# before running load the data in like this:
# datarecords <- read.table("MitoTable.csv", sep=",", header=TRUE)
# and then tell it what genotype to select like this:
# target_genotype <- "DD"
# to run the script just click the Source button!

# define the range of [p,n] to do the calculation for
max_p <- 15
max_n <- 25
number_graphs <- 50  # how many graphs to generate for each [p,n]

# fill an array to store the values
# at the end of value of 0 means an impossible combination of p and n.
result_array <- matrix(0, nrow = max_p+1, ncol = max_n+1)
result_array_planar <- matrix(0, nrow = max_p+1, ncol = max_n+1)


# first make a complete grid that will contain all the possible points plus impossible ones
# make a dataframe storing all the grid locations
# have a column that is the total count
# and another column that is the planar count
# go through and increment the count everytime a [p,n] combination comes up
# and then at the end eliminate any rows that have a count of zero.
pvals = seq(from = 0, to = max_p, by = 1)
nvals = seq(from = 0, to = max_n, by = 1)

# make a dataframe to hold a grid of n and p values for plotting
grid_data <- expand.grid(n = 0:(max_n), p = 0:(max_p))

# append a column for values that will store the plot avlues
# which are the fraction of random graphs that were planar
# -1 is a flag to indicate grid values that never occur
# since a value of zero could mean a legitimate grid point p,n but for which there are no planar graphs found
grid_data$value <- rep(-1, times = (max_n + 1)*(max_p + 1))










# iterate over the range of p and n values
# for each value, generate random mitographs and check how many are planar
for (nval in 0:max_n) 
  {
  print(nval)
  # generate valid p and n combination
  if (nval %% 2 == 0)
  {
    if (nval == 0)
    { p_min <- 2}
    else
    {
      p_min <- 0
    }
  }
  else
  {
    p_min <- 1
  }   
  p_upper_bound = min(max_p, (2 + nval))
  for (pval in seq(from = p_min, to = p_upper_bound, by = 2))
    {
    # set up degree vector for generating random graphs inside the loop
    nvector = rep(3, times=nval)
    pvector = rep(1, times=pval)
    degree_vector = c(nvector, pvector)
    cat('p=', pval, '  n=', nval, 'degree vector = ',degree_vector, '\n')
    
    graphcounter <- 0
    while (graphcounter < number_graphs)
      {
        # generate random graphs using configuration model
        # use method simple which allows self-loops and multiple edges
        rand_graph <- sample_degseq(out.deg = degree_vector, method = "simple")
        
        # make sure it just has one component
        total_graph_size = vcount(rand_graph)
        comps <- decompose(rand_graph)
        giant_component <- comps[[which.max(lengths(comps))]]
        largest_component_size = vcount(giant_component)
        
        if (largest_component_size == total_graph_size)
        {
          graphcounter <- graphcounter + 1
          # increment graph count in array
          result_array[pval+1,nval+1] = result_array[pval+1,nval+1] + 1
          
          
          # test for planarity
          # one limitation is that we need to use the graph_nel data format 
          # in order to use the test provided by RBGL
          # but graph_nel does not allow multiple edges.
          # fortunatey, multiple edges won't affect planarity so we can delete them
          # also delete self-loops since they also don't affect planarity
          rand_graph_simple <- simplify(rand_graph, remove.multiple = TRUE, remove.loops = TRUE)
          g_nel <- as_graphnel(rand_graph_simple) # requires RBGL
          is_planar <- boyerMyrvoldPlanarityTest(g_nel) # requires RBGL
          if (is_planar)
          {
            result_array_planar[pval+1,nval+1] = result_array_planar[pval+1,nval+1] + 1
            
          }
          #cat('generated a graph    value = ', result_array[pval+1,nval+1], '\n')
        } # test if a single component
      } # while loop to generate different random graphs
    
    # store the fraction of the random graphs that were planar which is what we will plot
    
    current_value = result_array[pval+1,nval+1]
    current_valueplanar = result_array_planar[pval+1,nval+1]
    cat('total = ', current_value, '   planar = ', current_valueplanar, '\n')
    value_fraction <- current_valueplanar/current_value
    grid_data$value[grid_data$n == nval & grid_data$p == pval] <- value_fraction
    
    
    
    } # loop over nval
  } # loop over pval




# now the next step is to generate a plot
# erase any entries that never occur 
corrected_grid_data <- grid_data %>% filter(value > -1)

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















