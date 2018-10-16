# Geochemical compositional group economic network statistics

library(tidyverse)
library(readxl)
library(broom)
library(igraph)
library(cowplot)


#----------------------------Economic BR Network Stats----------------------------####

# Read in finalized, undirected plate BR edgelist 
BReco_el_un <- read_csv("Econet_BR_UNDIRECTED_edgelist_complete_.csv")

# Read in finalized, undirected pre-migration BR edgelist
BReco_el_un_pre <- read_csv("Econet_BR_UNDIRECTED_edgelist_pre-migration_.csv")

# Read in finalized, undirected post-migration BR edgelist
BReco_el_un_post <- read_csv("Econet_BR_UNDIRECTED_edgelist_post-migration_.csv")

# Convert to igraph graph
BReco_g <- graph_from_edgelist(as.matrix(BReco_el_un[, c(1:2)]), directed = FALSE)
BReco_g_pre <- graph_from_edgelist(as.matrix(BReco_el_un_pre[, c(1:2)]), directed = FALSE)
BReco_g_post <- graph_from_edgelist(as.matrix(BReco_el_un_post[, c(1:2)]), directed = FALSE)

# Assign edge weights to graph
E(BReco_g)$weight <- BReco_el_un$weight
E(BReco_g_pre)$weight <- BReco_el_un_pre$weight
E(BReco_g_post)$weight <- BReco_el_un_post$weight

# Function to calculate degree, betweenness, closeness, and eigenvector centrality 
# for a graphand return a data frame with the scores
centr_all <- function(graph, g_name = "Score") {
  
  # Check that graph is an igraph object
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  
  # Name of graph
  g_name <- as.character(g_name)
  
  # Degree centralization
  res_centr <- centr_degree(graph)$centralization
  
  # Betweenness centralization
  res_centr[2] <- centr_betw(graph)$centralization
  
  # Closeness centralization
  res_centr[3] <- centr_clo(graph)$centralization
  
  # Eigenvector centralization
  res_centr[4] <- centr_eigen(graph)$centralization
  
  res_centr <- t(as.data.frame(res_centr))
  
  # Table of scores
  colnames(res_centr) <- c("Degree", "Betweenness", "Closeness", "Eigenvector")
  rownames(res_centr) <- g_name
  
  res_centr
}

# Calculate centralization scores for each graph
all_centr <- centr_all(BReco_g, g_name = "Flattened Across Time")
pre_centr <- centr_all(BReco_g_pre,  g_name = "Pre-Migration")
post_centr <- centr_all(BReco_g_post,  g_name = "Post-Migration")
rbind(pre_centr, post_centr, all_centr)

# Calculated Mean Weighted Degree (or strength)
mean(strength(BReco_g))
mean(strength(BReco_g_pre))
mean(strength(BReco_g_post))

#---------------------------Edge Betweenness Community Detection-------------------####

# Edge betweenness extends the concept of vertex betweenness centrality to edges by 
# assigning each edge a score that reflects the number of shortest paths that move 
# through that edge. 
# You might ask the question, which ties in a social network are the most important in 
# the spread of information?

# Calculated edge betweenness score for each network
ecopre_eb <- cluster_edge_betweenness(BReco_g_pre)
ecopost_eb <- cluster_edge_betweenness(BReco_g_post)
ecoall_eb <- cluster_edge_betweenness(BReco_g)

# Edge betweenness correctly assigns the pre- and post-migration
# sites to clusters, but with some interesting intricacies - Buckeye in pre and
# Lawrenz in post

# The pre- and post-migration eb communities are interesting as well

# Community detection via edge betweenness plot_across time
plot(ecoall_eb, BReco_g, col = membership(ecoall_eb), vertex.label.cex = c(1), 
     edge.arrow.size = .1, edge.curved = .1)
title(main = "Edge Betweenness Community Detection in \n the Economic Network", 
      cex.main = 1.5)

# Community detection via edge betweenness plot_pre-migration
plot(ecopre_eb, BReco_g_pre, col = membership(ecopre_eb), vertex.label.cex = c(1), 
     edge.arrow.size = .1, edge.curved = .1)
title(main = "Edge Betweenness Community Detection in \n the Economic Network", 
      cex.main = 1.5)

# Community detection via edge betweenness plot_post-migration
plot(ecopost_eb, BReco_g_post, col = membership(ecopost_eb), vertex.label.cex = c(1), 
     edge.arrow.size = .1, edge.curved = .1)
title(main = "Edge Betweenness Community Detection in \n the Economic Network", 
      cex.main = 1.5)

#-----------------Pre Randonmization for Pre-Migration Period Economic BR--------------####
#----------------------------------PRE_MIGRATION------------------------------------------#

# Initiate empty list for assessing BR pre-migration average path length and transitivity
gecopre <- vector('list', 5000)

# Initiate empty list for assessing BR pre-migration density density and mean degree
gecopre.d <- vector('list', 5000)

# Populate gpre list with random graphs of same order and size
for(i in 1:5000){
  gecopre[[i]] <- erdos.renyi.game(n = gorder(BReco_g_pre), p.or.m = gsize(BReco_g_pre),
                                   directed = FALSE, type = "gnm")
}

# Populate gecopre.d list with random graphs of same order and approximate density. A separate list 
# of 5000 randon graphs is necessary for density and mean degree because these statistics would 
# identical in random graphs of the same order and size as our observed graph. 
# Instead, a probability of edge creation equal to the observed density is used. Further, 
# only mean degree (as opposed to mean weighted degree) is used because Erdos-Renyi 
# random graphs do not support weights. 
for(i in 1:5000){
  gecopre.d[[i]] <- erdos.renyi.game(n = gorder(BReco_g_pre), p.or.m = edge_density(BReco_g_pre), 
                                     directed = FALSE, type = "gnp")
}

# Calculate average path length, transitivity (custering coefficient), density, and degree across 
# the 5000 random pre-migration graphs
ecopre.pl <- lapply(gecopre.d, mean_distance, directed = FALSE)
ecopre.trans <- lapply(gecopre, transitivity)
ecopre.density <- lapply(gecopre.d, edge_density)
ecopre.degree <- lapply(gecopre.d, function(x){
  y <- degree(x)
  mean(y)
}
)

# Unlist and change to a data frame for vizualizations
ecopre.pl <- as.data.frame(unlist(ecopre.pl))
ecopre.trans <- as.data.frame(unlist(ecopre.trans))
ecopre.density <- as.data.frame(unlist(ecopre.density))
ecopre.degree <- as.data.frame(unlist(ecopre.degree))

# Plot the distribution of random graph's average shortest path lengths with the pre-migration 
# BR network's ave. shortest path as line
p.ecopre.pl <- ggplot(ecopre.pl, aes(x = unlist(ecopre.pl))) + 
  geom_histogram(aes(y = ..density..), bins = 24) + 
  geom_vline(xintercept = (mean_distance(BReco_g_pre, directed = FALSE)), 
             linetype = "dashed", color = "red") +
  geom_density() +
  ggtitle("Distribution of 5000 Random Graph Average Shortest Path Lengths & \nPre-Migration Period Average Shortest Path Length") + 
  xlab("Average Shortest Path Length") +
  ylab("")

# Plot the distribution of random graph's transitivity with the pre-migration BR network's 
# transitivity path as line
p.ecopre.trans <- ggplot(ecopre.trans, aes(x = unlist(ecopre.trans))) + 
  geom_histogram(aes(y = ..density..), bins = 22) + 
  geom_vline(xintercept = (transitivity(BReco_g_pre)), linetype = "dashed", color = "red") +
  geom_density() +
  ggtitle("Distribution of Transitivity in 5000 Random Models & \nPre-Migration Period Network Transitivity") + 
  xlab("Transitivity (or Clustering Coefficient)") +
  ylab("")

# Plot the distribution of random graph's average density with the pre-migration jar network's
# ave. shortest path as line
p.ecopre.density <- ggplot(ecopre.density, aes(x = unlist(ecopre.density))) + 
  geom_histogram(aes(y = ..density..), bins = 22) + 
  geom_vline(xintercept = (edge_density(BReco_g_pre)), linetype = "dashed", color = "red") +
  geom_density() +
  ggtitle("Distribution of 5000 Random Graph Average Densities &\nPre-Migration Preiod Network Average Density") + 
  xlab("Average Density") +
  ylab("")

# Plot the distribution of random graph's mean degree with the pre-migration BR network's mean
# degree path as line
p.ecopre.degree <- ggplot(ecopre.degree, aes(x = unlist(ecopre.degree))) + 
  geom_histogram(aes(y = ..density..), bins = 22) + 
  geom_vline(xintercept = (mean(degree(BReco_g_pre, mode = "all"))), 
             linetype = "dashed", color = "red") +
  geom_density() +
  ggtitle("Distribution of Mean Degree in 5000 Random Models & \nPre-Migration Period Network Mean Degree") + 
  xlab("Mean Degree") +
  ylab("")

# Use plot_grid to plot all four graphs on the same grid
plot_grid(p.ecopre.pl, p.ecopre.trans, p.ecopre.density, p.ecopre.degree)

# Calculate the proportion of graphs with an average path length lower than observed
sum(ecopre.pl < mean_distance(BReco_g_pre, directed = False))/5000*100

# Calculate the proportion of graphs with a transitivity (mean clustering coefficient) 
# lower than our observed
sum(ecopre.trans < transitivity(BReco_g_pre))/5000*100

# Calculate the proportion of graphs with a density lower than our observed
sum(ecopre.density < edge_density(BReco_g_pre))/5000*100

# Calculate the proportion of graphs with a mean degree lower than observed
sum(ecopre.degree < mean(degree(BReco_g_pre)))/5000*100

#-----------------Post Randonmization for Post-Migration Period Economic BR--------------####
#----------------------------------POST_MIGRATION------------------------------------------#

# Initiate empty list for assessing BR pre-migration average path length and transitivity
gecopost <- vector('list', 5000)

# Initiate empty list for assessing BR pre-migration density density and mean degree
gecopost.d <- vector('list', 5000)

# Populate gpre list with random graphs of same order and size
for(i in 1:5000){
  gecopost[[i]] <- erdos.renyi.game(n = gorder(BReco_g_post), p.or.m = gsize(BReco_g_post),
                                   directed = FALSE, type = "gnm")
}

# Populate gecopre.d list with random graphs of same order and approximate density. A separate list 
# of 5000 randon graphs is necessary for density and mean degree because these statistics would 
# identical in random graphs of the same order and size as our observed graph. 
# Instead, a probability of edge creation equal to the observed density is used. Further, 
# only mean degree (as opposed to mean weighted degree) is used because Erdos-Renyi 
# random graphs do not support weights. 
for(i in 1:5000){
  gecopost.d[[i]] <- erdos.renyi.game(n = gorder(BReco_g_post), p.or.m = edge_density(BReco_g_post), 
                                     directed = FALSE, type = "gnp")
}

# Calculate average path length, transitivity (custering coefficient), density, and degree across 
# the 5000 random pre-migration graphs
ecopost.pl <- lapply(gecopost.d, mean_distance, directed = FALSE)
ecopost.trans <- lapply(gecopost, transitivity)
ecopost.density <- lapply(gecopost.d, edge_density)
ecopost.degree <- lapply(gecopost.d, function(x){
  y <- degree(x)
  mean(y)
}
)

# Unlist and change to a data frame for vizualizations
ecopost.pl <- as.data.frame(unlist(ecopost.pl))
ecopost.trans <- as.data.frame(unlist(ecopost.trans))
ecopost.density <- as.data.frame(unlist(ecopost.density))
ecopost.degree <- as.data.frame(unlist(ecopost.degree))

# Plot the distribution of random graph's average shortest path lengths with the pre-migration 
# BR network's ave. shortest path as line
p.ecopost.pl <- ggplot(ecopost.pl, aes(x = unlist(ecopost.pl))) + 
  geom_histogram(aes(y = ..density..), bins = 24) + 
  geom_vline(xintercept = (mean_distance(BReco_g_post, directed = FALSE)), 
             linetype = "dashed", color = "red") +
  geom_density() +
  ggtitle("Distribution of 5000 Random Graph Average Shortest Path Lengths & \nPre-Migration Period Average Shortest Path Length") + 
  xlab("Average Shortest Path Length") +
  ylab("")

# Plot the distribution of random graph's transitivity with the pre-migration BR network's 
# transitivity path as line
p.ecopost.trans <- ggplot(ecopost.trans, aes(x = unlist(ecopost.trans))) + 
  geom_histogram(aes(y = ..density..), bins = 10) + 
  geom_vline(xintercept = (transitivity(BReco_g_post)), linetype = "dashed", color = "red") +
  geom_density() +
  ggtitle("Distribution of Transitivity in 5000 Random Models & \nPre-Migration Period Network Transitivity") + 
  xlab("Transitivity (or Clustering Coefficient)") +
  ylab("")

# Plot the distribution of random graph's average density with the pre-migration jar network's
# ave. shortest path as line
p.ecopost.density <- ggplot(ecopost.density, aes(x = unlist(ecopost.density))) + 
  geom_histogram(aes(y = ..density..), bins = 19) + 
  geom_vline(xintercept = (edge_density(BReco_g_post)), linetype = "dashed", color = "red") +
  geom_density() +
  ggtitle("Distribution of 5000 Random Graph Average Densities &\nPre-Migration Preiod Network Average Density") + 
  xlab("Average Density") +
  ylab("")

# Plot the distribution of random graph's mean degree with the pre-migration BR network's mean
# degree path as line
p.ecopost.degree <- ggplot(ecopost.degree, aes(x = unlist(ecopost.degree))) + 
  geom_histogram(aes(y = ..density..), bins = 19) + 
  geom_vline(xintercept = (mean(degree(BReco_g_post, mode = "all"))), 
             linetype = "dashed", color = "red") +
  geom_density() +
  ggtitle("Distribution of Mean Degree in 5000 Random Models & \nPre-Migration Period Network Mean Degree") + 
  xlab("Mean Degree") +
  ylab("")

# Use plot_grid to plot all four graphs on the same grid
plot_grid(p.ecopost.pl, p.ecopost.trans, p.ecopost.density, p.ecopost.degree)

# Calculate the proportion of graphs with an average path length lower than observed
sum(ecopost.pl < mean_distance(BReco_g_post, directed = False))/5000*100

# Calculate the proportion of graphs with a transitivity (mean clustering coefficient) 
# lower than our observed
sum(ecopost.trans < transitivity(BReco_g_post))/5000*100

# Calculate the proportion of graphs with a density lower than our observed
sum(ecopost.density < edge_density(BReco_g_post))/5000*100

# Calculate the proportion of graphs with a mean degree lower than observed
sum(ecopost.degree < mean(degree(BReco_g_post)))/5000*100