## Turning membership in LA-ICP-MS compositional groups into networks of economic relationships

#'  The basic idea here is that, leveraging the criterion of abundance (Bishop et al., 1982), 
#'  as similarities in membership in different compositional groups increases between archaeological
#'  communities, so does the likelihood that individuals from those communities are engaging in 
#'  direct economic interactions. As used here, economic interactions are built around the concept of
#'  weak ties (Granovetter 1973). In contrast to ties that are built on deep affinity such as close
#'  friendships, family or marriage relationships, weak ties might be acquaintances or a stranger
#'  with a common cultural background. Weak ties emanating from economic relationships related to
#'  ceramic industry are constructed through such behaviors as exchange relationships, overlapping 
#'  resource exploitation areas, or similar production regimes.  
#'  
#'  Using the Brainerd Robinson coefficient of similarity, it is possible to create networks of 
#'  economic relationships through community-based membership in compositional groups. That is, 
#'  beginning with raw elemental data produced from LA-ICP-MS of 543 ceramic artifacts, several 
#'  compositional groups were identified in the Late Prehistoric central Illinois River valley. 
#'  The Brainerd Robinson coefficient of similarity assesses how similar any two given sites 
#'  are based on similarities in the number of individual sherd assignments in different 
#'  compositional groups. This method provides means to model relational economic interactions
#'  in an archaeological region. 

# Statistically robust compositional groups were identified in `Ceramic Analysis.R`. Beginning with 
# those counts, we'll clean the data and apply the Brainerd Robinson coefficient of similarity. 
# Networks are then constructed and analyzed to reveal changing patterns of economic relationships 
# in Illinois' archaeological heritage. 

# First, we'll load in some package libraries
library(tidyverse)
library(igraph)
library(corrplot)
library(reshape2)

# Then read in compositional group count data by site
comp_group <- read_csv("group assignments by site.csv")

# Two of the eight compositional groups are based on equivocal membership probabilities
# While a core group was extracted and refined, we need to drop the sherds that were 
# unable to be assigned to a core sub-group as well as those that were not able to be
# assigned to any other group. 
comp_group_refined <- comp_group %>%
                        select(-`Core A`, -unassigned)

# Sum up all of the retained sherds for compositional group construction
comp_group_refined %>%
  gather(key = Site, value = `Core A1`:`Outgroup 2`) %>%
  rename(count = "\`Core A1\`:\`Outgroup 2\`") %>%
  summarize(total_count = sum(count, na.rm = TRUE))
# Total is 314 out of the original 543, or 63% of the original data set

# Look at total number of sherds from each site
comp_group_refined %>%
  gather(key = group, value = `Core A`:`Outgroup 2`, -Site) %>%
  rename(count = "\`Core A\`:\`Outgroup 2\`") %>%
  group_by(Site) %>%
  summarize(total = sum(count, na.rm = TRUE))
# Perhaps the most problematic site here is Star Bridge, which had a massive drop from 
# ~30 or so sherds analyzed but only 9 placed within compositional groups. 
# Nevertheless, all 18 sites are represented by at least 8 sherds - not too bad. 

# The Brainerd-Robinson coefficient is a similarity metric that is unique to archaeology,
# and is used to compare assemblages based on proportions of categorical data such as
# vessel or point types. 

# The Brainerd-Robinson coefficient has been coded in R by Matt Peeples 
# (http://www.mattpeeples.net/BR.html) and by Gianmarco Alberti 
# (http://cainarchaeology.weebly.com/r-function-for-brainerd-robinson-similarity-coefficient.html). 
# Here, I follow Matt Peeple's BRsim implementation because it is substantially less resource
# intensive. However, I include a rescaling feature to rescale the BR coefficients
# from 0 - 200 to 0 - 1, which makes the output amenable for the construction
# of network graphs. 

# The input for the function is a dataframe with assemblages to be compared are found in 
# the rows and the categorical variables (such as pottery/lithic types, objects, 
# compositional groups, etc.) comprise the columns. Each variable is the numerical 
# amount of a particular categorical variable found at each site/sample/discrete 
# observation unit. 

# Here is the BRsim function as coded by Gianmarco
BRsim <- function(x, correction, rescale) {
  if(require(corrplot)){
    print("corrplot package already installed. Good!")
  } else {
    print("trying to install corrplot package...")
    install.packages("corrplot", dependencies=TRUE)
    suppressPackageStartupMessages(require(corrplot))
  }
  rd <- dim(x)[1]
  results <- matrix(0, rd, rd)
  if (correction == T){
    for (s1 in 1:rd) {
      for (s2 in 1:rd) {
        zero.categ.a <-length(which(x[s1,] == 0))
        zero.categ.b <-length(which(x[s2,] == 0))
        joint.absence <-sum(colSums(rbind(x[s1,], x[s2,])) == 0)
        if(zero.categ.a == zero.categ.b) {
          divisor.final <- 1
        } else {
          divisor.final <- max(zero.categ.a, zero.categ.b) - joint.absence+0.5
        }
        results[s1,s2] <- round((1 - (sum(abs(x[s1,] / sum(x[s1,]) - x[s2,] / sum(x[s2,]))))/2)/divisor.final,
                                digits=3)
      }
    } 
  } else {  
    for (s1 in 1:rd) {
      for (s2 in 1:rd) {
        results[s1,s2] <- round(1 - (sum(abs(x[s1,] / sum(x[s1,]) - x[s2, ] / sum(x[s2,]))))/2, digits=3)
      }
    }
  }
  rownames(results) <- rownames(x)
  colnames(results) <- rownames(x)
  col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white", "cyan", "#007FFF", "blue", "#00007F"))
  if (rescale == F) {
    upper <- 200
    results <- results * 200
  } else { 
    upper <- 1.0
  }
  corrplot(results, method="square", addCoef.col="red", is.corr=FALSE, cl.lim = c(0, upper), col = col1(100), tl.col="black", tl.cex=0.8) 
  return(results)
}

# Here is a more simplified version from Matt Peeples 
# Function for calculating Brainerd-Robinson (BR) coefficients
# *Note there is data pre-processing for Matt's script not included here
BR <- function(x) {
  rd <- dim(x)[1]
  results <- matrix(0,rd,rd)
  for (s1 in 1:rd) {
    for (s2 in 1:rd) {
      x1Temp <- as.numeric(x[s1, ])
      x2Temp <- as.numeric(x[s2, ])
      br.temp <- 0
      results[s1,s2] <- 200 - (sum(abs(x1Temp - x2Temp)))}}
  row.names(results) <- row.names(x)
  colnames(results) <- row.names(x)
  return(results)}

# My editing of the two
BR_au <- function(x, rescale = FALSE, counts = TRUE) {
  if (counts == T){
    x <- prop.table(as.matrix(x), 1) * 100
  } else {
  }
  rd <- dim(x)[1]
  results <- matrix(0,rd,rd)
  for (s1 in 1:rd) {
    for (s2 in 1:rd) {
      x1Temp <- as.numeric(x[s1, ])
      x2Temp <- as.numeric(x[s2, ])
      br.temp <- 0
      results[s1,s2] <- 200 - (sum(abs(x1Temp - x2Temp)))
    }
  }  
  row.names(results) <- row.names(x)
  colnames(results) <- row.names(x)
  if (rescale == F) {
    return(results)
  } else { 
    results <- results / 200
    return(results)
  }
}

# Before we run the BR functions, the data frame needs to have the Sites become a row name
# because the BR functions all take as inputs counts or percentages only. 
rownames(comp_group_refined) <- comp_group_refined$Site 
comp_group_refined <- comp_group_refined[, -1]

# Also need to change NAs into 0 (two methods provided below)
comp_group_refined <- comp_group_refined %>%
                        mutate_all(funs(replace(., is.na(.), 0)))

#    comp_group_refined %>%
#      mutate_all(funs(coalesce(., 0L)))

# Lost the rownames during manipulation, need to add them again
rownames(comp_group_refined) <- comp_group$Site

# A big advantage of Gianmarco's BR function is a succinct correlation plot. It can be thought
# of as a "heat-map" for BR similarities. 
BRsim(comp_group_refined, correction = FALSE, rescale = TRUE)

eco_BR <- BR_au(comp_group_refined, rescale = TRUE)

# The results of the BRsim function come in the form of an adjacency matrix. igraph 
# can easily handle this kind of data to create a network graph. Because the adjacency
# matrix is between 0 and 1, we need to tell igraph that the resulting network graph is
# weighted. Otherwise an edge will only be given for the relationship between each site 
# and itself. 
ecoBRgraph <- graph_from_adjacency_matrix(eco_BR, weighted = T)
BRel <- as_edgelist(ecoBRgraph) # convert to 2 column edgelist
BRw <- as.data.frame(E(ecoBRgraph)$weight) # extract edge weights
BRwel <- cbind(BRel, BRw) # append edge weights to edgelist
BRwel <- rename(BRwel, weight = `E(ecoBRgraph)$weight`) # rename weight column

# Assessing the distribution of the BR coefficients
BRwel %>%
  filter(`1` != `2`) %>% # drop recursive edges
  ggplot(aes(x = weight)) + 
  geom_histogram(aes(y = ..density..), bins = 25, colour = "black", fill = "white") +
  geom_density(alpha = 0.2) +
  geom_vline(aes(xintercept = mean(weight, na.rm = T)),   # Ignore NA values for mean
             color = "red", linetype = "dashed", size = 1) +
  xlab("Rescaled BR Coefficients") +
  ylab("Density") +
  theme_minimal() 

# Mean of BR coefficients (this will be used as a cutoff point for giving edges)
BRwel %>%
  filter(`1` != `2`) %>%
  filter(weight != 1) %>%
  filter(weight != 0) %>%
  summarise(Mean = mean(weight)) 
# Looks like the mean is 0.556. This can be round down to 0.55 for edge cutoffs


# But before we apply that cutoff, let's explore the range and frequency of BR
# scores if they were produced purely by chance based on our data set

# First, we will row and column randomize the BR input 10,000 times and create a list 
# of the results
# This means that we'll shuffle the order of row and column data with replacement
BReco_rand_list <- replicate(10000, comp_group_refined[sample(1:nrow(comp_group_refined), replace = T),
                                            sample(1:ncol(comp_group_refined), replace = T)], simplify = F)

# Setup an empty list to hold the BR coefficients for the randomized data
BReco_rand_result <- list()

# Number of simulations
nsim <- 10000

# Now we can iterate the BR algorithm over the randomized lists
for (i in 1:nsim) {
  BReco_rand_result[[i]] <- BR_au(BReco_rand_list[[i]], rescale = T)
}

# Turn adjacency matrices into three column data frames
for (i in 1:nsim) {
  BReco_rand_result[[i]] <- setNames(melt(BReco_rand_result[[i]]), c('1', '2', 'values'))
}

# Now we can extract the BR values from the data frames in the list
BReco_rand_result_vals <- lapply(BReco_rand_result, '[[', 3)

# And collapse that list into one long vector and turn into a tibble data frame
BReco_rand_vals <- tbl_df(unlist(BReco_rand_result_vals))

# Add a column to indicate these are simulated data
BReco_rand_vals <- BReco_rand_vals %>%
  mutate(Type = "Randomized BR")

# Append the actual data
BRwel <- tbl_df(BRwel)
BReco_vals_all <- BRwel %>%
  select(weight) %>%
  mutate(value = weight) %>%
  select(value) %>%
  mutate(Type = "Actual BR") %>%
  bind_rows(., BReco_rand_vals)

# Drop 0's and 1's since no sites are perfectly dissimilar or similar
BReco_vals_all <- BR_vals_all %>%
  filter(value != 1) %>%
  filter(value != 0)

# Plot density histograms of the observed and simulated BR coefficients
ggplot(BReco_vals_all, aes(x = value)) + 
  geom_histogram(data = subset(BReco_vals_all, Type == "Randomized BR"), aes(y=..density..),
                 alpha = 0.5, bins = 30, colour = "black", fill = "#2ca02c") +
  geom_density(data = subset(BReco_vals_all, Type == "Randomized BR"), 
               alpha = 0.1, color = "#2ca02c" , fill = "#2ca02c" , adjust = 2.5) + 
  geom_vline(data = subset(BReco_vals_all, Type == "Randomized BR"),
             aes(xintercept = mean(value, na.rm = T)),   # Ignore NA values for mean
             color = "#2ca02c" , linetype = "dashed", size = 1) +
  geom_histogram(data = subset(BReco_vals_all, Type == "Actual BR"),
                 aes(y = ..density..), bins = 30, colour = "black", fill = "#1f77b4" , alpha = 0.4) +
  geom_density(data = subset(BReco_vals_all, Type == "Actual BR"),
               alpha = 0.1, color = "#1f77b4" , fill = "#1f77b4" ) +
  geom_vline(data = subset(BReco_vals_all, Type == "Actual BR"),
             aes(xintercept = mean(value, na.rm = T)), color = "#1f77b4" , linetype = "dashed", size = 1) +
  xlab("Rescaled BR Coefficients") +
  ylab("Density") +
  theme_minimal() 

# "#1f77b4" = d3 blue
# "#2ca02c" = d3 green

# Looks like the simulated and observed data actually share similar distributions. 
# Nevertheless, there are significant nuances seen in the observed data, suggesting 
# deviations from random chance and a slightly lower than expected mean BR
# coefficient. This could reflect the small number of compositional groups (6), 
# limited number of samples from a few sites (some have 8 or 9 samples), or 
# simply a reflection of the limited geological diversity present in the CIRV. 
# However, applying the the > 0.55 cutoff indicates that edges will be 
# given in situations where the proportional similarity between two assemblages is 
# greater than the average proportional similarity across economic relationships in 
# the Late Prehistoric CIRV.

# Let's apply the 0.55 threshold
BRel_t <- BRwel %>% 
  filter(weight > 0.55 & `1` != `2`)

# Change column names to be suitable for Gephi
colnames(BRel_t) <- c("Source", "Target", "weight")

# Add columns with additional node information
# Read in tables of site names, geographic coords., and time distinction
# For time, 1 is a primary occupation prior to Oneota in-migration
# and 2 is a primary occupation succeeding Oneota in-migration
node_table <- read_csv("Jar_node_table.csv")
colnames(node_table) <- c("Source", "Label", "Long", "Lat", "Time")

# Join the node table columns to the edgelist by the Source node
econet_t1 <- left_join(BRel_t, node_table[-2], by = "Source")

# Prepare node tables to join time designation for the target node
colnames(node_table) <- c("Target", "Label", "Long", "Lat", "Time2")

# Join Time 2 column to Target node 
econet_edgelist_complete <- left_join(econet_t1, node_table[c(-2:-4)], by = "Target")

# Create Pre- and Post-Migration Edgelists
econet_pre_el_need_dist <-  econet_edgelist_complete %>%
                              filter(Time == Time2) %>%
                              filter(Time == 1)

econet_post_el_need_Law <- econet_edgelist_complete %>%
                              filter(Time == Time2) %>%
                              filter(Time == 2)

# Two sites have extended or multi-component occupations in both time periods
# So we need to include their connections in both time periods
Law_econet_post <- econet_edgelist_complete %>%
                      filter(Time == 2 & Target == "Lawrenz Gun Club" | 
                               Source == "Lawrenz Gun Club" & Time2 == 2) %>%
                      mutate(Time = replace(Time, Time == 1, 2)) %>%
                      mutate(Time2 = replace(Time2, Time2 == 1, 2))

Buck_econet_post <- econet_edgelist_complete %>%
                      filter(Time == 2 & Target == "Buckeye Bend" | 
                               Source == "Buckeye Bend" & Time2 == 2) %>%
                      mutate(Time = replace(Time, Time == 1, 2)) %>%
                      mutate(Time2 = replace(Time2, Time2 == 1, 2))

# Bind the LCG & Buckeye post-migration edges to the post-migration edgelists
econet_post_el_need_dist <- rbind(econet_post_el_need_Law, Law_econet_post, Buck_econet_post)

# Adding geographic coordinates
# Read in matrix of site distances
site_distances <- read_csv("Site Distances Matrix in km.csv")

#first column of site names to rownames 
site_distances <- column_to_rownames(site_distances, var = "X1")

# Convert geographic distance matrix to graph object
distance_g <- graph_from_adjacency_matrix(as.matrix(site_distances), weighted = TRUE, 
                                          mode = "directed")

# Convert geo distance graph object to edgelist
distance_el <- as_edgelist(distance_g)
distance_el_weight <- as.numeric(E(distance_g)$weight)
distance_el <- tbl_df(cbind(distance_el, distance_el_weight))
colnames(distance_el) <- c("Source", "Target", "weight")
distance_el$Distance <- as.numeric(distance_el$weight)

# Merge the geographic distance edgelist with directed plate edgelists
econet_pre_el_complete <- merge(econet_pre_el_need_dist, distance_el[-3])
econet_post_el_complete <- merge(econet_post_el_need_dist, distance_el[-3])

# Combine the pre- and post-migration data sets into a single edgelist
econet_el_BR_all_time_complete <- rbind(econet_pre_el_complete, econet_post_el_complete)

# Finally, we can export the complete edgelist for visualization in Gephi
write_csv(econet_el_BR_all_time_complete, "Economic_network_BR_edgelist_complete_.csv")

#### Undirected Economic Networks ####



# The edgelists created thus far have been directed. Since we are disregarding 
# directionality, it is imporant to account for duplicate edges. 
BRgraph_un <- graph_from_adjacency_matrix(eco_BR, weighted = T, mode = "undirected")

# Create undirected edgelist
BRel_un <- as_edgelist(BRgraph_un)

# Create the weights and format as a data frame for column binding
BRw_un <- E(BRgraph_un)$weight
BRw_un <- as.data.frame(BRw_un)

# Add the weights, and viola we have a weighted, directed edgelist for proportional
# stylistic similarity between sites. 
BRel_un <- cbind(BRel_un, BRw_un)

# Write out complete Brainerd Robinson edgelist
write_csv(BRel_un, "complete_ECO_BR_UNDIRECTED_edgelist.csv")

# Apply our threshold of > 0.4 so that we only give UNDIRECTED edges to the strongest
# proportional relationship. We can use dplyr to wrangle the edgelist and also drop
# recursive edges. 
BRel_t_un <- BRel_un %>% 
  filter(BRw_un > 0.55 & BRel_un[1] != BRel_un[2])

# Change column names to be suitable for Gephi
colnames(BRel_t_un) <- c("Source", "Target", "weight")

colnames(node_table) <- c("Source", "Label", "Long", "Lat", "Time")

# Join the node table columns to the edgelist by the Source node
eco_t1_un <- left_join(BRel_t_un, node_table[-2], by = "Source")

# Prepare node tables to join time designation for the target node
colnames(node_table) <- c("Target", "Label", "Long", "Lat", "Time2")

# Join Time 2 column to Target node 
econet_edgelist_complete_un <- left_join(eco_t1_un, node_table[c(-2:-4)], by = "Target")

# Create Pre- and Post-Migration Edgelists
econet_pre_el_need_dist_un <-  econet_edgelist_complete_un %>%
                                  filter(Time == Time2) %>%
                                  filter(Time == 1)

econet_post_el_need_Law_un <- econet_edgelist_complete_un %>%
                                filter(Time == Time2) %>%
                                filter(Time == 2)

# Two sites have extended or multi-component occupations in both time periods
# So we need to include their connections in both time periods
Law_econet_post_un <- econet_edgelist_complete_un %>%
                        filter(Time == 2 & Target == "Lawrenz Gun Club" | 
                                 Source == "Lawrenz Gun Club" & Time2 == 2) %>%
                        mutate(Time = replace(Time, Time == 1, 2)) %>%
                        mutate(Time2 = replace(Time2, Time2 == 1, 2))

Buck_econet_post_un <- econet_edgelist_complete_un %>%
                        filter(Time == 2 & Target == "Buckeye Bend" | 
                                 Source == "Buckeye Bend" & Time2 == 2) %>%
                        mutate(Time = replace(Time, Time == 1, 2)) %>%
                        mutate(Time2 = replace(Time2, Time2 == 1, 2))

# Bind the LCG & Buckeye post-migration edges to the post-migration edgelists
econet_post_el_need_dist_un <- rbind(econet_post_el_need_Law_un, 
                                    Law_econet_post_un, Buck_econet_post_un)

# Merge the geographic distance edgelist with undirected plate edgelists
econet_pre_el_complete_un <- merge(econet_pre_el_need_dist_un, distance_el[-3])
econet_post_el_complete_un <- merge(econet_post_el_need_dist_un, distance_el[-3])

# Combine the pre- and post-migration data sets into a single edgelist
econet_el_BR_all_time_complete_un <- rbind(econet_pre_el_complete_un, econet_post_el_complete_un)

# Finally, we can export the complete undirected edgelist for visualization in Gephi
write_csv(econet_el_BR_all_time_complete_un, "Econet_BR_UNDIRECTED_edgelist_complete_.csv")
write_csv(econet_pre_el_complete_un, "Econet_BR_UNDIRECTED_edgelist_pre-migration_.csv")
write_csv(econet_post_el_complete_un, "Econet_BR_UNDIRECTED_edgelist_post-migration_.csv")


