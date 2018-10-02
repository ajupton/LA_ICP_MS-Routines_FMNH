##' Analysis of ceramic LA-ICP-MS data


library(tidyverse)
library(infer)
library(broom)
library(stringr)
library(plotly)
library(rebus)
library(xlsx)
library(readxl)
library(plotly)
library(ggpubr)
library(cluster)
library(dendextend)
library(factoextra)
library(stats)
library(ICSNP)
library(shiny)
library(shinydashboard)
library(ggsci)


##### Data Import and Cleaning #####
samps <- read_csv("Upton_results_samples_shell_corrected_August_21_2018.csv")

# Remove samples that hit shell to the point of being unusable or four samples
# that were victim to a chamber leakage issue on 2017-Oct-06
removes <- str_detect(tolower(samps$Sample), "remove")
samples <- samps[!removes,]
rm_samps <- samps[removes,]

# Validate removed samples
rm_samps[,c(1:2)]

# Pull out clay samples (we'll add them back in later on)
clay_rows <- str_detect(tolower(samples$Sample), "c[:digit:][:digit:]")
clay_samps <- samples[clay_rows,]
samples <- samples[!clay_rows,]

# Clean up clay sample names
clay_samps$Sample <- str_replace(clay_samps$Sample, pattern = "_1" %R% END, "")

# Add clay i.d.'s to a separate column
clay_samps <- clay_samps %>%
                mutate(id = parse_number(clay$Sample)) %>%
                arrange(Sample)

##### Add features to ceramic samples #####
# Extract sample unique sherd i.d. number 
# First remove the run information from sample names
samples$Sample <- str_replace(samples$Sample, pattern = "_1" %R% END, "")
samples$Sample <- str_replace_all(samples$Sample, c("_run[:digit:][:digit:]" = "", 
                                                    "_run[:digit:]" = "", 
                                                    "__" = "", 
                                                    "_run [:digit:]" = "", 
                                                    "run1" = "",
                                                    "_r" = "", 
                                                    " run 1" = "", 
                                                    "_" = " "))

# Now extract the sample i.d.'s to a separate column
samples <- samples %>%
            mutate(id = parse_number(samples$Sample)) 

# Read in contextual data for the ceramic samples
ceramic_features_by_id <- read_xlsx(path = "ceramic features.xlsx", sheet = 1)
ceramic_features_by_site <- read_xlsx(path = "ceramic features.xlsx", sheet = 2)

# Join ceramic features to sample data
samples <- left_join(samples, ceramic_features_by_id)
samples <- left_join(samples, ceramic_features_by_site)

# Number of sherds by site and by vessel type
samples %>%
  group_by(Site, Vessel_Class) %>%
  summarise(num = n())  # %>%
  # write_csv("Number of sherds by site and by vessel type.csv")

# Number of sherds by site and by cultural group
samples %>%
  group_by(Site, Cultural_Group) %>%
  summarize(num = n())

# Check for any linear relationships between Calcium and the other elements
# Looks like there are some significant at a 0.05 alpha, but there is a significant
# amount of heteroscedasticity and residual variation in all but Sr, which
# expectedly does highly correlate with Ca
summary(lm(Ca ~ ., data = samples[,3:length(samples)]))

# Plotting to show how strong the linear relationships are for some elements
p <- ggplot(samples, aes(x = Sr, y = Ca)) + geom_smooth() + geom_point()
#ggplotly(p)

# Remove problem elements. Some elements are known to be unreliably measured using the ICP-MS
# at the EAF. Following Golitko (2010), these include the following elements. 
problem_elements <- c("P", "Sr", "Ba", "Ca", "Hf", "As", "Cl")

# Other elements such as Ca and Sr are affected by shell tempering. Want to drop those as well. 

# Overall these are the Elements retained - 44 in all.
elems_retained <- c("Al","B", "Be", "Ce", "Co", "Cr", "Cs", "Dy", "Er", "Eu", "FeO",
                    "Gd", "Ho", "In", "K", "La", "Li", "Lu", "Mg", "MnO", "Mo", "Na", "Nb", 
                    "Nd", "Ni", "Pb", "Pr", "Rb", "Sc", "Si", "Sm", "Sn", "Ta", "Tb", "Th", "Ti",
                    "Tm", "U", "V", "W", "Y", "Yb", "Zn", "Zr")

ceramic.names.use <- names(samples)[(names(samples) %in% elems_retained)]
#length(ceramic.names.use) == length(elems_retained) # check that all elements are retained 
samples_good <- samples %>% select(ceramic.names.use)

# Check to ensure the elements were removed are supposed to be removed
anti_join(data.frame(names(samples)), 
          data.frame(names(samples_good)), by = c("names.samples." = "names.samples_good."))

# Need to drop the "O" for oxide after elements measured as %oxide composition since they have 
# already been converted to ppm
names(samples_good) <- str_remove(names(samples_good), "O")

# Bind sample id and other data to the logged chemical concentrations 
sample_pcaready <- bind_cols(samples[,c(1:2, 71:80)], samples_good)

# Some ceramic samples were run on an older ICP-MS machine during an initial pilot study. 
# I need to tease these out pending quality control from Laure Dussubieux, a chemist at the 
# Field Museum. 
sample_old_machine <- sample_pcaready %>% filter(Date < 2016)
sample_new_pcaready <- sample_pcaready %>% filter(Date > 2016)

######## End of data cleaning, beginnging of statistical analysis ########

# First step is to take the log base 10 of all samples to account for scalar differences
# in the magnitude of chemical compositions across the elements, from major to minor to trace
sample_new_pcaready[,13:56] <- log10(sample_new_pcaready[,13:56])

##### PCA #####

# Exploring PCA
sample_pca <- sample_new_pcaready %>% 
                nest() %>% 
                mutate(pca = map(data, ~ prcomp(.x %>% select(Si:Th))),
                       pca_aug = map2(pca, data, ~augment(.x, data = .y)))

# Check variance explained by each model
var_exp_sample <- sample_pca %>% 
  unnest(pca_aug) %>% 
  summarize_at(.vars = vars(contains("PC")), .funs = funs(var)) %>% 
  gather(key = pc, value = variance) %>% 
  mutate(var_exp = variance/sum(variance),
         cum_var_exp = cumsum(var_exp),
         pc = str_replace(pc, ".fitted", ""))

# Looks like we need to retain the first 12 PC's to hit 90% of the data's variability
# Graphing this out might help
var_exp_sample %>% 
  rename(`Variance Explained` = var_exp,
         `Cumulative Variance Explained` = cum_var_exp) %>% 
  gather(key = key, value = value, 
         `Variance Explained`:`Cumulative Variance Explained`) %>% 
  mutate(pc = str_replace(pc, "PC", "")) %>%
  mutate(pc = as.numeric(pc)) %>%
  ggplot(aes(reorder(pc, sort(as.numeric(as.character(pc)))), value, group = key)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~key, scales = "free_y") +
  theme_bw() +
  lims(y = c(0, 1)) +
  labs(y = "Variance", x = "",
       title = "Variance explained by each principal component")

# Check number of PCs to retain to reach 90% of the variability in the original dataset
var_exp_sample %>% filter(cum_var_exp < 0.909) # Need to retain the first 12 PCs. 12 PCs is much less than 44 elements

# Plot the first two PCs with Geography_2 as group separation
geo2_pc1pc2 <-sample_pca %>%
              mutate(
                pca_graph = map2(
                  .x = pca,
                  .y = data,
                  ~ autoplot(.x, loadings = TRUE, loadings.label = TRUE,
                             #loadings.label.repel = TRUE,
                             loadings.label.colour = "black",
                             loadings.colour = "gray85",
                             loadings.label.alpha = 0.5,
                             loadings.label.size = 3,
                             loadings.label.hjust = 1.1,
                             frame = TRUE,
                             frame.type = "norm",
                             data = .y, 
                             colour = "Geography_2", 
                             shape = "Geography_2",
                             frame.level = .9, 
                             frame.alpha = 0.001, 
                             size = 2) +
                    theme_bw() + 
                    #geom_text(label = .y$Sample) +
                    labs(x = "Principal Component 1",
                         y = "Principal Component 2",
                         title = "First two principal components of PCA on CIRV Ceramic dataset")
                )
              ) %>%
              pull(pca_graph)


geo2_pc1pc2[[1]] + scale_fill_manual(values = c("black","black")) + 
  scale_color_manual(values = c("black","black","black")) 

# This shows significant overlap but a general trend that follows the clay: in general there is
# less elemental enrichment in clay resources in the southern portion of the CIRV compared to the northern part
# with the north-south line of demarcation being the Spoon-Illinois River confluence (clay along the Spoon 
# is included in the north)

# Check the first two PCs with Sites as group separation
site_pc1pc2 <- sample_pca %>%
  mutate(
    pca_graph = map2(
      .x = pca,
      .y = data,
      ~ autoplot(.x, loadings = TRUE, loadings.label = TRUE,
                 #loadings.label.repel = TRUE,
                 loadings.label.colour = "black",
                 loadings.colour = "gray85",
                 loadings.label.alpha = 0.5,
                 loadings.label.size = 3,
                 loadings.label.hjust = 1.1,
                 frame = TRUE,
                 frame.type = "norm",
                 data = .y, 
                 colour = "Site", 
                 shape = "Geography_2",
                 frame.level = .9, 
                 frame.alpha = 0.001, 
                 size = 2) +
        theme_bw() + 
        #geom_text(label = .y$Sample) +
        labs(x = "Principal Component 1",
             y = "Principal Component 2",
             title = "First two principal components of PCA on CIRV Ceramic dataset")
    )
  ) %>%
  pull(pca_graph)

# Interact with the chart above
ggplotly(site_pc1pc2[[1]]) # plotly drops the stat_ellipse frames for some reason
# This is a challenge to interpret, but it doesn't seem as though there is meaningful patterning
# when considering the different sites on PC1-PC2 aside from some outliers in Walsh/Crable. 

# Check the first two PCs with Vessel Class as group separation
vessel_pc1pc2 <- sample_pca %>%
  mutate(
    pca_graph = map2(
      .x = pca,
      .y = data,
      ~ autoplot(.x, loadings = TRUE, loadings.label = TRUE,
                 #loadings.label.repel = TRUE,
                 loadings.label.colour = "black",
                 loadings.colour = "gray85",
                 loadings.label.alpha = 0.5,
                 loadings.label.size = 3,
                 loadings.label.hjust = 1.1,
                 frame = TRUE,
                 frame.type = "norm",
                 data = .y, 
                 colour = "Vessel_Class", 
                 shape = "Vessel_Class",
                 frame.level = .9, 
                 frame.alpha = 0.001, 
                 size = 2) +
        theme_bw() + 
        #geom_text(label = .y$Sample) +
        labs(x = "Principal Component 1",
             y = "Principal Component 2",
             title = "First two principal components of PCA on CIRV Ceramic dataset")
    )
  ) %>%
  pull(pca_graph)

ggplotly(vessel_pc1pc2[[1]])
# The vessel graph is interesting. At first glance, it doesn't seem as though there is much in the way
# of separation by vessel class, but there appears to be some nuances to that upon futher 
# consideration. There are some plates that are low on both PC1 and PC2 axes as well as jars
# that are significantly more enriched on PC1

# How about separation by time?
time_pc1pc2 <- sample_pca %>%
  mutate(
    pca_graph = map2(
      .x = pca,
      .y = data,
      ~ autoplot(.x, loadings = TRUE, loadings.label = TRUE,
                 #loadings.label.repel = TRUE,
                 loadings.label.colour = "black",
                 loadings.colour = "gray85",
                 loadings.label.alpha = 0.5,
                 loadings.label.size = 3,
                 loadings.label.hjust = 1.1,
                 frame = TRUE,
                 frame.type = "norm",
                 data = .y, 
                 colour = "Time", 
                 shape = "Time",
                 frame.level = .9, 
                 frame.alpha = 0.001, 
                 size = 2) +
        theme_bw() + 
        #geom_text(label = .y$Sample) +
        labs(x = "Principal Component 1",
             y = "Principal Component 2",
             title = "First two principal components of PCA on CIRV Ceramic dataset")
    )
  ) %>%
  pull(pca_graph)

ggplotly(time_pc1pc2[[1]])
# Again, this appears similar to the prior PC biplot separated by vessel class - there is
# no general trend of group separation but some interesting insights when considering 
# outliers. 

# Perhaps Oneota presence may be more revealing
oneota_pc1pc2 <- sample_pca %>%
  mutate(
    pca_graph = map2(
      .x = pca,
      .y = data,
      ~ autoplot(.x, loadings = TRUE, loadings.label = TRUE,
                 #loadings.label.repel = TRUE,
                 loadings.label.colour = "black",
                 loadings.colour = "gray85",
                 loadings.label.alpha = 0.5,
                 loadings.label.size = 3,
                 loadings.label.hjust = 1.1,
                 frame = TRUE,
                 frame.type = "norm",
                 data = .y, 
                 colour = "Oneota_Present", 
                 shape = "Oneota_Present",
                 frame.level = .9, 
                 frame.alpha = 0.001, 
                 size = 2) +
        theme_bw() + 
        #geom_text(label = .y$Sample) +
        labs(x = "Principal Component 1",
             y = "Principal Component 2",
             title = "First two principal components of PCA on CIRV Ceramic dataset")
    )
  ) %>%
  pull(pca_graph)

ggplotly(oneota_pc1pc2[[1]])
# I certainly can't see any meaningful trends here. This suggests that Oneota and Mississippian 
# potters are almost undoubtedly using similar (or the same) clay. However, more work is needed
# to confirm this hypothesis. 

# Does temper percent matter?
tempperc_pc1pc2 <- sample_pca %>%
  mutate(
    pca_graph = map2(
      .x = pca,
      .y = data,
      ~ autoplot(.x, loadings = TRUE, loadings.label = TRUE,
                 #loadings.label.repel = TRUE,
                 loadings.label.colour = "black",
                 loadings.colour = "gray85",
                 loadings.label.alpha = 0.5,
                 loadings.label.size = 3,
                 loadings.label.hjust = 1.1,
                 frame = TRUE,
                 frame.type = "norm",
                 data = .y, 
                 colour = "Temper_Perc", 
                 shape = "Geography_2",
                 frame.level = .9, 
                 frame.alpha = 0.001, 
                 size = 2) +
        theme_bw() + 
        #geom_text(label = .y$Sample) +
        labs(x = "Principal Component 1",
             y = "Principal Component 2",
             title = "First two principal components of PCA on CIRV Ceramic dataset")
    )
  ) %>%
  pull(pca_graph)

ggplotly(tempperc_pc1pc2[[1]])
# Can't really discern anything here

# Does temper size matter?
tempsize_pc1pc2 <- sample_pca %>%
  mutate(
    pca_graph = map2(
      .x = pca,
      .y = data,
      ~ autoplot(.x, loadings = TRUE, loadings.label = TRUE,
                 #loadings.label.repel = TRUE,
                 loadings.label.colour = "black",
                 loadings.colour = "gray85",
                 loadings.label.alpha = 0.5,
                 loadings.label.size = 3,
                 loadings.label.hjust = 1.1,
                 frame = TRUE,
                 frame.type = "norm",
                 data = .y, 
                 colour = "Temper_Size", 
                 shape = "Geography_2",
                 frame.level = .9, 
                 frame.alpha = 0.001, 
                 size = 2) +
        theme_bw() + 
        #geom_text(label = .y$Sample) +
        labs(x = "Principal Component 1",
             y = "Principal Component 2",
             title = "First two principal components of PCA on CIRV Ceramic dataset")
    )
  ) %>%
  pull(pca_graph)

ggplotly(tempsize_pc1pc2[[1]])
# Interestingly, it appears that the smallest temper size only appears in the northern part of the valley. 
# That might suggest that there is either a preference for smaller temper grains there or it is a 
# response to the clay available in the north. 

# Finally, let's check Cultural Group
culture_pc1pc2 <- sample_pca %>%
  mutate(
    pca_graph = map2(
      .x = pca,
      .y = data,
      ~ autoplot(.x, loadings = TRUE, loadings.label = TRUE,
                 #loadings.label.repel = TRUE,
                 loadings.label.colour = "black",
                 loadings.colour = "gray85",
                 loadings.label.alpha = 0.5,
                 loadings.label.size = 3,
                 loadings.label.hjust = 1.1,
                 frame = TRUE,
                 frame.type = "norm",
                 data = .y, 
                 colour = "Cultural_Group", 
                 shape = "Cultural_Group",
                 frame.level = .9, 
                 frame.alpha = 0.001, 
                 size = 2) +
        theme_bw() + 
        #geom_text(label = .y$Sample) +
        labs(x = "Principal Component 1",
             y = "Principal Component 2",
             title = "First two principal components of PCA on CIRV Ceramic dataset")
    )
  ) %>%
  pull(pca_graph)

ggplotly(culture_pc1pc2[[1]])
# Very little to no discernible patterning here. This again indicates that both cultural
# groups are likely to be using the same clays. 

# Based on the initial inspection of PCs 1 and 2, it looks like three elements in particular 
# are driving some of the group separation (subtle as it is): Mo, Mn, and Si
# Let's plot those three elements in a 3D scatter plot
Mo_Mn_Si <- plot_ly(sample_new_pcaready, x = ~Mo, y = ~Mn, z = ~Si, color = ~Geography_2)

# Explore Samples by date
ggplotly(ggplot(sample_new_pcaready, aes(x = Mo, y = Mg, color = Date)) + 
  stat_ellipse(aes(color = Geography_2)) + geom_text(aes(label = Date), size = 2))

# All in all, only a general trend of the north-south distinction holds when considering prior information
# on PC1-PC2 biplots. That distinction is marked by significant overlap. We'll consider that when
# running group membership probabilities. But first, it's necessary to explore how a variety of 
# statistical methods will group the data. We'll append that group information to our PC list such that
# we can consider both prior information and statistical infomation in groups before moving on to 
# group refinement. 



########################## Cluster Analysis ########################


##### Hierarchical Cluster Analysis #####

# Now that I have a sense of the structure of the ceramic data set based on PCA, the next step
# in compositional analysis is to see how the groups defined from prior information compare
# to groups constructed using statistical clustering methods such as HCA, kmeans, and kmediods

# Let's start with some tree-based methods (aka Hierarchical cluster analysis or HCA)
# We'll use agglomerative methods here (bottom up) as opposed to divisive methods (top down)

# Prep the dataset
# Set rownames to aid in interpretations of dendrograms and other plots
rownames(sample_new_pcaready) <- sample_new_pcaready$Sample

# Drop the prior known information features
sample_new_distready <- sample_new_pcaready %>% select(c(-1:-12))

# First make a dissimilarity matrix based on Euclidean distance
euc_dist_ceramics <- dist(sample_new_distready, method = "euclidean")

# We can check agglomerative coefficients with agnes to see which method(s) might
# work best with the ceramic compositional dataset
clustmethods <- c( "average", "single", "complete", "ward")
names(clustmethods) <- c( "average", "single", "complete", "ward")

# function to compute agglomerative coefficient
ac <- function(x) {
  agnes(euc_dist_ceramics, method = x)$ac
}

map_dbl(clustmethods, ac)

# Looks like complete and Ward linkage methods will work best. We'll run those

# Hierarchical clustering using Ward's Linkage
wardhc1 <- hclust(euc_dist_ceramics, method = "ward.D")
wardhc1_dend <- as.dendrogram(wardhc1) # create dendrogram object

# Plot Ward dendrogram
plot(wardhc1_dend, nodePar = list(lab.cex = 0.15, pch = NA))
# Looks like there are three well defined clusters at a height of 20. 
# We can color the dendrogram at that height
wardhc1_dend_20 <- color_branches(wardhc1_dend, h = 20)
plot(wardhc1_dend_20, cex.axis = 0.75, cex.lab = 0.75, nodePar = list(lab.cex = 0.15, pch = NA))

# This looks like a good hypothetical groupings to add to our original dataset
# We'll add all statistical clusters to a dataset sample_new_stat_clusters
ward_dist_groups <- cutree(wardhc1_dend_20, h = 20)
table(ward_dist_groups) # How many samples are in each cluster>?

sample_new_stat_clusters <- sample_new_pcaready %>%
                              select(Sample) %>%
                              mutate(Ward_HCA_Cluster = ward_dist_groups)

# Visualize the clusters from HCA using Ward's linkage
fviz_cluster(list(data = sample_new_distready, cluster = ward_dist_groups))

# Complete linkage also has a high agglomerative coefficient, let's model it
completehc1 <- hclust(euc_dist_ceramics, method = "complete")
completehc1_dend <- as.dendrogram(completehc1)

# Plot Complete linkage dendrogram cut at 2.4, which results in 6 clusters (3 main and 3 minor)
completehc1_dend_2.4 <- color_branches(completehc1_dend, h = 2.4)
plot(completehc1_dend_2.4, cex.axis = 0.75, cex.lab = 0.75, nodePar = list(lab.cex = 0.15, pch = NA))
complete_dist_groups <- cutree(completehc1_dend, h = 2.4)
table(complete_dist_groups)
sample_new_stat_clusters <- sample_new_stat_clusters %>%
                              mutate(Complete_HCA_Cluster = complete_dist_groups)

# Visualize the clusters from HCA using Complete linkage
fviz_cluster(list(data = sample_new_distready, cluster = complete_dist_groups))

# Let's compare the Ward's and Complete Linkage dendrograms with a tanglegram 
# (this is very resource intensive, so I'm commenting it out)
# tanglegram(wardhc1_dend, completehc1_dend)

# Now let's see how these HCA groups correspond to other clustering methods 


##### K-means Cluster Analysis #####


# First, it's a good idea to use a few methods to assess the number of clusters to model
# Elbow Method
fviz_nbclust(sample_new_distready, kmeans, method = "wss") # 3 - 8 optimal clusters; 3-4 looks good
# Silhouette Method
fviz_nbclust(sample_new_distready, kmeans, method = "silhouette") # 3 optimal clusters
# Gap Stat
fviz_nbclust(sample_new_distready, kmeans, method = "gap_stat") # 1 optimal cluster

# Based on the optimal cluster methods, it looks like we should run kmeans twice, once with 
# 3 clusters and once with 4 clusters

# 3 Cluster K-means
k3 <- kmeans(sample_new_distready, 
             centers = 3, # number of clusters
             nstart = 50, # number of random initial configurations out of which the best one is chosen
             iter.max = 500) # number of allowable iterations allowed 

# Visualize 3 cluster kmeans 
fviz_cluster(k3, data = sample_new_distready)

# Assign to clustering assignments data frame
sample_new_stat_clusters <- sample_new_stat_clusters %>%
                              mutate(Kmeans_3 = k3$cluster)


# 4 Cluster K-means
k4 <- kmeans(sample_new_distready, centers = 4, nstart = 50, iter.max = 500)

# Visualize 4 cluster kmeans
fviz_cluster(k4, data = sample_new_distready)

# Assign to clustering assignments data frame
sample_new_stat_clusters <- sample_new_stat_clusters %>%
                              mutate(Kmeans_4 = k4$cluster)



##### K-mediods Cluster Analysis #####

# For k-mediods, we'll be using the pam function from the cluster package. pam stands for 
# "partitioning around mediods"

# As with k-means, it's a good idea to use a few methods to assess the number of clusters to model
# Elbow Method
fviz_nbclust(sample_new_distready, pam, method = "wss") # 5 looks optimal here
# Silhouette Method
fviz_nbclust(sample_new_distready, pam, method = "silhouette") # 2 optimal clusters
# Gap Stat
fviz_nbclust(sample_new_distready, pam, method = "gap_stat") # 1 optimal cluster

# We'll run two clusters - one with 2 and one with 5
# 2 cluster K-mediods
pam2 <- pam(sample_new_distready, 2)

# Plot 2 cluster k-mediods
fviz_cluster(pam2, data = sample_new_distready)

# 5 cluster K-mediods
pam5 <- pam(sample_new_distready, 5)

# Plot 5 cluster k-mediods
fviz_cluster(pam5, data = sample_new_distready)

# Assign k-mediods results to clustering assignments data frame
sample_new_stat_clusters <- sample_new_stat_clusters %>%
                              mutate(Kmediods_2 = pam2$clustering, 
                                     Kmediods_5 = pam5$clustering)


# One last exploratory metric would be to take the most often occuring group assignment number, 
# the mode
# Little function to calculate the mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Apply this row-wise to the data
mode_assignment <- apply(sample_new_stat_clusters, 1, Mode)


####_________________Begin Mahalanobis distance and membership assignments________________###

# First, create a data frame of the first 12 PC's, which account for 90% of the variability in 
# the elemental data set. This will allow group membership probability assessments with a
# group as small as 14 (or perhaps 13)
pc1to12 <- sample_pca[['pca_aug']][[1]] %>% 
  select(.fittedPC1, .fittedPC2, .fittedPC3, .fittedPC4, .fittedPC5, .fittedPC6, .fittedPC7, 
         .fittedPC8, .fittedPC9, .fittedPC10, .fittedPC11, .fittedPC12)

# This function written by Matt Peeples allows for for calculating group membership probabilities
# by chemical compositional distance using Mahalanobis distances and Hotellings T^2 statistic
# This is identical to the procedure used in MURRAP GAUSS routines for the same purpose and has
# been cross referenced against that routine to ensure accuracy for data presented in this 
# analysis
group.mem.probs <- function(x2.l,attr1.grp,grps) {
  
  # x2.l = transformed element data
  # attr1 = group designation by sample
  # grps <- vector of groups to evaluate
  
  probs <- list()
 for (m in 1:length(grps)) {
    x <- x2.l[which(attr1.grp == grps[m]),]
    probs[[m]] <- matrix(0,nrow(x),length(grps))
    colnames(probs[[m]]) <- grps
    rownames(probs[[m]]) <- rownames(x)
    
    grps2 <- grps[-m]
    
    p.val <- NULL
    for (i in 1:nrow(x)) {p.val[i] <- HotellingsT2(x[i,], x[-i,])$p.value}
    probs[[m]][,m] <- round(p.val,5)*100
    
    for (j in 1:length(grps2)) {
      p.val2 <- NULL
      for (i in 1:nrow(x)) {p.val2[i] <- HotellingsT2(x[i,],x2.l[which(attr1.grp == grps2[j]),])$p.value}
      probs[[m]][,which(grps == grps2[j])] <- round(p.val2, 5)*100}}
  return(probs)
}

########### WARD HCA #################
# Calculate group membership probabilities for the HCA Ward group assignments based on PCA data
ward_group_mem <- group.mem.probs(pc1to12, sample_new_stat_clusters$Ward_HCA_Cluster, 
                    unique(sample_new_stat_clusters$Ward_HCA_Cluster)) 

# Create list of data that is grouped the same as the group probability list
ward_samp_list <- split(sample_new_stat_clusters[, c(1:2)], 
                        f = sample_new_stat_clusters$Ward_HCA_Cluster)

# Convert the list of matrices of group membership probabilities to data frames 
# and bind rows into one data frame
ward_group_mem <- map(ward_group_mem, as.data.frame) %>% bind_rows()

# Convert the list of matrices of sample names to data frames and bind into one data frame
ward_samp_df <- map(ward_samp_list, as.data.frame) %>% bind_rows()

# Bind to initial sample id and group assignment from Ward HCA 
# and convert to data frame for easier handling
ward_group_mem <- as.data.frame(bind_cols(ward_group_mem, ward_samp_df))

# New column of membership probability for initially assigned group
ward_group_mem$assigned_val <- ward_group_mem[1:3][cbind(seq_len(nrow(ward_group_mem)), 
                                                         ward_group_mem$Ward_HCA_Cluster)]

# Set the initial roup assignment value to zero to allow for comparisons
ward_group_mem[cbind(seq_len(nrow(ward_group_mem)), ward_group_mem$Ward_HCA_Cluster)] <- 0

# The heuristic I am using to assess group membership asks whether or not the probability of
# group membership in the original assigned cluster is greater than 10% and that the 
# probability of membership in any other cluster is less that 10%. This follows 
# Peeples (2010) in part and is a fairly conservative threshold. 
ward_group_mem %>%
 # mutate(out_group_sum = `1` + `2` + `3`) %>%
  mutate(new_assign = ifelse(assigned_val > 10 & (`1` < 10 & `2` < 10 & `3` < 10), 
                            Ward_HCA_Cluster, "unassigned")) %>% 
#   filter(new_assign != "unassigned")  
  summarize(perc_unassigned = sum(new_assign == "unassigned")/n() * 100)
# Applying the heuristic to the initial group assignments for the Ward HCA clusters results 
# in an 77.16% unassignment rate. This is quite high. Let's check other methods to see how they fair. 


########### Kmeans 4 #################
# Group probabilities for the kmeans 4 cluster solution on transformed PCA data 
kmean4_group_mem <- group.mem.probs(pc1to12, sample_new_stat_clusters$Kmeans_4, 
                                      unique(sample_new_stat_clusters$Kmeans_4)) 

# Create list of data that is grouped the same as the group probability list
kmean4_samp_list <- split(sample_new_stat_clusters[, c("Sample", "Kmeans_4")], 
                        f = sample_new_stat_clusters$Kmeans_4)

# Convert the matrices of group membership probabilities to data frames and bind rows into one data frame
kmean4_group_mem <- map(kmean4_group_mem, as.data.frame) %>% bind_rows()

# Convert the list of matrices of sample names to data frames and bind into one data frame
kmean4_samp_df <- map(kmean4_samp_list, as.data.frame) %>% bind_rows()

# Bind to initial sample id and group assignment from Kmean 4
# and convert to data frame for easier handling
kmean4_group_mem <- as.data.frame(bind_cols(kmean4_group_mem, kmean4_samp_df))

# Convert to tibble data frame for easier handling
kmean4_group_mem <- as.data.frame(kmean4_group_mem)

# New column of membership probability for initially assigned group
kmean4_group_mem$assigned_val <- kmean4_group_mem[1:4][cbind(seq_len(nrow(kmean4_group_mem)), 
                                                           kmean4_group_mem$Kmeans_4)]

# Set the initial roup assignment value to zero to allow for comparisons
kmean4_group_mem[cbind(seq_len(nrow(kmean4_group_mem)), kmean4_group_mem$Kmeans_4)] <- 0

# Assess membership probabilities using my heuristic
kmean4_group_mem %>% 
  mutate(new_assign = ifelse(assigned_val > 10 & (`1` < 10 & `2` < 10 & `3` < 10 & `4` < 10), 
                             Kmeans_4, "unassigned")) %>% 
  #   filter(new_assign != "unassigned")  
  summarize(perc_unassigned = sum(new_assign == "unassigned")/n() * 100)
# At an 89.69%, it doesn't seem like kmeans 4 group clusters faired much better than Ward HCA


########### Kmediods (pam) 5 #################
# Group probabilities for the kmediods (pam) 5 cluster solution on PC's 1 to 12 (90% of variability)
kmed5_group_mem <- group.mem.probs(pc1to12, sample_new_stat_clusters$Kmediods_5, 
                                    unique(sample_new_stat_clusters$Kmediods_5)) 

# Create list of data that is grouped the same as the group probability list
kmed5_samp_list <- split(sample_new_stat_clusters[, c("Sample", "Kmediods_5")], 
                          f = sample_new_stat_clusters$Kmediods_5)

# Convert the matrices of group membership probabilities to data frames and bind rows into one data frame
kmed5_group_mem <- map(kmed5_group_mem, as.data.frame) %>% bind_rows()

# Convert the list of matrices of sample names to data frames and bind into one data frame
kmed5_samp_df <- map(kmed5_samp_list, as.data.frame) %>% bind_rows()

# Bind to initial sample id and group assignment from Kmed 5
# and convert to data frame for easier handling
kmed5_group_mem <- as.data.frame(bind_cols(kmed5_group_mem, kmed5_samp_df))

# New column of membership probability for initially assigned group
kmed5_group_mem$assigned_val <- kmed5_group_mem[1:5][cbind(seq_len(nrow(kmed5_group_mem)), 
                                                           kmed5_group_mem$Kmediods_5)]

# Set the initial roup assignment value to zero to allow for comparisons
kmed5_group_mem[cbind(seq_len(nrow(kmed5_group_mem)), kmed5_group_mem$Kmediods_5)] <- 0

# Assess membership probabilities using my heuristic
kmed5_group_mem %>% 
  mutate(new_assign = ifelse(assigned_val > 10 & (`1` < 10 & `2` < 10 & `3` < 10 & `4` < 10 & `5` < 10), 
                             Kmediods_5, "unassigned")) %>% 
  #   filter(new_assign != "unassigned")  
  summarize(perc_unassigned = sum(new_assign == "unassigned")/n() * 100)
# Ouch, at 95.58% unassigned using the heuristic criteria, this doesn't hold up



########### Kmediods (pam) 2 #################
# Group probabilities for the kmediods (pam) 2 cluster solution on PC's 1 to 12 (90% of variability)
kmed2_group_mem <- group.mem.probs(pc1to12, sample_new_stat_clusters$Kmediods_2, 
                                   unique(sample_new_stat_clusters$Kmediods_2)) 

# Create list of data that is grouped the same as the group probability list
kmed2_samp_list <- split(sample_new_stat_clusters[, c("Sample", "Kmediods_2")], 
                         f = sample_new_stat_clusters$Kmediods_2)

# Convert the matrices of group membership probabilities to data frames and bind rows into one data frame
kmed2_group_mem <- map(kmed2_group_mem, as.data.frame) %>% bind_rows()

# Convert the list of matrices of sample names to data frames and bind into one data frame
kmed2_samp_df <- map(kmed2_samp_list, as.data.frame) %>% bind_rows()

# Bind to initial sample id and group assignment from Kmed 5
# and convert to data frame for easier handling
kmed2_group_mem <- as.data.frame(bind_cols(kmed2_group_mem, kmed2_samp_df))

# New column of membership probability for initially assigned group
kmed2_group_mem$assigned_val <- kmed2_group_mem[1:2][cbind(seq_len(nrow(kmed2_group_mem)), 
                                                           kmed2_group_mem$Kmediods_2)]

# Set the initial roup assignment value to zero to allow for comparisons
kmed2_group_mem[cbind(seq_len(nrow(kmed2_group_mem)), kmed2_group_mem$Kmediods_2)] <- 0

# Assess membership probabilities using my heuristic
kmed2_group_mem %>% 
  mutate(new_assign = ifelse(assigned_val > 2.5 & (`1` < 10 & `2` < 10), 
                             Kmediods_2, "unassigned")) %>% 
 # filter(assigned_val < `1` | assigned_val < `2`)  
  summarize(perc_unassigned = sum(new_assign == "unassigned")/n() * 100)
# A 75.32% unassigned using the heuristic criteria is better, but still doesn't hold up


########### Mahalanobis-first route ##########################################################
########### Core and Unassigned Group Assignments ############################################
# Another common method used for constructing core chemical compositional groups in
# archaeology is to initially treat the entire data set as one large group and iteratively
# removing samples with a membership probability of less than 1%

# Double the PC data so the group can be compared to itself
pc1to12_twice <- bind_rows(pc1to12, pc1to12)

# Double the stat cluster assignment data
sample_new_stat_clusters_twice <- bind_rows(sample_new_stat_clusters, sample_new_stat_clusters)

# Create vector of group assignments 
one_two <- c(rep(1, 543), rep(2, 543))

# Bind group assignments to cluster data
sample_new_stat_clusters_twice <- cbind(sample_new_stat_clusters_twice, one_two)

# Group probabilities for the group as one data set on PC's 1 through 12
one_group_mem <- group.mem.probs(pc1to12_twice, sample_new_stat_clusters_twice$one_two, 
                                    unique(sample_new_stat_clusters_twice$one_two)) 

# Create list of data that is grouped the same as the group probability list
one_samp_list <- split(sample_new_stat_clusters_twice[, c("Sample", "one_two")], 
                         f = sample_new_stat_clusters_twice$one_two)

# Convert the matrices of group membership probabilities to data frames 
# and bind rows into one data frame
one_group_mem <- map(one_group_mem, as.data.frame) %>% bind_rows()

# Convert the list of matrices of sample names to data frames and bind into one data frame
one_samp_df <- map(one_samp_list, as.data.frame) %>% bind_rows()

# Bind to initial sample id and group assignment from 
# and convert to data frame for easier handling
one_group_mem <- as.data.frame(bind_cols(one_group_mem, one_samp_df))

# Create data frame of sample to retain after first iteraction
iter1 <- one_group_mem %>%
          filter(one_two == 1) %>%
          filter(`1` > 1) %>%
          select(Sample) 

# Create data frame of unassigned samples after first iteraction
iter1_unassigned <- one_group_mem %>%
                      filter(one_two == 1) %>%
                      filter(`1` < 1) %>%
                      select(Sample) 
  
# Subset initial groups
one_group_mem1 <- one_group_mem %>%
                    filter(Sample %in% iter1$Sample)


### Iteration two 
# Bind samples list to PCA data, filter out the unassigned samples after iteration one 
# and select PC data only for group membership probability calculation
pc1to12_twice_iter2 <- bind_cols(sample_new_stat_clusters_twice[, c("Sample", "one_two")], 
                                  pc1to12_twice) %>%
                          filter(Sample %in% iter1$Sample) %>%
                          select(-Sample, -one_two)

# Prep the sample names and assignments for iteration 2                
sample_new_stat_clusters_twice_iter2 <- sample_new_stat_clusters_twice[, c("Sample", "one_two")] %>%
                                          filter(Sample %in% iter1$Sample)

# Group probabilities for iteration 2 of the group as one data set on PC's 1 through 12
one_group_mem_iter_2 <- group.mem.probs(pc1to12_twice_iter2, 
                                        sample_new_stat_clusters_twice_iter2$one_two, 
                                 unique(sample_new_stat_clusters_twice_iter2$one_two)) 

# Create list of data that is grouped the same as the group probability list
one_samp_list_iter2 <- split(sample_new_stat_clusters_twice_iter2, 
                          f = sample_new_stat_clusters_twice_iter2$one_two)

# Convert the matrices of group membership probabilities to data frames 
# and bind rows into one data frame
one_group_mem_iter2 <- map(one_group_mem_iter_2, as.data.frame) %>% bind_rows()

# Convert the list of matrices of sample names to data frames and bind into one data frame
one_samp_df_iter2 <- map(one_samp_list_iter2, as.data.frame) %>% bind_rows()

# Bind to initial sample id and group assignment from 
# and convert to data frame for easier handling
one_group_mem_iter2 <- as.data.frame(bind_cols(one_group_mem_iter2, one_samp_df_iter2))

# Create data frame of sample to retain after first iteraction
iter2 <- one_group_mem_iter2 %>%
            filter(one_two == 1) %>%
            filter(`1` > 1) %>%
            select(Sample) 

# Create data frame of unassigned samples after first iteraction
iter2_unassigned <- one_group_mem_iter2 %>%
                      filter(one_two == 1) %>%
                      filter(`1` < 1) %>%
                      select(Sample) 

# Subset initial groups
one_group_mem2 <- one_group_mem_iter2 %>%
                     filter(Sample %in% iter2$Sample)


### Iteration 3
# Bind samples list to PCA data, filter out the unassigned samples after iteration two 
# and select PC data only for group membership probability calculation
pc1to12_twice_iter3 <- bind_cols(sample_new_stat_clusters_twice_iter2, 
                                 pc1to12_twice_iter2) %>%
                            filter(Sample %in% iter2$Sample) %>%
                            select(-Sample, -one_two)

# Prep the sample names and assignments for iteration 3               
sample_new_stat_clusters_twice_iter3 <- sample_new_stat_clusters_twice_iter2 %>%
                                            filter(Sample %in% iter2$Sample)

# Group probabilities for iteration 3 of the group as one data set on PC's 1 through 12
one_group_mem_iter_3 <- group.mem.probs(pc1to12_twice_iter3, 
                                        sample_new_stat_clusters_twice_iter3$one_two, 
                                        unique(sample_new_stat_clusters_twice_iter3$one_two)) 

# Create list of data that is grouped the same as the group probability list
one_samp_list_iter3 <- split(sample_new_stat_clusters_twice_iter3, 
                             f = sample_new_stat_clusters_twice_iter3$one_two)

# Convert the matrices of group membership probabilities to data frames 
# and bind rows into one data frame
one_group_mem_iter3 <- map(one_group_mem_iter_3, as.data.frame) %>% bind_rows()

# Convert the list of matrices of sample names to data frames and bind into one data frame
one_samp_df_iter3 <- map(one_samp_list_iter3, as.data.frame) %>% bind_rows()

# Bind to initial sample id and group assignment 
# and convert to data frame for easier handling
one_group_mem_iter3 <- as.data.frame(bind_cols(one_group_mem_iter3, one_samp_df_iter3))

# Create data frame of sample to retain after third iteraction
iter3 <- one_group_mem_iter3 %>%
            filter(one_two == 1) %>%
            filter(`1` > 1) %>%
            select(Sample) 

# Create data frame of unassigned samples after third iteraction
iter3_unassigned <- one_group_mem_iter3 %>%
                      filter(one_two == 1) %>%
                      filter(`1` < 1) %>%
                      select(Sample) 

# Subset initial groups
one_group_mem3 <- one_group_mem_iter3 %>%
                    filter(Sample %in% iter3$Sample)

### Iteration 4
# Bind samples list to PCA data, filter out the unassigned samples after iteration three 
# and select PC data only for group membership probability calculation
pc1to12_twice_iter4 <- bind_cols(sample_new_stat_clusters_twice_iter3, 
                                 pc1to12_twice_iter3) %>%
                          filter(Sample %in% iter3$Sample) %>%
                          select(-Sample, -one_two)

# Prep the sample names and assignments for iteration 4             
sample_new_stat_clusters_twice_iter4 <- sample_new_stat_clusters_twice_iter3 %>%
                                           filter(Sample %in% iter3$Sample)

# Group probabilities for iteration 4 of the group as one data set on PC's 1 through 12
one_group_mem_iter_4 <- group.mem.probs(pc1to12_twice_iter4, 
                                        sample_new_stat_clusters_twice_iter4$one_two, 
                                        unique(sample_new_stat_clusters_twice_iter4$one_two)) 

# Create list of data that is grouped the same as the group probability list
one_samp_list_iter4 <- split(sample_new_stat_clusters_twice_iter4, 
                             f = sample_new_stat_clusters_twice_iter4$one_two)

# Convert the matrices of group membership probabilities to data frames 
# and bind rows into one data frame
one_group_mem_iter4 <- map(one_group_mem_iter_4, as.data.frame) %>% bind_rows()

# Convert the list of matrices of sample names to data frames and bind into one data frame
one_samp_df_iter4 <- map(one_samp_list_iter4, as.data.frame) %>% bind_rows()

# Bind to initial sample id and group assignment 
# and convert to data frame for easier handling
one_group_mem_iter4 <- as.data.frame(bind_cols(one_group_mem_iter4, one_samp_df_iter4))

# Create data frame of sample to retain after fourth iteraction
iter4 <- one_group_mem_iter4 %>%
            filter(one_two == 1) %>%
            filter(`1` > 1) %>%
            select(Sample) 

# Create data frame of unassigned samples after fourth iteraction
iter4_unassigned <- one_group_mem_iter4 %>%
                      filter(one_two == 1) %>%
                      filter(`1` < 1) %>%
                      select(Sample) 

# Subset initial groups
one_group_mem4 <- one_group_mem_iter4 %>%
                    filter(Sample %in% iter4$Sample)


### Iteration 5
# Bind samples list to PCA data, filter out the unassigned samples after iteration four 
# and select PC data only for group membership probability calculation
pc1to12_twice_iter5 <- bind_cols(sample_new_stat_clusters_twice_iter4, 
                                 pc1to12_twice_iter4) %>%
                          filter(Sample %in% iter4$Sample) %>%
                          select(-Sample, -one_two)

# Prep the sample names and assignments for iteration 5             
sample_new_stat_clusters_twice_iter5 <- sample_new_stat_clusters_twice_iter4 %>%
                                           filter(Sample %in% iter4$Sample)

# Group probabilities for iteration 5 of the group as one data set on PC's 1 through 12
one_group_mem_iter_5 <- group.mem.probs(pc1to12_twice_iter5, 
                                        sample_new_stat_clusters_twice_iter5$one_two, 
                                        unique(sample_new_stat_clusters_twice_iter5$one_two)) 

# Create list of data that is grouped the same as the group probability list
one_samp_list_iter5 <- split(sample_new_stat_clusters_twice_iter5, 
                             f = sample_new_stat_clusters_twice_iter5$one_two)

# Convert the matrices of group membership probabilities to data frames 
# and bind rows into one data frame
one_group_mem_iter5 <- map(one_group_mem_iter_5, as.data.frame) %>% bind_rows()

# Convert the list of matrices of sample names to data frames and bind into one data frame
one_samp_df_iter5 <- map(one_samp_list_iter5, as.data.frame) %>% bind_rows()

# Bind to initial sample id and group assignment 
# and convert to data frame for easier handling
one_group_mem_iter5 <- as.data.frame(bind_cols(one_group_mem_iter5, one_samp_df_iter5))

# Create data frame of sample to retain after fifth iteraction
iter5 <- one_group_mem_iter5 %>%
          filter(one_two == 1) %>%
          filter(`1` > 1) %>%
          select(Sample) 

# Create data frame of unassigned samples after fifth iteraction
iter5_unassigned <- one_group_mem_iter5 %>%
                      filter(one_two == 1) %>%
                      filter(`1` < 1) %>%
                      select(Sample) 

# Subset initial groups
one_group_mem5 <- one_group_mem_iter5 %>%
                    filter(Sample %in% iter5$Sample)

### Iteration 6
# Bind samples list to PCA data, filter out the unassigned samples after iteration five 
# and select PC data only for group membership probability calculation
pc1to12_twice_iter6 <- bind_cols(sample_new_stat_clusters_twice_iter5, 
                                 pc1to12_twice_iter5) %>%
                                filter(Sample %in% iter5$Sample) %>%
                                select(-Sample, -one_two)

# Prep the sample names and assignments for iteration 6            
sample_new_stat_clusters_twice_iter6 <- sample_new_stat_clusters_twice_iter5 %>%
                                          filter(Sample %in% iter5$Sample)

# Group probabilities for iteration 6 of the group as one data set on PC's 1 through 12
one_group_mem_iter_6 <- group.mem.probs(pc1to12_twice_iter6, 
                                        sample_new_stat_clusters_twice_iter6$one_two, 
                                        unique(sample_new_stat_clusters_twice_iter6$one_two)) 

# Create list of data that is grouped the same as the group probability list
one_samp_list_iter6 <- split(sample_new_stat_clusters_twice_iter6, 
                             f = sample_new_stat_clusters_twice_iter6$one_two)

# Convert the matrices of group membership probabilities to data frames 
# and bind rows into one data frame
one_group_mem_iter6 <- map(one_group_mem_iter_6, as.data.frame) %>% bind_rows()

# Convert the list of matrices of sample names to data frames and bind into one data frame
one_samp_df_iter6 <- map(one_samp_list_iter6, as.data.frame) %>% bind_rows()

# Bind to initial sample id and group assignment 
# and convert to data frame for easier handling
one_group_mem_iter6 <- as.data.frame(bind_cols(one_group_mem_iter6, one_samp_df_iter6))

# Create data frame of sample to retain after sixth iteraction
iter6 <- one_group_mem_iter6 %>%
            filter(one_two == 1) %>%
            filter(`1` > 1) %>%
            select(Sample) 

# Create data frame of unassigned samples after sixth iteraction
iter6_unassigned <- one_group_mem_iter6 %>%
                      filter(one_two == 1) %>%
                      filter(`1` < 1) %>%
                      select(Sample) 

# Subset initial groups
one_group_mem6 <- one_group_mem_iter6 %>%
                    filter(Sample %in% iter6$Sample)

### Iteration 7
# Bind samples list to PCA data, filter out the unassigned samples after iteration six 
# and select PC data only for group membership probability calculation
pc1to12_twice_iter7 <- bind_cols(sample_new_stat_clusters_twice_iter6, 
                                 pc1to12_twice_iter6) %>%
                          filter(Sample %in% iter6$Sample) %>%
                          select(-Sample, -one_two)

# Prep the sample names and assignments for iteration 6            
sample_new_stat_clusters_twice_iter7 <- sample_new_stat_clusters_twice_iter6 %>%
                                          filter(Sample %in% iter6$Sample)

# Group probabilities for iteration 7 of the group as one data set on PC's 1 through 12
one_group_mem_iter_7 <- group.mem.probs(pc1to12_twice_iter7,
                                        sample_new_stat_clusters_twice_iter7$one_two, 
                                        unique(sample_new_stat_clusters_twice_iter7$one_two)) 

# Create list of data that is grouped the same as the group probability list
one_samp_list_iter7 <- split(sample_new_stat_clusters_twice_iter7, 
                             f = sample_new_stat_clusters_twice_iter7$one_two)

# Convert the matrices of group membership probabilities to data frames 
# and bind rows into one data frame
one_group_mem_iter7 <- map(one_group_mem_iter_7, as.data.frame) %>% bind_rows()

# Convert the list of matrices of sample names to data frames and bind into one data frame
one_samp_df_iter7 <- map(one_samp_list_iter7, as.data.frame) %>% bind_rows()

# Bind to initial sample id and group assignment 
# and convert to data frame for easier handling
one_group_mem_iter7 <- as.data.frame(bind_cols(one_group_mem_iter7, one_samp_df_iter7))

# Create data frame of sample to retain after seventh iteraction
iter7 <- one_group_mem_iter7 %>%
          filter(one_two == 1) %>%
          filter(`1` > 1) %>%
          select(Sample) 

# Create data frame of unassigned samples after seventh iteraction
iter7_unassigned <- one_group_mem_iter7 %>%
                      filter(one_two == 1) %>%
                      filter(`1` < 1) %>%
                      select(Sample) 

# Subset initial groups
one_group_mem7 <- one_group_mem_iter7 %>%
                     filter(Sample %in% iter7$Sample)

### Iteration 8
# Bind samples list to PCA data, filter out the unassigned samples after iteration seven 
# and select PC data only for group membership probability calculation
pc1to12_twice_iter8 <- bind_cols(sample_new_stat_clusters_twice_iter7, 
                                 pc1to12_twice_iter7) %>%
                        filter(Sample %in% iter7$Sample) %>%
                        select(-Sample, -one_two)

# Prep the sample names and assignments for iteration 7           
sample_new_stat_clusters_twice_iter8 <- sample_new_stat_clusters_twice_iter7 %>%
                                           filter(Sample %in% iter7$Sample)

# Group probabilities for iteration 8 of the group as one data set on PC's 1 through 12
one_group_mem_iter_8 <- group.mem.probs(pc1to12_twice_iter8, 
                                        sample_new_stat_clusters_twice_iter8$one_two, 
                                        unique(sample_new_stat_clusters_twice_iter8$one_two)) 

# Create list of data that is grouped the same as the group probability list
one_samp_list_iter8 <- split(sample_new_stat_clusters_twice_iter8, 
                             f = sample_new_stat_clusters_twice_iter8$one_two)

# Convert the matrices of group membership probabilities to data frames 
# and bind rows into one data frame
one_group_mem_iter8 <- map(one_group_mem_iter_8, as.data.frame) %>% bind_rows()

# Convert the list of matrices of sample names to data frames and bind into one data frame
one_samp_df_iter8 <- map(one_samp_list_iter8, as.data.frame) %>% bind_rows()

# Bind to initial sample id and group assignment 
# and convert to data frame for easier handling
one_group_mem_iter8 <- as.data.frame(bind_cols(one_group_mem_iter8, one_samp_df_iter8))

# Create data frame of sample to retain after fifth iteraction
iter8 <- one_group_mem_iter8 %>%
          filter(one_two == 1) %>%
          filter(`1` > 1) %>%
          select(Sample) 

# Create data frame of unassigned samples after fifth iteraction
iter8_unassigned <- one_group_mem_iter8 %>%
                      filter(one_two == 1) %>%
                      filter(`1` < 1) %>%
                      select(Sample) 

# Subset initial groups
one_group_mem8 <- one_group_mem_iter8 %>%
                     filter(Sample %in% iter8$Sample)


### Iteration 9
# Bind samples list to PCA data, filter out the unassigned samples after iteration eight 
# and select PC data only for group membership probability calculation
pc1to12_twice_iter9 <- bind_cols(sample_new_stat_clusters_twice_iter8, 
                                 pc1to12_twice_iter8) %>%
                          filter(Sample %in% iter8$Sample) %>%
                          select(-Sample, -one_two)

# Prep the sample names and assignments for iteration 8      
sample_new_stat_clusters_twice_iter9 <- sample_new_stat_clusters_twice_iter8 %>%
                                          filter(Sample %in% iter8$Sample)

# Group probabilities for iteration 9 of the group as one data set on PC's 1 through 12
one_group_mem_iter_9 <- group.mem.probs(pc1to12_twice_iter9, 
                                        sample_new_stat_clusters_twice_iter9$one_two, 
                                        unique(sample_new_stat_clusters_twice_iter9$one_two)) 

# Create list of data that is grouped the same as the group probability list
one_samp_list_iter9 <- split(sample_new_stat_clusters_twice_iter9, 
                             f = sample_new_stat_clusters_twice_iter9$one_two)

# Convert the matrices of group membership probabilities to data frames 
# and bind rows into one data frame
one_group_mem_iter9 <- map(one_group_mem_iter_9, as.data.frame) %>% bind_rows()

# Convert the list of matrices of sample names to data frames and bind into one data frame
one_samp_df_iter9 <- map(one_samp_list_iter9, as.data.frame) %>% bind_rows()

# Bind to initial sample id and group assignment 
# and convert to data frame for easier handling
one_group_mem_iter9 <- as.data.frame(bind_cols(one_group_mem_iter9, one_samp_df_iter9))

# Create data frame of sample to retain after fifth iteraction
iter9 <- one_group_mem_iter9 %>%
            filter(one_two == 1) %>%
            filter(`1` > 1) %>%
            select(Sample) 

# Create data frame of unassigned samples after fifth iteraction
iter9_unassigned <- one_group_mem_iter9 %>%
                      filter(one_two == 1) %>%
                      filter(`1` < 1) %>%
                      select(Sample) 

# Subset initial groups
one_group_mem9 <- one_group_mem_iter9 %>%
                   filter(Sample %in% iter9$Sample)

# Data frame of unassigned samples
maha_unassigned <- bind_rows(iter1_unassigned, iter2_unassigned, iter3_unassigned, 
                             iter4_unassigned, iter5_unassigned, iter6_unassigned, 
                             iter7_unassigned, iter8_unassigned) %>% 
                      arrange(Sample) %>%
                      mutate(one_two = 2)

###### End of Core-Unassigned membership iterations #####

# Now that I have a core group and an unassigned group, it's important to assess whether or not
# any of the unassigned samples might warrant inclusion back into the core group.
# To do this, the unassigned samples will be projected against the core group as before. 
# Defined PC loadings for core and unassigned samples
pc1to12_core_unassigned <-  sample_new_stat_clusters_twice_iter9 %>%
                              filter(one_two == 1) %>%
                              bind_rows(maha_unassigned) %>%
                              left_join(sample_pca[['pca_aug']][[1]], by = "Sample") %>% 
                              select(.fittedPC1:.fittedPC12)
          
# Prep the sample names and assignments for core|unassigned evaluation 
sample_core_unassigned_clusters <-  sample_new_stat_clusters_twice_iter9 %>%
                                      filter(one_two == 1) %>%
                                      bind_rows(maha_unassigned) 

# Group probabilities for iteration 9 of the group as one data set on PC's 1 through 12
core_unassigned_group_prob <- group.mem.probs(pc1to12_core_unassigned, 
                                              sample_core_unassigned_clusters$one_two, 
                                        unique(sample_core_unassigned_clusters$one_two)) 

# Create list of data that is grouped the same as the group probability list
core_unassigned_list <- split(sample_core_unassigned_clusters, 
                             f = sample_core_unassigned_clusters$one_two)

# Convert the matrices of group membership probabilities to data frames and bind rows into one data frame
core_unassigned_group_prob <- map(core_unassigned_group_prob, as.data.frame) %>% bind_rows()

# Convert the list of matrices of sample names to data frames and bind into one data frame
core_unassigned_df <- map(core_unassigned_list, as.data.frame) %>% bind_rows()

# Bind to initial sample id and group assignment 
# and convert to data frame for easier handling
core_unassigned_group_prob <- as.data.frame(bind_cols(core_unassigned_group_prob, core_unassigned_df))

# Check to see if there are any unassigned above the 1% threshold for membership in the core
core_unassigned_group_prob %>%
  filter(one_two == 2 & `1` > 1)
# Does not appear to be the case

# Check to see if there are any core samples below 1% threshold of being assigned to the core 
core_unassigned_group_prob %>%
  filter(one_two == 1 & `1` < 1)
# Also does not appear to be the case. This confirms that we have statistically robust core and 
# unassigned groups. 

# Create interactive 3D scatter plot showing first three PC's and the core and unassigned samples      
p <-  sample_new_stat_clusters_twice_iter9 %>%
        filter(one_two == 1) %>%
        bind_rows(maha_unassigned) %>%
        left_join(sample_pca[['pca_aug']][[1]], by = "Sample") %>%
        mutate(one_two = factor(one_two, labels = c("Core", "Unassigned"))) %>%
        mutate(symbols1 = ifelse(one_two == "Core", "plus", "triangle-up")) %>%
       # ggplot(aes(x = .fittedPC1, y = .fittedPC3, color = one_two)) + geom_point() 
        plot_ly(type = "scatter3d", x = ~.fittedPC1, y = ~.fittedPC2, z = ~.fittedPC3, 
                color = ~as.factor(one_two), size = 3, colors = c('grey40', 'black'), alpha = 0.8,
                text = ~(paste("Sample ID", Sample, '<br>Site:', Site, "<br>Geography_2:", 
                               Geography_2, "<br>Time:", Time, "<br>Cultural Group:", Cultural_Group)),
             #  marker = list(symbol = ~I(symbols1)), size = .3,
                symbol = ~one_two, #symbols = ~symbols1,
              mode = "markers") %>%
        layout(scene = list(xaxis = list(title = 'Principal Component 1'),
                            yaxis = list(title = 'Principal Component 2'),
                            zaxis = list(title = 'Principal Component 3')))

# Adjust plot features
pb <- plotly_build(p)
pb$x$data[[1]]$marker$symbol <- 'diamond-open'
pb$x$data[[2]]$marker$symbol <- 'circle-open'
pb # Display interactive 3D scattergram

# Table of core and unassigned group membership
sample_new_stat_clusters_twice_iter9 %>%
  filter(one_two == 1) %>%
  bind_rows(maha_unassigned) %>%
  left_join(sample_pca[['pca_aug']][[1]], by = "Sample") %>%
  mutate(one_two = factor(one_two, labels = c("Core", "Unassigned"))) %>%
  select(one_two) %>%
  table()


############################## Unassigned Group Structure #########################################
# There are 127 unassigned samples (or 23.4% of the original ceramic sample)
unassigned <- maha_unassigned %>%
                left_join(sample_pca[['pca_aug']][[1]], by = "Sample") %>%
                select(-starts_with(".")) # Drop PC's from full data set PCA

# In taking a look at a table of the sites from where the outliers were recovered, it looks like
# five sites in particular have outlier vessels: Crable, Morton Village, Orendorf C, Ten Mile Creek, 
# and Walsh
table(unassigned$Site)
# Looking through the other pieces of a prior information, there don't appear to be any "smoking-gun"
# trends that may help guide cluster analysis of the Unassigned group

# Prepare samples for distance and clustering methods
unassigned_distready <- unassigned %>%
                          select(Si:Th)

##### Kmeans of Unassigned #####
# First, it's a good idea to use a few methods to assess the number of clusters to model
# Elbow Method
fviz_nbclust(unassigned_distready, kmeans, method = "wss") # 4 - 8 optimal clusters; 4-5 looks good
# Silhouette Method
fviz_nbclust(unassigned_distready, kmeans, method = "silhouette") # 2 optimal clusters
# Gap Stat
fviz_nbclust(unassigned_distready, kmeans, method = "gap_stat") # 1 optimal cluster

# Based on the optimal cluster methods, it looks like we should run kmeans twice, once with 
# 2 clusters and once with 5 clusters

# 2 Cluster K-means
unassigned_k2 <- kmeans(unassigned_distready, 
                         centers = 2, # number of clusters
                         nstart = 50, # number of random initial configs out of which best is chosen
                         iter.max = 500) # number of allowable iterations allowed 

# Visualize 2 cluster kmeans 
fviz_cluster(unassigned_k2, data = unassigned_distready)

# Assign to clustering assignments data frame
unassigned_stat_clusters <- maha_unassigned %>%
                                select(Sample) %>%
                                mutate(Kmeans_2 = unassigned_k2$cluster)


# 5 Cluster K-means
unassigned_k5 <- kmeans(unassigned_distready, centers = 5, nstart = 50, iter.max = 500)

# Visualize 5 cluster kmeans
fviz_cluster(unassigned_k5, data = unassigned_distready)

# Assign to clustering assignments data frame
unassigned_stat_clusters <- unassigned_stat_clusters %>%
                              mutate(Kmeans_5 = unassigned_k5$cluster, 
                                     Kmeans_2 = unassigned_k2$cluster)


##### K-mediods of Unassigned #####

# For k-mediods, we'll be using the pam function from the cluster package. pam stands for 
# "partitioning around mediods"

# As with k-means, it's a good idea to use a few methods to assess the number of clusters to model
# Elbow Method
fviz_nbclust(unassigned_distready, pam, method = "wss") # 5 looks optimal here
# Silhouette Method
fviz_nbclust(unassigned_distready, pam, method = "silhouette") # 2 optimal clusters
# Gap Stat
fviz_nbclust(unassigned_distready, pam, method = "gap_stat") # 1 optimal cluster

# We'll run two clusters - one with 2 and one with 5
# 2 cluster K-mediods
pam2_unassigned <- pam(unassigned_distready, 2)

# Plot 2 cluster k-mediods
fviz_cluster(pam2_unassigned, data = unassigned_distready)

# 5 cluster K-mediods
pam5_unassigned <- pam(unassigned_distready, 5)

# Plot 5 cluster k-mediods
fviz_cluster(pam5_unassigned, data = unassigned_distready)

# Assign k-mediods results to clustering assignments data frame
unassigned_stat_clusters <- unassigned_stat_clusters %>%
                                mutate(Kmediods_2 = pam2_unassigned$clustering, 
                                       Kmediods_5 = pam5_unassigned$clustering)

# There appears to be fairly broad agreement between kmeans and kmediods about the different
# clusters present, but it is important to see how these hold up to comparison using visual inspection

# Convert all unassigned statistical cluster assignments to character for joining
unassigned_stat_clusters[,2:5] <- sapply(unassigned_stat_clusters[,2:5], as.character)

# Make data frame with core sample assignments and unassigned cluster assignments
core_and_unassigned_clusters <- sample_new_stat_clusters_twice_iter9 %>%
                                  filter(one_two == 1) %>%
                                  mutate(Kmeans_2 = "Core", 
                                         Kmeans_5 = "Core", 
                                         Kmediods_2 = "Core",
                                         Kmediods_5 = "Core") %>%
                                  select(-one_two) %>%
                                  bind_rows(unassigned_stat_clusters)

# Join the core assignments to the original PCA data, which is stored in a nested prcomp list object
sample_pca[["data"]][[1]] <- left_join(sample_pca[["data"]][[1]], core_and_unassigned_clusters, by = "Sample")

# Join the core assignments to the augmented PCA data, which is stored in a nested prcomp list object
sample_pca[["pca_aug"]][[1]] <- left_join(sample_pca[["pca_aug"]][[1]], 
                                          core_and_unassigned_clusters, by = "Sample")

# Create column to apply alpha to core group points in biplots for easier interpretation
sample_pca[["data"]][[1]] <- sample_pca[["data"]][[1]] %>%
                                 mutate(alpha = ifelse(Kmeans_5 == "Core", 0.25, 1)) %>%
                                 mutate(alpha = as.vector(alpha))

# Vectorize the alpha column
core_alpha <- as.vector(sample_pca[["data"]][[1]]$alpha)

# Create plot of PC 1 and PC 2 with the 90% conf intervals around the core and outgroups
unass_pc1pc2_kmean2 <- sample_pca %>%
  mutate(
    pca_graph = map2(
      .x = pca,
      .y = data,
      ~ autoplot(.x, loadings = TRUE, loadings.label = TRUE,
                 loadings.label.repel = TRUE,
                 loadings.label.colour = "black",
                 loadings.colour = "gray45",
                 loadings.label.alpha = 0.9,
                 loadings.label.size = 3.5,
                 loadings.label.hjust = -0.5,
                 frame = TRUE,
                 frame.type = "norm",
                 data = .y, 
                 colour = "Kmeans_5", 
                 shape = "Kmeans_5",
                 frame.level = .9, 
                 frame.alpha = 0.001, 
                 size = 2,
                 alpha = core_alpha) +
        theme_bw() + 
       # geom_text(label = .y$Sample) +
        labs(x = "Principal Component 1",
             y = "Principal Component 2")
    )
  ) %>%
  pull(pca_graph)


unass_pc1pc2_kmean2[[1]] + scale_fill_manual(values = c("black","black", "black", 
                                                        "black", "black", "black")) + 
  scale_color_manual(values = c("black","black","black", "black", "black", "black")) +
  scale_shape_manual(values=c(3, 18, 16, 2, 43, 1)) 


###### Shiny app to biplot the various elements and PCs against one another #####
##   UI   ##
ui_sample <- fluidPage(
  pageWithSidebar (
    headerPanel('Bivariate Plotting'),
    sidebarPanel(
      selectInput('x', 'X Variable', names(sample_pca[["pca_aug"]][[1]]), 
                  selected = names(sample_pca[["pca_aug"]][[1]])[[14]]),
      selectInput('y', 'Y Variable', names(sample_pca[["pca_aug"]][[1]]),
                  selected = names(sample_pca[["pca_aug"]][[1]])[[15]]),
      selectInput('color', 'Color', names(sample_pca[["pca_aug"]][[1]]),
                  selected = names(sample_pca[["pca_aug"]][[1]])[[103]]),
      #Slider for plot height
      sliderInput('plotHeight', 'Height of plot (in pixels)', 
                  min = 100, max = 2000, value = 550)
    ),
    mainPanel(
      plotlyOutput('plot1')
    )
  )
)

## Server ##
server_sample <- function(input, output, session) {
  
  # Combine the selected variables into a new data frame
  selectedData <- reactive({
    sample_pca[["pca_aug"]][[1]][, c(input$x, input$y, input$color)]
  })
  
  output$plot1 <- renderPlotly({
    
    #Build plot with ggplot syntax 
    p <- ggplot(data = sample_pca[["pca_aug"]][[1]], aes_string(x = input$x, 
                                                 y = input$y, 
                                                 color = input$color, 
                                                 shape = input$color)) + 
      geom_point() + 
      theme(legend.title = element_blank()) + 
      stat_ellipse(level = 0.9) +
      scale_color_igv() + 
      theme_bw() +
      xlab(paste0(input$x, " (log base 10 ppm)")) +
      ylab(paste0(input$y, " (log base 10 ppm)"))
    
    ggplotly(p) %>%
      layout(height = input$plotHeight, autosize = TRUE, 
             legend = list(font = list(size = 12))) 
  })
  
}

shinyApp(ui_sample, server_sample)


# Assess membership probabilities of the outgroup samples
# Out-groups 1, 2, and 5 are large enough to be assessed for Mahalanobis distance probabilities
table(sample_pca[["pca_aug"]][[1]]$Kmeans_5)

kmeans125_pcs <- sample_pca[["pca_aug"]][[1]] %>%
                    filter(Kmeans_5 == 1 | Kmeans_5 == 2 | Kmeans_5 == 5) %>%
                    select(.fittedPC1:.fittedPC12)

kmeans125_samps <- sample_pca[["pca_aug"]][[1]] %>%
                     filter(Kmeans_5 == 1 | Kmeans_5 == 2 | Kmeans_5 == 5) %>%
                     select(Kmeans_5)

group.mem.probs(kmeans125_pcs, kmeans125_samps$Kmeans_5, 
               unique(kmeans125_samps$Kmeans_5))

################################## Core Group Structure #########################################





