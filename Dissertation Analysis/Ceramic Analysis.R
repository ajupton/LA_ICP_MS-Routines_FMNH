# Checking for linear relationships with calcium

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

###_______________ Add features to ceramic samples ________________________###
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

########___________End of data cleaning, beginnging of statistical analysis__________########

# First step is to take the log base 10 of all samples to account for scalar differences
# in the magnitude of chemical compositions across the elements, from major to minor to trace
sample_new_pcaready[,13:56] <- log10(sample_new_pcaready[,13:56])

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
                             frame = FALSE,
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


geo2_pc1pc2[[1]]

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

###____________________________________Cluster Analysis_____________________________________###

#####################################
### Hierarchical Cluster Analysis ###
#####################################
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

################################
### K-means Cluster Analysis ###
################################

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


##################################
### K-mediods Cluster Analysis ###
##################################

# For k-mediods, we'll be using the pam function from the cluster package. pam stands for 
# "partitioning acount mediods"

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



####_________________Beging Mahalanobis distance and membership assignments________________###

# Let's now look at group membership probabilities. This function written by Matt Peeples allows for 
# for calculating group membership probabilities by chemical compositional distance using 
# Mahalanobis distances and Hotellings T^2 statistic
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


