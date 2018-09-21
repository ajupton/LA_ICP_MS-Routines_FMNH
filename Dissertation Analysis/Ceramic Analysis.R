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
summary(lm(Ni ~ ., data = samples[,3:length(samples)]))

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

# Plot the first two PCs with Geography_2 as group separation
sample_pca %>%
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

# This shows significant overlap but a general trend that follows the clay: in general there is
# less elemental enrichment in clay resources in the southern portion of the CIRV

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

# Interact with the chart above
#ggplotly(site_pc1pc2[[1]]) # plotly drops the stat_ellipse frames for some reason

# Based on the initial inspection of PCs 1 and 2, it looks like three elements in particular 
# are driving some of the group separation (subtle as it is): Mo, Mn, and Si
# Let's plot those three elements in a 3D scatter plot
Mo_Mn_Si <- plot_ly(sample_new_pcaready, x = ~Mo, y = ~Mn, z = ~Si, color = ~Geography_2)

# Explore Samples by date
ggplotly(ggplot(sample_new_pcaready, aes(x = Mo, y = Mg, color = Date)) + 
  stat_ellipse(aes(color = Geography_2)) + geom_text(aes(label = Date), size = 2))


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


