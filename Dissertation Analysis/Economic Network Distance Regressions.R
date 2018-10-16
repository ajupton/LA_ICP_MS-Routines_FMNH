library(infer)
library(tidyverse)
library(igraph)
library(reshape2)
library(stringr)
library(cowplot)
library(broom)

# Read in finalized, undirected plate BR edgelist 
BReco_el_un <- read_csv("Econet_BR_UNDIRECTED_edgelist_complete_.csv")

# Read in finalized, undirected pre-migration BR edgelist
BReco_el_un_pre <- read_csv("Econet_BR_UNDIRECTED_edgelist_pre-migration_.csv")

# Read in finalized, undirected post-migration BR edgelist
BReco_el_un_post <- read_csv("Econet_BR_UNDIRECTED_edgelist_post-migration_.csv")

# Inference testing with linear models
# Take 100 samples of half the network size each from the economic BR data sets
# The idea is to explore regression trends on the slope coefficient using samples 
# from each data set. Does the trend with the entire data hold true when 
# sub-samples are taken from the data?
# This is a two-tailed test to see if a linear relationship (positive or negative) exists
# between distance (explanatory variable) and weight (response variable)
BRecopresamples <- rep_sample_n(BReco_el_un_pre[, c(3, 8)], size = 21, reps = 100)
BRecopostsamples <- rep_sample_n(BReco_el_un_post[, c(3, 8)], size = 5, reps = 100)
BRecoallsamples <- rep_sample_n(BReco_el_un[, c(3, 8)], size = 26, reps = 100)

# Add replicate col to align observed trends with random samples
ecopre_observed <- BReco_el_un_pre[, c(3, 8)] %>%
                      mutate(replicate = 200) 

ecopost_observed <- BReco_el_un_post[, c(3, 8)] %>%
                      mutate(replicate = 200) 

ecoall_observed <- BReco_el_un[, c(3, 8)] %>%
                      mutate(replicate = 200) 

# Model showing proportional similarity across time
BReco_lm_all <- ggplot(BRecoallsamples, aes(x = Distance, y = weight, group = replicate)) + 
                  geom_point(size = 2, shape = 20) + 
                  stat_smooth(geom = "line", se = FALSE, alpha = 0.4, method = "lm") + 
                  ggtitle("Ceramic Industry Economic Network Across Time") + 
                  background_grid(major = 'y', minor = "none") +
                  xlab("Distance (km)") +
                  ylab("Degree of Proportional Similarity in Compositional Groups") +
                  theme(strip.background = element_blank(),
                        strip.text.x = element_blank()) +
                  stat_smooth(data = ecoall_observed, aes(x = Distance, y = weight), color ="red3", 
                              linetype = "twodash", method = "lm", se = FALSE) 

# Model showing proportional similarity in the pre-migration CIRV
BReco_lm_pre <- ggplot(BRecopresamples, aes(x = Distance, y = weight, group = replicate)) + 
                  geom_point(size = 2, shape = 20) + 
                  stat_smooth(geom = "line", se = FALSE, alpha = 0.4, method = "lm") + 
                  ggtitle("Pre-Migration Ceramic Industry Economic Network") + 
                  background_grid(major = 'y', minor = "none") +
                  xlab("Distance (km)") +
                  ylab("Degree of Proportional Similarity in Compositional Groups") +
                  theme(strip.background = element_blank(),
                        strip.text.x = element_blank()) +
                  stat_smooth(data = ecopre_observed, aes(x = Distance, y = weight), color ="red3", 
                              linetype = "twodash", method = "lm", se = FALSE) 

# Model showing proportional similarity in the post-migration CIRV
BReco_lm_post <- ggplot(BRecopostsamples, aes(x = Distance, y = weight, group = replicate)) + 
                  geom_point(size = 2, shape = 20) + 
                  stat_smooth(geom = "line", se = FALSE, alpha = 0.4, method = "lm") + 
                  ggtitle("Post-Migration Ceramic Industry Economic Network") + 
                  background_grid(major = 'y', minor = "none") +
                  xlab("Distance (km)") +
                  ylab("Degree of Proportional Similarity in Compositional Groups") +
                  theme(strip.background = element_blank(),
                        strip.text.x = element_blank()) +
                  stat_smooth(data = ecopost_observed, aes(x = Distance, y = weight), color ="red3", 
                              linetype = "twodash", method = "lm", se = FALSE) 

# Inference
# First, let's calculate the observed slope of the lm in the jar and plate attribute networks
BReco_all_slope <- lm(weight ~ Distance, data = BReco_el_un) %>%
                        tidy() %>%   
                        filter(term == "Distance") %>%
                        pull(estimate)    

BReco_pre_slope <- lm(weight ~ Distance, data = BReco_el_un_pre) %>%
                        tidy() %>%   
                        filter(term == "Distance") %>%
                        pull(estimate)    

BReco_post_slope <- lm(weight ~ Distance, data = BReco_el_un_post) %>%
                        tidy() %>%   
                        filter(term == "Distance") %>%
                        pull(estimate)    

# Simulate 500 slopes with a permuted dataset for economic network - this will allow us to 
# develop a sampling distribution of the slop under the hypothesis that there is no 
# relationship between the explanatory (Distance) and response (weight) variables. 
set.seed(1568)
BReco_all_perm_slope <- BReco_el_un %>%
                          specify(weight ~ Distance) %>%
                          hypothesize(null = "independence") %>%
                          generate(reps = 500, type = "permute") %>%
                          calculate(stat = "slope") 

BReco_pre_perm_slope <- BReco_el_un_pre %>%
                          specify(weight ~ Distance) %>%
                          hypothesize(null = "independence") %>%
                          generate(reps = 500, type = "permute") %>%
                          calculate(stat = "slope") 

BReco_post_perm_slope <- BReco_el_un_post %>%
                          specify(weight ~ Distance) %>%
                          hypothesize(null = "independence") %>%
                          generate(reps = 500, type = "permute") %>%
                          calculate(stat = "slope") 

ggplot(BReco_all_perm_slope, aes(x = stat)) + geom_density() + theme_classic()
ggplot(BReco_pre_perm_slope, aes(x = stat)) + geom_density() + theme_classic()
ggplot(BReco_post_perm_slope, aes(x = stat)) + geom_density() + theme_classic()

mean(BReco_all_perm_slope$stat)
mean(BReco_pre_perm_slope$stat)
mean(BReco_post_perm_slope$stat)
sd(BReco_all_perm_slope$stat)
sd(BReco_pre_perm_slope$stat)
sd(BReco_post_perm_slope$stat)

# Calculate the absolute value of the slope
abs_BRco_all_obs_slope <- lm(weight ~ Distance, data = BReco_el_un) %>%
                            tidy() %>%   
                            filter(term == "Distance") %>%
                            pull(estimate) %>%
                            abs()

abs_BReco_pre_obs_slope <- lm(weight ~ Distance, data = BReco_el_un_pre) %>%
                              tidy() %>%   
                              filter(term == "Distance") %>%
                              pull(estimate) %>%
                              abs()

abs_BReco_post_obs_slope <- lm(weight ~ Distance, data = BReco_el_un_post) %>%
                              tidy() %>%   
                              filter(term == "Distance") %>%
                              pull(estimate) %>%
                              abs()

# Compute the p-value  
BReco_all_perm_slope %>% 
  mutate(abs_BReco_all_obs_slope = abs(stat)) %>%
  summarize(p_value = mean(abs_BReco_all_obs_slope > BReco_all_perm_slope))

BReco_pre_perm_slope %>% 
  mutate(abs_BReco_pre_obs_slope = abs(stat)) %>%
  summarize(p_value = mean(abs_BReco_pre_obs_slope > BReco_pre_perm_slope))

BReco_post_perm_slope %>% 
  mutate(abs_BReco_post_obs_slope = abs(stat)) %>%
  summarize(p_value = mean(abs_BReco_post_obs_slope > BReco_post_perm_slope))

# Linear models sans visualization
# explore residuals 
BReco_all_lm <- augment(lm(weight ~ Distance, data = BReco_el_un))
BReco_pre_lm <- augment(lm(weight ~ Distance, data = BReco_el_un_pre))
BReco_post_lm <- augment(lm(weight ~ Distance, data = BReco_el_un_post))

# Check SSE - how well do the models fit?
augment(lm(weight ~ 1, data = BReco_el_un)) %>% summarize(SSE = var(.resid)) # null
BReco_all_lm %>% summarize(SSE = var(.resid))

augment(lm(weight ~ 1, data = BReco_el_un_pre)) %>% summarize(SSE = var(.resid)) # null
BReco_pre_lm %>% summarize(SSE = var(.resid))

augment(lm(weight ~ 1, data = BReco_el_un_post)) %>% summarize(SSE = var(.resid)) # null
BReco_post_lm %>% summarize(SSE = var(.resid))

# Looks like the models do fit very well

# Breakdown of linear model results for plate attribute networks 
summary(lm(weight ~ Distance, data = BReco_el_un)) # for each 1 km increase in distance, weight drops 0.0007723 and at 0 distance, a weight of 0.7378 is expected
summary(lm(weight ~ Distance, data = BReco_el_un))$coefficients # all = p-value of 0.1454, fail to  
# reject null hypothesis - no significant linear relationship b/t distance and weight across time

summary(lm(weight ~ Distance, data = BReco_el_un_pre)) # for each 1 km increase in distance, weight drops 0.0003959 and at 0 distance, a weight of 0.7385 is expected
summary(lm(weight ~ Distance, data = BReco_el_un_pre))$coefficients # pre p-value of 0.6918, fail to
#reject the null hypothesis - no significant linear relationship b/t distance and weight in pre

summary(lm(weight ~ Distance, data = BReco_el_un_post)) # for each 1 km increase in distance, weight drops 0.0004263 and at 0 distance, a weight of 0.6835 is expected
summary(lm(weight ~ Distance, data = BReco_el_un_post))$coefficients # post p-value of 0.5499, fail 
# to reject null - no significant linear relationship b/t distance and weight in post

# Check correlations
cor(BReco_el_un$Distance, BReco_el_un$weight)
cor(BReco_el_un_pre$Distance, BReco_el_un_pre$weight)
cor(BReco_el_un_post$Distance, BReco_el_un_post$weight)


### No relationship between distance and degree of economic interactions is able to be identified
#   this is interesting, as it would be expected that sites closer in proximity would exhibit 
#   stronger economic relationships via a higher degree of exchange of finished vessels, overlapping
#   resource exploitation areas, or similar paste preparation and production regimes. 