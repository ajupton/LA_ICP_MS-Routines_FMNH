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
        results[s1,s2] <- round((1 - (sum(abs(x[s1, ] / sum(x[s1,]) - x[s2, ] / sum(x[s2,]))))/2)/divisor.final, digits=3)
      }
    } 
  } else {  
    for (s1 in 1:rd) {
      for (s2 in 1:rd) {
        results[s1,s2] <- round(1 - (sum(abs(x[s1, ] / sum(x[s1,]) - x[s2, ] / sum(x[s2,]))))/2, digits=3)
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

