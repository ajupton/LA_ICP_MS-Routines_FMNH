# Clay analysis

library(tidyverse)
library(stringr)
library(plotly)
library(shiny)
library(shinydashboard)
library(ggsci)
library(broom)
library(knitr)
library(ggfortify)
library(stats)
library(ICSNP)
library(factoextra)
library(dendextend)

# Read in full dataset
all_samples <- read_csv("Upton_results_samples_shell_corrected_August_21_2018.csv")

# Read in clay context data 
clay_context <- read_csv("clay context data.csv")

# Each clay sample is named "C##_1"
# An expedient way of isolating the clay samples
clay <- arrange(all_samples[str_detect(all_samples$Sample, "C[:digit:]"), ], Sample)

# Clean up clay and clay context sample names
clay_context$Sample <- str_replace(clay_context$Sample, pattern = "_1" %R% END, "")
clay$Sample <- str_replace(clay$Sample, pattern = "_1" %R% END, "")

# Join clay data to clay context
clay <- left_join(clay_context, clay)

# Take the log of the elemental composition data since they are on very different 
# scales
claylog <- log10(clay[,8:ncol(clay)])
claylog <- bind_cols(clay[, 1:7], claylog)

# Number of clay samples analyzed 
claylog %>%
  summarise(num = n())

# Remove problem elements. Some elements are known to be unreliably measured using the ICP-MS
# at the EAF. Following Golitko (2010), these include the following elements. 
problem_elements <- c("P", "Sr", "Ba", "Ca", "Hf", "As", "Cl")

# Other elements such as Ca and Sr are affected by shell tempering. Want to drop those as well. 

# Overall these are the Elements retained - 44 in all.
elems_retained <- c("Al","B", "Be", "Ce", "Co", "Cr", "Cs", "Dy", "Er", "Eu", "FeO",
                    "Gd", "Ho", "In", "K", "La", "Li", "Lu", "Mg", "MnO", "Mo", "Na", "Nb", 
                    "Nd", "Ni", "Pb", "Pr", "Rb", "Sc", "Si", "Sm", "Sn", "Ta", "Tb", "Th", "Ti",
                    "Tm", "U", "V", "W", "Y", "Yb", "Zn", "Zr")

names.use <- names(claylog)[(names(claylog) %in% elems_retained)]
# length(names.use) == length(elems_retained) # check that all elements are retained 
claylog_good <- claylog[, names.use]

# Check to ensure the elements were removed are supposed to be removed
anti_join(data.frame(names(clay)), 
          data.frame(names(claylog_good)), by = c("names.clay." = "names.claylog_good."))

# Need to drop the "O" for oxide after elements measured as %oxide composition since they have 
# already been converted to ppm
names(claylog_good) <- c("Si","Na","Mg","Al","K","Mn","Fe","Ti","Li", "Be","B","Sc","V","Cr","Ni",
                         "Co","Zn","Rb","Zr","Nb","In","Sn","Cs","La","Ce","Pr","Ta","Y","Pb","U",
                         "W","Mo","Nd","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Th") 

# Bind sample id and other data to the logged chemical concentrations 
clay_pcaready <- bind_cols(claylog[,c(1:7)], claylog_good)

# Remove two non-clay sample
clay_pcaready <- filter(clay_pcaready, Sample != "C26") %>% filter(Sample != "C31")

#write_csv(clay_pcaready, "Clay PCA Ready.csv")

# Exploring PCA
clay_pca <- clay_pcaready %>% 
  nest() %>% 
  mutate(pca = map(data, ~ prcomp(.x %>% select(Si:Th))),
         pca_aug = map2(pca, data, ~augment(.x, data = .y)))

# Check variance explained by each model
var_exp <- clay_pca %>% 
  unnest(pca_aug) %>% 
  summarize_at(.vars = vars(contains("PC")), .funs = funs(var)) %>% 
  gather(key = pc, value = variance) %>% 
  mutate(var_exp = variance/sum(variance),
         cum_var_exp = cumsum(var_exp),
         pc = str_replace(pc, ".fitted", ""))

# Looks like we need to retain the first 7 PC's to hit 90% of the data's variability
# Graphing this out might help
var_exp %>% 
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


# Plot PCs 1 & 2 against each other
cp1p2_plot <- clay_pca %>%
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
                   title = "First two principal components of PCA on CIRV Clay dataset")
          )
        ) %>%
      pull(pca_graph)

# autoplot is lazy with color. In order to make this publication friendly, have to 
# manually edit the color scales
cp1p2_plot[[1]] + scale_fill_manual(values = c("black","black")) + 
  scale_color_manual(values = c("black","black","black")) 

# Plot PCs 1 & 3 against each other
cp1p3_plot <- clay_pca %>%
  mutate(
    pca_graph = map2(
      .x = pca,
      .y = data,
      ~ autoplot(.x, x = 1, y = 3, loadings = TRUE, loadings.label = TRUE,
                 loadings.label.repel = TRUE,
                 loadings.label.colour = "black",
                 loadings.colour = "gray85",
                 loadings.label.alpha = 0.5,
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
             y = "Principal Component 3",
             title = "First two principal components of PCA on CIRV Clay dataset")
    )
  ) %>%
  pull(pca_graph)

# autoplot is lazy with color. In order to make this publication friendly, have to 
# manually edit the color scales
cp1p3_plot[[1]] + scale_fill_manual(values = c("black","black")) + 
  scale_color_manual(values = c("black","black","black")) 

# Shiny app to biplot the various elements against one another
# With 44 elements, there are p(p-1)/2 or 946 biplots to investigate!
# Therefore, it's a lot easier to make an app to easily and quickly run through the options

############
##   UI   ##
############
ui <- fluidPage(
pageWithSidebar (
  headerPanel('Bivariate Plotting'),
  sidebarPanel(
    selectInput('x', 'X Variable', names(clay_pcaready), 
                selected = names(clay_pcaready)[[8]]),
    selectInput('y', 'Y Variable', names(claylog_good),
                selected = names(clay_pcaready)[[9]]),
    selectInput('color', 'Color', names(clay_pcaready)),
    #Slider for plot height
    sliderInput('plotHeight', 'Height of plot (in pixels)', 
                min = 100, max = 2000, value = 550)
  ),
  mainPanel(
    plotlyOutput('plot1')
  )
)
)

############
## Server ##
############


server <- function(input, output, session) {
  
  # Combine the selected variables into a new data frame
  selectedData <- reactive({
    claylog_good[, c(input$x, input$y)]
  })

  
  output$plot1 <- renderPlotly({
    
    #Build plot with ggplot syntax 
    p <- ggplot(data = clay_pcaready, aes_string(x = input$x, 
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

shinyApp(ui, server)

# Based on the biplots, it looks like there is good separation for the most part 
# in the north and south portions of the valley when comparing Lithium to Vanadium or Beryllium

# Biplot of Li and V
ggplot(data = clay_pcaready, aes(x = Li, y = V, color = Geography_2, shape = Geography_2)) + 
  geom_point() + 
  theme(legend.title = element_blank()) + 
  stat_ellipse(level = 0.9) +
  scale_color_igv() + 
  theme_bw() +
  xlab("Lithium (log base 10 ppm)") +
  ylab("Vanadium (log base 10 ppm)") +
  scale_fill_manual(values = c("black","black")) + 
  scale_color_manual(values = c("black","black","black")) 

# Biplot of Li and Be
ggplot(data = clay_pcaready, aes(x = Li, y = Be, color = Geography_2, shape = Geography_2)) + 
  geom_point() + 
  theme(legend.title = element_blank()) + 
  stat_ellipse(level = 0.9) +
  scale_color_igv() + 
  theme_bw() +
  xlab("Lithium (log base 10 ppm)") +
  ylab("Beryllium (log base 10 ppm)") +
  scale_fill_manual(values = c("black","black")) + 
  scale_color_manual(values = c("black","black","black")) 

# Biplot of Li and Be
ggplot(data = clay_pcaready, aes(x = Li, y = Be, color = Geography_2, shape = Geography_2)) + 
  geom_point() + 
  theme(legend.title = element_blank()) + 
  stat_ellipse(level = 0.9) +
  scale_color_igv() + 
  theme_bw() +
  xlab("Lithium (log base 10 ppm)") +
  ylab("Beryllium (log base 10 ppm)") +
  scale_fill_manual(values = c("black","black")) + 
  scale_color_manual(values = c("black","black","black")) 

# Biplot of Ni and Cs
ggplot(data = clay_pcaready, aes(x = Ni, y = Cs, color = Geography_2, shape = Geography_2)) + 
  geom_point() + 
  theme(legend.title = element_blank()) + 
  stat_ellipse(level = 0.9) +
  scale_color_igv() + 
  theme_bw() +
  xlab("Nickel (log base 10 ppm)") +
  ylab("Cesium (log base 10 ppm)") +
  scale_fill_manual(values = c("black","black")) + 
  scale_color_manual(values = c("black","black","black")) 

# A table of average element concentrations and standard deviations between the two 
# groups may be instructive of their differences numerically as opposed to visually
clay_group_ave_std <- clay_pcaready %>%
                        select(Geography_2, Si:Th) %>%
                        gather(Element, Si:Th, -Geography_2) %>%
                        mutate(`Si:Th` = 10^`Si:Th`) %>% # convert from log 10
                        group_by(Geography_2, Element) %>%
                        summarize(mean = mean(`Si:Th`, na.rm = TRUE), std = sd(`Si:Th`, na.rm = TRUE))

# Count number of clay in the different groups
clay_pcaready %>%
  group_by(Geography_2) %>%
  summarize(count = n())

write_csv(clay_group_ave_std, "Clay group averages and stds.csv")

# Almost every element is enriched in northerly clays and as a result depleted in the 
# southerly clays, taking a look at that via a histogram is instructive
clay_pcaready %>%
  select(Geography_2, Si:Th) %>%
  gather(Element, Si:Th, -Geography_2) %>%
  mutate(`Si:Th` = 10^`Si:Th`) %>% 
  filter(Element == "Sn") %>%
  ggplot(aes(x = Element, group = Geography_2, y = `Si:Th`)) + geom_boxplot()
  
# It looks like there is a good deal of separation in the geochemistry of clays between the 
# Northern portion of the central Illinois River Valley (including the Spoon/Illinois confluence) 
# and the Southern portion of the CIRV, south of the Spoon River
# But let's check to see if statistical techniques come to a similar conclusion

###____________________________________HCA________________________________________________###
# First, we create a data frame for distance calculations including the elemental data only
clay_for_dist <- claylog_good
rownames(clay_for_dist) <- claylog$Sample

# Now let's perform some hierarchical clustering using Euclidean distance
clay_hca <- hclust(dist(clay_for_dist))

# Create dendrogram object
dend_clay <- as.dendrogram(clay_hca)

# Plot dendogram object to look for good cut-off heights - 2.5 seems to be a good height
plot(dend_clay, nodePar = list(lab.cex = .75, pch = NA))

# Looks like the hierarchical clustering doesn't group precisely as the geographic/geologic
# prior knowledge would suggest. This is an indication of the hetergeneous nature of clay as 
# well as the complex geological processes that have resulted in clay availability in the CIRV.


###____________________________Mahalanobis Distance________________________________________###
# Since HCA wasn't overly insightful, we can at least check membership probabilities between
# the north and south groups statistically. The standard method of doing this in 
# archaeology is via Mahalanobis distance, which is commonly used for outlier detection. 


# Extract the first 7 PC's (accounting for 90% of variability) and bind to 
# sample/geography data
clay_pc1to7 <- clay_pca %>% 
                  unnest(pca_aug) %>% 
                  select(starts_with(".fitted")) %>%
                  bind_cols(clay_pcaready[, c(1,3)], .) %>%
                  select(c(1:9))

clay1to7_north <- clay_pc1to7 %>% filter(Geography_2 == "North")

clay1to7_south <- clay_pc1to7 %>% filter(Geography_2 == "South")
      
# Edit colnames
colnames(clay_pc1to7) <- str_remove(colnames(clay_pc1to7), ".fitted")

# Mahalanobis distance of North to North
mahalanobis(clay1to7_north[,3:9], colMeans(clay1to7_north[,3:9]), cov(clay1to7_north[,3:9]))

# With 7 predictor variables (PCs 1-7), the critical chi-square value is 24.32
# Given that the highest MD value among the northerly clays is 20.92, it doesn't
# look like there are any outliers

# Have to pair down the number of predictors to 5 for the South, since there are only 7 
# samples. The critical chi-square value is 20.52 for that many, looking good for the south.
mahalanobis(clay1to7_south[,3:7], colMeans(clay1to7_south[,3:7]), cov(clay1to7_south[,3:7]))

# Let's now look at group membership probabilities. This function written by Matt Peeples allows for 
# for calculating group membership probabilities by chemical compositional distance using 
# Mahalanobis distances and Hotellings T^2 statistic
group.mem.probs <- function(x2.l,attr1.grp,grps) {
  
  # x2.l = transformed element data
  # attr1 = group designation by sample
  # grps <- vector of groups to evaluate
  
  probs <- list()
  for (m in 1:length(grps)) {
    x <- x2.l[which(attr1.grp==grps[m]),]
    probs[[m]] <- matrix(0,nrow(x),length(grps))
    colnames(probs[[m]]) <- grps
    rownames(probs[[m]]) <- rownames(x)
    
    grps2 <- grps[-m]
    
    p.val <- NULL
    for (i in 1:nrow(x)) {p.val[i] <- HotellingsT2(x[i,],x[-i,])$p.value}
    probs[[m]][,m] <- round(p.val,5)*100
    
    for (j in 1:length(grps2)) {
      p.val2 <- NULL
      for (i in 1:nrow(x)) {p.val2[i] <- HotellingsT2(x[i,],x2.l[which(attr1.grp==grps2[j]),])$p.value}
      probs[[m]][,which(grps==grps2[j])] <- round(p.val2,5)*100}}
  return(probs)
}

# But how do the samples compare to each other on the first 5 PCs (85% ov observed variability)?
group.mem.probs(clay_pc1to7[3:5], clay_pc1to7$Geography_2, unique(clay_pc1to7$Geography_2)) 

# How about using some elements that show good separation between the groups?
group.mem.probs(clay_pcaready[, c("Ni", "Cs")], clay_pc1to7$Geography_2, unique(clay_pc1to7$Geography_2)) 

# In both cases, there is a marked lack of clear group separation in statistical space for 
# samples in both groups. That is, there are samples defined as North that have a higher 
# probability of grouping with the Southerly sherds and vice versa. 
# To a certain degree, this is expected - this is an experimental analysis looking within a 
# single river valley, and indeed there is not statistically significant separation between 
# the groups as a result. 
# Nevertheless, it is instructive that chemical differences do appear as one moves from the 
# northeast to the southwest in the CIRV, conforming to geologic patterns of exposing parent 
# material of older ages. As a result, an argument can be made that pottery would likely 
# follow this patterning based on raw material availability. 


## Exploratory cluster analysis

# Optimal number of clusters based on the elbow method using the total within sum of squares
fviz_nbclust(clay_pc1to7[3:9], kmeans, method = "wss")

clay_dist <- hclust(dist(clay_pc1to7[3:9]))

View(clay_pc1to7)

# Create dendrogram object
clay_dend_df_com <- as.dendrogram(clay_dist)

# Plot dendogram object to look for good cut-off heights - 2.5 seems to be a good height
plot(clay_dend_df_com, nodePar = list(lab.cex = 0.15, pch = NA))

dend_2.5 <- color_branches(clay_dend_df_com, h = 1.950)
plot(dend_2.5, cex.axis = 0.75, cex.lab = 0.75, nodePar = list(lab.cex = .85, pch = NA))
