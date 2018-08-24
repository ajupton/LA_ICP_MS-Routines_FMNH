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

# Read in full dataset
all_samples <- read_csv("Upton_results_samples_shell_corrected_August_21_2018.csv")

# Read in clay context data 
clay_context <- read_csv("clay context data.csv")

# Clean up clay and clay context sample names
clay_context$Sample <- str_replace(clay_context$Sample, pattern = "_1" %R% END, "")
clay$Sample <- str_replace(clay$Sample, pattern = "_1" %R% END, "")

# Each clay sample is named "C##_1"
# An expedient way of isolating the clay samples
clay <- arrange(all_samples[str_detect(all_samples$Sample, "C[:digit:]"), ], Sample)

# Join clay data to clay context

# Take the log of the elemental composition data since they are on very different 
# scales
claylog <- log(clay[,6:ncol(clay)])
claylog <- bind_cols(clay[, 1:5], claylog)
names(clay)

# Remove problem elements following Golitko (2010) 
problem_elements <- c("P", "Sr", "Ba", "Ca", "Hf", "As", "Cl")
names.use <- names(claylog)[!(names(claylog) %in% problem_elements)]
claylog_good <- claylog[, names.use]
names(claylog_good)

# Check to ensure the elements were removed
anti_join(data.frame(names(clay)), 
          data.frame(names(claylog_good)), by = c("names.clay." = "names.claylog_good."))

# Remove one non-clay sample
claylog_good <- filter(claylog_good, Sample != "C26_1")

# Exploring PCA
clay_pca <- claylog_good %>% 
  nest() %>% 
  mutate(pca = map(data, ~ prcomp(.x %>% select(Si:Th), 
                                  center = TRUE, scale = TRUE)),
         pca_aug = map2(pca, data, ~augment(.x, data = .y)))

# Check variance explained by each model
var_exp <- clay_pca %>% 
  unnest(pca_aug) %>% 
  summarize_at(.vars = vars(contains("PC")), .funs = funs(var)) %>% 
  gather(key = pc, value = variance) %>% 
  mutate(var_exp = variance/sum(variance),
         cum_var_exp = cumsum(var_exp),
         pc = str_replace(pc, ".fitted", ""))

# Looks like we need to retain the first 9 PC's to hit 90% of the data's variability
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
  labs(y = "Variance",
       title = "Variance explained by each principal component")


# Plot PCs against each other
cp1p2_plot <- clay_pca %>%
        mutate(
          pca_graph = map2(
            .x = pca,
            .y = data,
            ~ autoplot(.x, loadings = TRUE, loadings.label = TRUE,
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


# Shiny app to biplot the various elements against one another
# With 46 elements, there are p(p-1)/2 or 1035 biplots to investigate!
# Therefore, it's a lot easier to make an app to easily and quickly run through the options

############
##   UI   ##
############
ui <- fluidPage(
pageWithSidebar (
  headerPanel('Bivariate Plotting'),
  sidebarPanel(
    selectInput('x', 'X Variable', names(claylog_good), 
                selected = names(claylog_good)[[7]]),
    selectInput('y', 'Y Variable', names(claylog_good),
                selected = names(claylog_good)[[6]]),
    selectInput('color', 'Color', names(claylog_good)),
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
    p <- ggplot(data = claylog, aes_string(x = input$x, 
                                          y = input$y, 
                                          color = input$color, 
                                          shape = input$color)) + 
      geom_point() + 
      theme(legend.title = element_blank()) + 
      stat_ellipse() +
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
