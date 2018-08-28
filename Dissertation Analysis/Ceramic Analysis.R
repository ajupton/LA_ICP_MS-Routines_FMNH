# Checking for linear relationships with calcium

library(tidyverse)
library(infer)
library(broom)
library(stringr)
library(plotly)
library(rebus)
library(xlsx)

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
                mutate(id = parse_number(clay$Sample)) 

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
p <- ggplot(samples, aes(x = Sr, y = CaO)) + geom_smooth() + geom_point()
#ggplotly(p)





