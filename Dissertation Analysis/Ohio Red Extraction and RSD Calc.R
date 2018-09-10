# Extract Ohio Reds and Calculate Relative Standard Deviation

library(tidyverse)
library(stringr)
library(plotly)

# Import data
dfall <- read_csv("Upton_results_OhioRed_August_21_2018.csv")

# Change Sample column to all lower case to ensure complete string detection
dfall$Sample <- tolower(dfall$Sample)

# Search the sample column for the word Ohio based on the abbreviation oh
ohio <- str_detect(dfall$Sample, "oh")

# Double check by searching same column for Red
red <- str_detect(dfall$Sample, "red")

# Check to see if the two detection methods are identical
sum(ohio == red) == nrow(dfall)

# Index to extract all Ohio Red Samples
ohioreds <- dfall[red,]

# Function to calculate RSD
RSD <- function(x){
  meann <- mean(x)
  relsd <- sd(x)/meann
  relsd
}

# Calculate RSD across the rows
redRSD <- sapply(ohioreds[, 3:ncol(ohioreds)], RSD)

# Calculate average and standard deviation of values across each of the Ohio Reds
redAVG <- ohioreds %>%
              gather(element, sample, Si:Th) %>%
              group_by(element) %>%
              summarize(Avg = mean(sample), SD = sd(sample))

# Plot to check average values
ohioreds %>%
  gather(element, sample, Si:Th) %>%
  group_by(element, Date) %>%
  summarize(AVG = mean(sample)) %>%
ggplot(aes(x = Date)) + 
  geom_line(aes(y = AVG, color = element, group = element)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Filter the plot to look at HREE and LREE average values
p <- ohioreds %>%
      gather(element, sample, Si:Th) %>%
      group_by(Date, element) %>%
      summarize(AVG = mean(sample)) %>%
      filter(AVG < 100) %>%
      ggplot(aes(x = Date)) + 
      geom_line(aes(y = AVG, color = element, group = element)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplotly(p)

# Check very high RSD elements
all_samples %>%
  gather(element, sample, Si:Th) %>%
  group_by(Date, element) %>%
  summarize(AVG = mean(sample)) %>%
  filter(element == "Bi") %>%
  ggplot(aes(x = Date)) + 
  geom_line(aes(y = AVG, color = element, group = element)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Bind Ohio Red averages to relative standard deviations
redRSD <- data.frame(redRSD)
redRSD <- rownames_to_column(redRSD, var = "element")
OHred_avg_rsd <- left_join(redRSD, redAVG, by = "element")

# Add RSD to the Ohio Red Samples
red_with_RSD <- bind_rows(ohioreds, redRSD)

write_csv(OHred_avg_rsd, "Ohio Red Averages and RSD_Aug-21-2018.csv")
write_csv(red_with_RSD, "Ohio Reds with RSD_Aug-21-2018.csv")

# Now check for any differences between samples run on different machines
ohioreds$Date <- as.POSIXct(paste(ohioreds$Date), format = "%Y-%b-%d", tz = "UTC") 
redAVG_group <- ohioreds %>%
                  mutate(Machine = ifelse(Date > as.POSIXct('2016-01-01', tz = "UTC"), 
                                          "New", "Old")) %>%
                  gather(element, sample, SiO2:Th) %>%
                  group_by(Machine, element) %>%
                  summarize(Avg = mean(sample), SD = sd(sample))

# Count the number of Ohio Red samples run on each machine
ohioreds %>%
  mutate(Machine = ifelse(Date > as.POSIXct('2016-01-01', tz = "UTC"), 
                          "New", "Old")) %>%
  select(Machine) %>%
  group_by(Machine) %>%
  summarize(num = n())

write_csv(redAVG_group, "Ohio Reds across Machines_Aug_22_2018.csv")

