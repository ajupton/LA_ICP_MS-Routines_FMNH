# Map of sites and clay resources

library(ggmap)
library(tidyverse)
library(ggrepel)

# Read in context data
cc_loc <- read_csv("Clay Ceramic lat long.csv")

# Set center for the map
lat_mid <- mean(cc_loc$lat)
lon_mid <- mean(cc_loc$lon)

# Get map from google map terrain without any labels (via the style argument)
b <- get_googlemap(center = c(lon = lon_mid, lat = lat_mid), zoom = 9, 
                   maptype = "terrain", source = "google", 
                   style = 'feature:all|element:labels|visibility:off')

# Create map of sites/clay samples and label to check for accuracy
ggmap(b) + geom_point(data = cc_loc, aes(x = lon, y = lat, shape = Type)) + 
  geom_text_repel(data = cc_loc[c(1:17),], aes(x = lon, y = lat, label = Site_Sample))

# Create transparent background map to overlay on geologic map
map <- ggplot() + geom_point(data = cc_loc, aes(x = lon, y = lat, shape = Type)) + 
        theme(
          panel.background = element_rect(fill = "transparent"), # bg of the panel
           plot.background = element_rect(fill = "transparent"), # bg of the plot
           panel.grid.major = element_blank(), # get rid of major grid
           panel.grid.minor = element_blank(), # get rid of minor grid
           legend.background = element_rect(fill = "transparent"), # get rid of legend bg
           legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
        ) +
        xlab("") + ylab("")

# Save map with transparent background
ggsave(map, filename = "site-clay map.png", bg = "transparent")

## After exporting, map of sites/clay samples was overlain on top of a bedrock geology map
#  of the state of Illinois

