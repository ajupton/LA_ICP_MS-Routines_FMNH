## Reading Geochemical data in an Excel Spreadsheet into R

# tidyverse and readxl
library(tidyverse)
library(readxl)
library(stringr)

# The routine below will read LA-ICP-MS data stored in .xlsx files into R
# It is expected that each day of LA-ICP-MS analysis is stored as a separate sheet 
# within the .xlsx workbook and each sheet is named the day of analysis. 

# This routine is structured for ceramic sample data, but with some minor adjustments, 
# can be applied to other sample types based on the standards used 
# (i.e. replace the search for "red" to the name of the standard used in line 90 below)

# Determine path to file 
path <- "~file_name.xlsx"

# Use map to iterate read_excel over each worksheet in the workbook
ld <- path %>%
        excel_sheets() %>%
        #set_names() %>%   
        set_names(.,.) %>% # Had error messages from the above line, this fixes the issues
        map(read_excel, path = path)

# Bind the columns in the lists together to form one dataframe
df <- bind_cols(ld)

# Change name of first column to element
names(df)[names(df) == 'X__1'] <- 'element'

# Now, let's get tidy!
# Grab the first column as rownames, which will become the variable names
rnames <- df[,1]

# Then grab the column names, which will become a new column "Sample" once transposed
Sample <- colnames(df[-1]) # drop the first name because it will become the rownames

# Transpose the dataframe
df <- t(df[, -1]) # drop the first column or it will convert the numbers to strings

# Set the column names
colnames(df) <- unlist(rnames) #rnames is stored as a list, so we have to unlist it 

# Convert to tibble dataframe
df <- tbl_df(df)

# Let's add the date as a column to our data frame so we know when each sample was run
# First we need to figure out how many samples were run each day
ld_lengths <- lapply(ld, length) # could also use map(ld, length)

# With that information we can create a simple for loop to replicate the 
# dates the appropriate number of times for the number of samples run each day
res1 <- as.data.frame(NULL)
for(i in names(ld_lengths)) {
  res <- rep(i, ld_lengths[[i]])
  res1 <- c(res1, res)
}

# Now we'll create a data frame of those dates and add it to our sample data
date_col <- tbl_df(sapply(res1, paste0, collapse = ""))
colnames(date_col) <- "Date"
df <- cbind(date_col[2:nrow(date_col),], df)

# Add column of samples names, which were the columns names before transposing
df <- cbind(Sample, df)

# Use stringr to get rid of repetitive element row names - which have an "X" 
# in them by default since they don't have a column name
dfnames <- df$Sample
x_detect <- str_detect(dfnames, "X")
df <- df[!x_detect, ]

# One pesky column name has a note in it, let's get rid of it too 
# (this method can be used to clean up any columns that have notes/extra text)

# note <- str_detect(df$Sample, "High")
# df <- df[!note, ]

# Finally we can convert the sample data to numeric to allow for calculations
df[,3:ncol(df)] <- sapply(df[,3:ncol(df)], as.numeric)

### At this point, you can jump to the "Percent Oxide to PPM.R" routine or just use the ###
#   data as is

# Write the full dataframe to a csv
write_csv(df, "full_data_frame.csv")

# Now lets get rid of the Ohio Red Standard Samples
dfsamps <- df$Sample
orows <- str_detect(tolower(dfsamps), "red")
df_samples <- df[!orows, ]

# Write csv with samples only, Ohio Red standards removed
write_csv(df_samples, "data_samples_only.csv")
