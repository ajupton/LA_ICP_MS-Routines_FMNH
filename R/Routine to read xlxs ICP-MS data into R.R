## Reading Geochemical data into R

# tidyverse and readxl
library(tidyverse)
library(readxl)
library(stringr)

# The routine below will read LA-ICP-MS data stored in .xlsx files into R
# It is expected that each day of LA-ICP-MS analysis is stored as a separate sheet within
# the .xlsx workbook. 

# This routine is structured for ceramic sample data, but with some minor adjustments, can
# be applied to other sample types based on the standards used (i.e. replace the search for 
# "Red" to the name of the standard used in line 85 below)

# determine path to file 
path <- "~file_name.xlsx"
path <- "~file_name.xlsx"

# use map to iterate read_excel over each worksheet in the workbook
ld <- path %>%
        excel_sheets() %>%
        set_names() %>%
        map(read_excel, path = path)

# Bind the columns in the lists together to form one dataframe
df <- bind_cols(ld)

# change name of first column to element
names(df)[names(df) == 'X__1'] <- 'element'

# now, let's get tidy!
# grab the first column as rownames, which will become the variable names
rnames <- df[,1]

# then grab the column names, which will become a new column "Sample" once transposed
Sample <- colnames(df[-1]) #we can drop the first name because it will become the rownames

# transpose the dataframe
df <- t(df[, -1]) #have to drop the first column or it will convert the numbers to strings

# set the column names
colnames(df) <- unlist(rnames) #rnames is stored as a list, so we have to unlist it 

# convert to tibble dataframe
df <- tbl_df(df)

# let's add the date as a column to our data frame so we know when each sample was run
# first we need to figure out how many samples were run each day
ld_lengths <- lapply(ld, length)

# with that information we can create a simple for loop to replicate the dates the appropriate
# number of times for the number of samples run each day
res1 <- as.data.frame(NULL)
for(i in names(ld_lengths)) {
  res <- rep(i, ld_lengths[[i]])
  res1 <- c(res1, res)
}

# now we'll create a data frame of those dates and add it to our sample data
date_col <- tbl_df(sapply(res1, paste0, collapse = ""))
colnames(date_col) <- "Date"
df <- cbind(date_col[2:nrow(date_col),], df)

# add column of samples names, which were the columns names before transposing
df <- cbind(Sample, df)

# use stringr to get rid of repetitive element row names - which have an "X" in them by default since they don't have a column name
dfnames <- df$Sample
x_detect <- str_detect(dfnames, "X")
df <- df[!x_detect, ]

# one pesky column name has a note in it, let's get rid of it too
note <- str_detect(df$Sample, "High")
df <- df[!note, ]

# finally we can convert the sample data to numeric to allow for calculations
df[,3:ncol(df)] <- sapply(df[,3:ncol(df)], as.numeric)

### At this point, you can jump to the "Percent Oxide to PPM.R" routine or just use the ###
#   data as i

# write the full dataframe to a csv
write_csv(df, "full_data_frame.csv")

#now lets get rid of the Ohio Red Standard Samples
dfsamps <- df$Sample
orows <- str_detect(dfsamps, "Red")
df_samples <- df[!orows, ]

#write csv with samples only, Ohio Red standards removed
write_csv(df_samples, "data_samples_only.csv")
