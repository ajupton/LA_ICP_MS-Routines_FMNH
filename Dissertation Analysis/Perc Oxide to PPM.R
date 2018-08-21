## Reading Geochemical data into R and converting from % Oxide to ppm 

# Load packages
library(tidyverse)
library(readxl)
library(stringr)
library(magrittr)

# Determine path to file
path <- "C:/Users/uptonaw/Dropbox/LA-ICP-MS Data Backup/___Calc_ALL RESULTS/all.xlsx"
path <- "all.xlsx"

# Use map to iterate read_excel over each worksheet in the workbook
ld <- path %>%
  excel_sheets() %>%
  set_names(., .) %>% # this was giving me problems, but the two dots is a workaround
  map(read_excel, path = path) 

# Bind the columns in the lists together to form one dataframe
df <- bind_cols(ld)

# Change name of first column to element
names(df)[names(df) == 'X__1'] <- 'element'

# Now, let's get tidy!
# Grab the first column as rownames, which will become the variable names
rnames <- df[,1]

# Then grab the column names, which will become a new column "Sample" once transposed
Sample <- colnames(df[-1]) #we can drop the first name because it will become the rownames

# Transpose the dataframe
df <- t(df[, -1]) #have to drop the first column or it will convert the numbers to strings

# Set the column names
colnames(df) <- unlist(rnames) #rnames is stored as a list, so we have to unlist it 

# Convert to tibble dataframe
df <- tbl_df(df)

# Add the date as a column to our data frame so we know when each sample was run
# first we need to figure out how many samples were run each day
ld_lengths <- lapply(ld, length)

# With that information we can create a simple for loop to replicate the dates the 
# appropriatenumber of times for the number of samples run each day
res1 <- as.data.frame(NULL)
for(i in names(ld_lengths)) {
  res <- rep(i, ld_lengths[[i]])
  res1 <- c(res1, res)
}

# Create a dataframe of those dates and add it to our sample data
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
note <- str_detect(df$Sample, "High")
df <- df[!note, ]

# Convert the sample data to numeric to allow for calculations
df[,3:ncol(df)] <- sapply(df[,3:ncol(df)], as.numeric)

#lets get to work on converting from %oxide to ppm
#each element has a unique coefficient to use when converting, so we'll make a function for each and apply them across the rows
sio2 <- function(x){
  x * 1000000/2.1393 #have to multiply by a million then divide by the coefficient
}

nao2 <- function(x){
  x * 1000000/1.348
}

mgo <- function(x){
  x * 1000000/1.6583
}

al2o3 <- function(x){
  x * 1000000/1.8895
}

p2o5 <- function(x){
  x * 1000000/2.2914
}

k2o <- function(x){
  x * 1000000/1.2046
}

cao <- function(x){
  x * 1000000/1.3992
}

mno <- function(x){
  x * 1000000/1.2912
}

fe2o3 <- function(x){
  x * 1000000/1.4298
}

ti <- function(x){
  x * 1000000/1.6681
}

bao <- function(x){
  x * 1000000/1.1165
}

# Now we can use apply these functions across the appropriate columns
df$SiO2 <- sio2(df$SiO2)
df$Na2O <- nao2(df$Na2O)
df$MgO <- mgo(df$MgO)
df$Al2O3 <- al2o3(df$Al2O3)
df$P2O3 <- p2o5(df$P2O3)
df$K2O <- k2o(df$K2O)
df$CaO <- cao(df$CaO)
df$MnO <- mno(df$MnO)
df$Fe2O3 <- fe2o3(df$Fe2O3)
df$Ti <- ti(df$Ti)
df$BaO <- bao(df$BaO)

#write the full dataframe to a csv
write_csv(df, "Upton_results_samples_and_OhioRed_August_19_2018.csv")

#now lets get rid of the Ohio Red Samples
dfsamps <- df$Sample
orows <- str_detect(dfsamps, "Red")
df_samples <- df[!orows, ]

#write csv with samples only, Ohio Red standards removed
write_csv(df_samples, "Upton_results_samples_calc_August_19_2018.csv")
