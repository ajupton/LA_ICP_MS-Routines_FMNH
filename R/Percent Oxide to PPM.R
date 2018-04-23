# This routine will allow for easy conversion from % oxide to ppm for LA-ICP-MS data 

# The input file should be a tidy data frame with rows of ICP-MS sample/standard observations
# and columns as the elements

# tidyverse and readxl
library(tidyverse)
library(readxl)
library(stringr)

# First, make sure all elemental/% oxide data are numeric

# convert the sample data to numeric to allow for calculations
df[,3:ncol(df)] <- sapply(df[,3:ncol(df)], as.numeric)

# lets get to work on converting from %oxide to ppm
# each element has a unique coefficient to use when converting, so we'll make a function 
# for each and apply them across the rows
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

# Now we can apply these functions across the appropriate columns
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

# write the full dataframe to a csv
write_csv(df, "full_data_frame_converted_to_ppm.csv")

#now lets get rid of the Ohio Red Standard Samples
dfsamps <- df$Sample
orows <- str_detect(dfsamps, "Red")
df_samples <- df[!orows, ]

#write csv with samples only, Ohio Red standards removed
write_csv(df_samples, "data_samples_only_converted_to_ppm.csv")