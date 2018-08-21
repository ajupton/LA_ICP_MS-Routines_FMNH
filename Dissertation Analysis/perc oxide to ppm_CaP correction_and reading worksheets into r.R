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

##----------------------------Correction for Shell Tempering Here-------------------------##

#In analysis, I need to correct the sherd samples for the presence of shell tempering.
#Shell is composed almost entirely of calcium which is in the same row in the periodic table
#as strontium and barium. 

#First step is to drop the Ohio Red standard samples because they don't need correcting
orows <- str_detect(tolower(df$Sample), "red")
df_samples1 <- df[!orows, ]

#Add up all elements calculated in percent oxide aside from Ca and Ba
CaP_correction <- df_samples1 %>% 
                    select(SiO2, Na2O, MgO, Al2O3, K2O, Sb2O5, 
                           MnO, Fe2O3, CuO, SnO2, Ti, PbO2, BaO, Bi, ZnO) %>% 
                    rowSums() %>%
                    tbl_df()

#Correct the elements by dividing their amount by the corrected percent oxide
df_samples1[, c(3:length(df_samples1))] <- sapply(df_samples1[, c(3:length(df_samples1))], 
                                 function(x){x/CaP_correction}) %>%
                            bind_cols()

# Bind the shell corrected ceramic samples with the Ohio Reds
df_shell_corrected <- tbl_df(bind_rows(df_samples1, df[orows, ]))

####------------------------------------------------------------------------------------####

# Converting from %oxide to ppm
# Each element has a unique coefficient to use when converting, so we'll make a 
# function for each and apply them across the rows
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

# Apply these functions across the appropriate columns
df_shell_corrected$SiO2 <- sio2(df_shell_corrected$SiO2)
df_shell_corrected$Na2O <- nao2(df_shell_corrected$Na2O)
df_shell_corrected$MgO <- mgo(df_shell_corrected$MgO)
df_shell_corrected$Al2O3 <- al2o3(df_shell_corrected$Al2O3)
df_shell_corrected$P2O3 <- p2o5(df_shell_corrected$P2O3)
df_shell_corrected$K2O <- k2o(df_shell_corrected$K2O)
df_shell_corrected$CaO <- cao(df_shell_corrected$CaO)
df_shell_corrected$MnO <- mno(df_shell_corrected$MnO)
df_shell_corrected$Fe2O3 <- fe2o3(df_shell_corrected$Fe2O3)
df_shell_corrected$Ti <- ti(df_shell_corrected$Ti)
df_shell_corrected$BaO <- bao(df_shell_corrected$BaO)

# Since we've converted from %oxide, it's a good idea to change the element names
# Some "O's" are left to differentiate the elements measured as both %oxide and not
names(df_shell_corrected) <- c("Sample", "Date","Si","Na","Mg","Al","P","Cl","K","Ca","SbO","Mn",
                               "Fe","CuO","Sn","Ti","Pb","Ba","Bi","ZnO","Li","Be", "B","P","Cl1",
                               "Sc","Ti1","V","Cr","Mn","Fe","Ni", "Co","Cu","Zn","As","Rb","Sr",
                               "Zr","Nb","Ag","In","Sn","Sb","Cs","Ba","La","Ce","Pr","Ta","Au",
                               "Y","Pb","Bi1","U","W","Mo","Nd","Sm","Eu","Gd","Tb","Dy","Ho",
                               "Er","Tm","Yb","Lu","Hf","Th")  

# Write the full dataframe to a csv
write_csv(df_shell_corrected, 
          "Upton_results_samples_and_OhioRed_shell_corrected_all_elements_August_21_2018.csv")

# Drop the Ohio Red Samples
dfsamps <- tolower(df_shell_corrected$Sample)
orows <- str_detect(dfsamps, "red")
df_samples <- df_shell_corrected[!orows, ]
df_reds <- df_shell_corrected[orows, ]

# Write csv with samples only, Ohio Red standards removed
write_csv(df_samples, "Upton_results_samples_shell_corrected_August_21_2018.csv")

# Write csv with Ohio Reds only
write_csv(df_reds, "Upton_results_OhioRed_August_21_2018.csv")
