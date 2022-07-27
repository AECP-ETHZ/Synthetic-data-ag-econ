# ---------------------------------------------------------------------------- #
#                                                                              #
# R code for the article: "A note on synthetic data for replication purposes   #
# in agricultural economics"                                                   #
#                                                                              #
# Authors: Stefan Wimmer and Robert Finger                                     #
#                                                                              #
# Corresponding author: Stefan Wimmer (swimmer@ethz.ch)                        #
#                                                                              #
# Citation: Wimmer, S., Finger, R. (2023). "A note on synthetic data for       #
# replication purposes in agricultural economics." Journal of Agricultural     #
# Economics.                                                                   #
#                                                                              #
# This program prepares the FADN data for production function and frontier     #
# estimation. The data (stored in the subfolder 'Data') are not publicly       #
# available but can be requested from DG-Agri. The script was successfully     #
# run with R version 4.2.0.                                                    #
#                                                                              #
# The codes are also published in the following github repository:             #
# https://github.com/AECP-ETHZ/Synthetic-data-ag-econ                          #
#                                                                              #
# ---------------------------------------------------------------------------- #

# Content:
# 0) Preparation (working directory, packages, data)
# 1) Generate output variable
# 2) Generate productive input variables
# 3) Make homogeneous sample
# 4) Save the sample

# -------------------- #
#### 0. Preparation ####
# -------------------- #

# Set working directory
my.wd <- getwd() # NOTE: getwd() can be replaced with any working directory
setwd(my.wd)

# Create subfolders for R output
dir.create(file.path("rOutput"))

# Install missing packages
list.of.packages <- c("dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load packages
library(dplyr) # for data management

# Load FADN data and keep crop farms only
dat <- read.csv("Data/DEU2018.csv") %>% 
  filter(TF8 == 1)

# -------------------------------------------- #
#### 1. Generate productive input variables ####
# -------------------------------------------- #
  
  #------------------#
  ##### 1.1 Land #####
  #------------------#

  # Define land as total utilised agricultural area  
  dat$land <- dat$SE025
  
  #-------------------#
  ##### 1.2 Labor #####
  #-------------------#
  
  # Define labor as the sum of unpaid and paid labor input
  dat$labor <- dat$SE015 + dat$SE020
  
  #----------------------#
  ##### 1.3 Material #####
  #----------------------#

  # Define material as variable costs
  dat$material <- dat$SE285 + dat$SE295 + dat$SE300 + dat$SE305 + # crops related inputs until here
                  dat$SE310 + dat$SE320 + dat$SE330 + # livestock related inputs until here 
                  dat$SE340 + dat$SE345 + dat$SE350 + dat$SE356 # other materials until here
                  
  #---------------------#
  ##### 1.4 Capital #####
  #---------------------#
  
  # Define capital stock as average farm capital (which excludes land)
  dat$capstock <- dat$SE510

# --------------------------------- #
#### 2. Generate output variable ####
# --------------------------------- #

# Define output as the sum of all crop outputs
dat$output <- dat$SE140 + dat$SE145 + dat$SE150 + dat$SE155 +
              dat$SE160 + dat$SE165 + dat$SE170 + dat$SE175 +
              dat$SE180 + dat$SE185 + dat$SE190 + dat$SE195 + 
              dat$SE200

  #note: this is "total crop output", i.e. revenues for cereals, protein, 
  #      potatoes, sugar beet, oilseed, industrial crops, vegetables, 
  #      fruit trees, citrus fruit, wine grapes, olives, forage crops, 
  #      and other crop output

# Drop farms with non-positive crop output
dat <- dat %>% 
  filter(output > 0) # --> only one observation dropped

# ------------------------ #
#### 3. Save the sample ####
# ------------------------ #

save(dat, file="rOutput/FADN.Rds")
  
# --------------------------------------- #
#            End of the script            #
# --------------------------------------- #
  

    