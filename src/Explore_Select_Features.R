# Plot Select Features of interest

source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/Metadata_prep_funcs.R")
source("src/miscellaneous_funcs.R")
source("src/Community_Composition_Funcs.R")


#-----------------------------------------------------------------------------------------------------
#           This script provides a scaffold to explore selected features of the dataset 



#-------------------------------------------------------------------------------
# Step 1) Select an object type 

# See github table
obj <- dat.genus

#-------------------------------------------------------------------------------
# Step 2) Browse data table 
# (command click on data.table)

data_table <- explore_table(obj)

#-------------------------------------------------------------------------------
# Step 3) Select a feature of interest -
# Make sure spelling is correct and feature is found
# in data table (Step 2)

feature <- "Adlercreutzia"

#-------------------------------------------------------------------------------
# Step 4) Plot a feature of interest 
# Note: uses ArcSin Sqrt Transformation

plot_feature(obj, feature)




