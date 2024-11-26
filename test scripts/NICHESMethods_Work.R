# Create NICHESMethods 2024-11-25

# Initial creation of package
library(devtools)
# create_package("/Users/msbr/GitHub/NICHESMethods") # only need to run once at the outset

# Step 1: Pull most recent NICHESMethods branch from cloud with GitHub desktop

# Step 2: Set WD to GitHub/NICHESMethods
setwd("/Users/msbr/GitHub/NICHESMethods")

# Step 3: Load the package from local copy
getwd() # confirm that in the right place
load_all() # load local package
check() # check that it works
