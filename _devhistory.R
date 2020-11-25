###################################################
###################################################
##### DEVELOPMENT HISTORY FOR jk.dusz.tarping #####
###################################################
###################################################

# ________________________
##### Project set up #####

### First, I created a package named "jk.dusz.tarping" (.Rproj) with the following command:
# usethis::create_package("/Users/francois/Desktop/Recherche & Environnement/Projets Mac/2020_Dusz_Tarping/Data analyses/jk.dusz.tarping")
# in which I specified the ABSOLUTE path to my package's folder (alternatively, I could have clicked the
# buttons in RStudio to create a package with devtools).
# IMPORTANT NOTE: to recreate this project (or reuse the code in another project), you do not need to re-run
# this line of code and recreate the package's folder.

# Then, I created a "dev history" file to keep track of everything I do in the project:
# usethis::edit_file("_devhistory.R")
# Here, I do not need to run it because the present script file (the one you are reading) ALREADY IS my
# development history file!
# From now, all the command lines will be written in this file.



### First thing first, we will tell R to ignore our _devhistory.R file as it's only for us!
# As this file is not part of a typical package structure, we need to tell R to ignore it when checking
# and installing the package:
usethis::use_build_ignore("_devhistory.R")

# Then, to produce our first commit:
usethis::use_git(message = ":tada: Initial commit") # For a reason I don't understand, there is an error
# saying that my user.name was not configured (then another one for my "user.email"). To solve this, you
# need to run the following commands in the Terminal window (and adapt it with your own user name and
# address):
# git config --global user.name "fanf"
# git config --global user.email "francois.martin74190@gmail.com"
# You can now, rerun the previous line, and restart RStudio (otherwise, the Git tab won't appear)!!!

# If, like me today, you are using a MacOS, you should also do this:
usethis::use_git_ignore(".DS_Store")
usethis::use_build_ignore(".DS_Store")
# It will ignore the .DS_Store hidden file which is specific to MacOS and should not be part of a package.
# And we commit this change:
usethis::use_git(message = ":see_no_evil: Ban .DS_Store files")

# Before we go any further, we will edit some information about our package using the DESCRIPTION file
usethis::edit_file("DESCRIPTION")
usethis::use_git(message = ":bulb: Edit package metadata")

# To create the package's documentation:
usethis::use_package_doc() # It creates a dummy file in the R folder that should NOT be modified!
devtools::document() # Creates the documentation and the man folder (for "manual").
usethis::use_git(message = ":bulb: Update documentation")

# To associate a licence to the package
usethis::use_mit_license(name = "Francois-Marie Martin") # Open-source license
usethis::use_git(message = ":page_facing_up: Add package license")





# ___________________________________________________
##### Writing custom functions for this package #####

### To create a R file (script) for my first custom (personalized) function:
usethis::use_r("import_raw_data")
# NOTE: it automatically places the R file in the R folder (as it should be). The R folder should ONLY
# contain R scripts for functions and NOTHING ELSE!!! Yet, RStudio may sometimes put other things in it,
# so it is a good idea to go and see once in a while.

# To use pipes everywhere in the package without loading the "magrittr" package:
usethis::use_pipe() # Creates automatically a pipe function (and associated R script in the R folder)

# Before writing my function to import my data, I need to add my data files in the project folder.
# So I create a "data" folder and copy-paste my data in it (manually):
dir.create("data")
usethis::use_git(message = ":zap: Created a data and imported data in it")

# Now I can write my function in the associated R file while keeping in mind that:
# *** I need to create a roxygen skeleton to write the documentation;
# *** I need to add the packages I used in my functions within the "Imports" field of the DESCRIPTION;
# *** I should NOT SOURCE my functions (i.e. I should save my function's file but do not run the
# code in it and create "manually" the function you just wrote)! Because we will do it like that
# to avoid conflicts:
devtools::load_all() # Now, all functions in the R folder are available!
usethis::use_git(message = ":metal: New functions: import_raw_data and pipes")

devtools::document() # To create the functions' documentation in the "man" folder, and to update the
# NAMESPACE file of the package (that should NEVER be edited manually).
usethis::use_git(message = ":bulb: Update documentation")



### It's time to check the integrity of our package:
devtools::check() # It gives me 2 warnings (one because I did not declared the "here" and"readr" packages
# in the dependencies; and one because my dataset is not documented). For now, I disregard them.
# IMPORTANT NOTE: I had a lot of PROBLEMS in my first attempts to create this function because the dataset
# contained comments in French and English, with special characters and punctuation (;,:[] etc.) and R
# thought that my punctuation was field separators! Do not put comment in .csv or .txt files!

# To finally install the package
devtools::install()



### To test our functions, we will use the "testthat" package:
usethis::use_testthat()
usethis::use_git(message = ":white_check_mark: Setup testthat")

# Now, we could create some "unit tests" to test our import_raw_data.R function, but we won't.
# usethis::use_test("import_raw_data")
# Here, we don't need to do that because we know our function works as we want it to. For other functions
# and purposes, we should look more closely into that (cf. lesson from N.Casajus FRB-Cesab on package
# building)!
# NOTE: All tests files are stored in tests/testthat/ and their names must start with test-*

usethis::use_git(message = ":bell: Update _devhistory")
