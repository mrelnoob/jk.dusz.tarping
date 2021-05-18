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
# From now on, all command lines will be written in this file.
# IMPORTANT NOTE: in this case, I started a package project from scratch. If I want to start from a
# GitHub repository, the workflow is different but some steps will be required so I encourage my future
# self to carefully review the steps of this workflow anyway.



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





# ________________________________________________________________________
##### How to create my first custom functions and install my package #####

### To create a R file (script) for my first custom (personalized) function:
usethis::use_r("00_import_raw_data")
# NOTE: it automatically places the R file in the R folder (as it should be). The R folder should ONLY
# contain R scripts for functions and NOTHING ELSE!!! Yet, RStudio may sometimes put other things in it,
# so it is a good idea to go and see once in a while.

# To use pipes (i.e. %>%) everywhere in the package without explicitly loading the "magrittr" package:
usethis::use_pipe() # Automatically creates a pipe function (and associated R file in the R folder),
# while updating the DESCRIPTION file to tell R that it should import the "magrittr" package.

# Before writing my function to import my data, I need to add my data files in the project folder.
# So I create a "mydata" folder and copy-paste my data in it (manually):
dir.create("mydata")
usethis::use_build_ignore("mydata") # Because it is not expected in a regular package root folder (if
# I don't ignore it, it will cause warnings in my package checks and all kind of craps. I just spent
# five hours to get rid of them, so trust me).
usethis::use_git(message = ":zap: Created a mydata folder and imported data in it")



### Now I can write my function in the associated R file while keeping in mind that:
#  --> I need to create a roxygen skeleton to write the documentation.
#  --> I should NOT SOURCE my functions (i.e. I should save my function's file but I should not run the
#    code in it and create "manually" the function you just wrote) to avoid conflicts, see below!
# IMPORTANT NOTE: if a function uses functions from other packages, you need to tell R about it by updating
# the NAMESPACE file. It will be done automatically by devtools when we produce our documentation if we
# have previously listed the dependencies (packages) in the Roxygen2 header of the functions thanks to the
# tags (#' @import package OR #' @importFrom package function)! So we do this and add in our function's
# Roxygen2 header the required tags and packages, and we DO NOT FORGET to also add these
# dependencies in the "Imports" field of the DESCRIPTION:
usethis::use_package("readr")
usethis::use_package("here")
usethis::use_package("stringr")
# REMINDER: The NAMESPACE controls what happens when our package is loaded but not when it's installed.
# This is the role of DESCRIPTION!

############################################################################################################
################# IMPORTANT PART: to load, check and install the functions of my package ###################
############################################################################################################

### The following lines (until the part about connecting to Github) are run iteratively every time a new
# function is created within the package or every time a function is modified. This is the central part
# that loads, inspects (i.e. run the R CMD check), test, and install my functions and my whole package!

# To load our functions, we will thus use:
devtools::load_all() # Now, all functions in the R folder are available!

devtools::document() # To create the functions' documentation in the "man" folder, and to update the
# NAMESPACE file of the package (that should NEVER be edited manually).
usethis::use_git(message = ":bulb: Update documentation")



### To test our functions, we will use the "testthat" package:
#usethis::use_testthat()
#usethis::use_git(message = ":white_check_mark: Setup testthat")

# Now, we could create some "unit tests" to test our 00_import_raw_data.R function, but we won't modify it.
#usethis::use_test("import_raw_data")
# Here, we don't want to do real tests because we know our function works as we want it to. For other
# functions and purposes, we should look more closely into that (cf. lesson from N.Casajus FRB-Cesab on
# package building)!
# NOTE: All tests files are stored in tests/testthat/ and their names must start with test-*



### When everything is ready, it's time to check the integrity of our package:
devtools::check() # Ok!
# IMPORTANT NOTE: I had a lot of PROBLEMS in my first attempts to create the loading function because the
# dataset contained comments in French and English, with special characters and punctuation (;,:[] etc.)
# and R thought that my punctuation was field separators! NEVER put comments in .csv or .txt files!!!

# To finally install the package:
devtools::install()
usethis::use_git(message = ":metal: Updated function: cond.boxplots()!")
system("git push")



### We will now update the version of our package:
usethis::use_version(which = "minor") # It automatically updates our package version.
usethis::use_news_md() # Creates a NEWS.md file, that I should maintain updated.
usethis::use_git(message = ":package: Release v0.1.0")





# _____________________________________________________________________________________
##### How to connect the package to my GitHub account (in a GitHub last workflow) #####

### Get a a GitHub Personal Access Token (PAT) which will allow R to control GitHub (only once):
usethis::browse_github_token() # First, I generate a token and copy its value, and follow the instructions
# displayed in the R console (edit the R environment etc.)...
usethis::edit_r_environ() # ... to manually add the PAT to the R environment, save, close, and restart R...
# Second, after restarting R, I verify:
usethis::github_token()
# OK! So we are ready to connect our local project to a new Github repository:
usethis::use_github(protocol = "ssh") # I use this option as I have already generated SSH keys
# for this computer (?). Yet, it will fail, because we need to run the following command:
system("git push --set-upstream origin master") # Or without system() but in the Terminal.
# Now our local project and Github are linked and synchronized (for now). After committing future
# changes we will be able to push normally (using RStudio git panel or git push).





# ____________________________________________________
##### How to add a README file to our repository #####

usethis::use_readme_rmd() # Creates a README.Rmd and adds it automatically to .Rbuildignore (and opens it).
# After manually editing the file, we need to compile it into a .md document (otherwise, GitHub and
# the CRAN won't be able to read it and display it on their websites):
rmarkdown::render("README.Rmd")
# As render() also produces a .html file that is not useful here, we will ignore it:
usethis::use_build_ignore("README.html")
usethis::use_git_ignore("README.html")

usethis::use_git(message = ":pencil: Edit README")
system("git push")
# IMPORTANT NOTE: Each time you edit the README.Rmd you will have to update the .md
# with rmarkdown::render("README.Rmd") and, of course, you should also commit+push it to update GitHub!





####################################
####################################
##### Using the Drake workflow #####
####################################
####################################

# Up to this point, I've created a working package called `jk.dusz.tarping` that, as of today (14/04/21),
# contains the functions used to produce my data preparation report. It is not very useful nor compulsory
# to have a proper _devhistory.R file nor to use Drake. But I wanted to have a very clean and well
# commented "development history" file to help me build more useful packages in the future. I'll delete it
# if I finally take the time to change this file into a proper tutorial (in RMarkdown).

# Now, I will try to carry on my data analysis process while and place it into a Drake workflow to ensure
# the reproducibility, consistency and robustness of my process. I am not sure if I will integrate the
# following part into my package (because it's fairly time consuming) but I will try.
# It should not be a problem not to transform everything as a package. As long as I don't use the commands
# meant for package building (such as devtools::check() or document()), there should be no reasons for
# errors (I guess)...



# _____________________________________________________
##### STEP 1 - Preparing the Drake project set up #####

### I have to create the scripts for the different elements of my plan:
usethis::use_r(name = "01_data_cleaning.R")
usethis::use_r(name = "02_data_preparation.R") # These first 2 scripts are what others call "wrangle"
usethis::use_r(name = "03_modelling_eff.R")
usethis::use_r(name = "04_modelling_high_eff.R")
usethis::use_r(name = "05_modelling_redges.R")
usethis::use_r(name = "06_modelling_rtarped.R")
usethis::use_r(name = "07_modelling_roverlaps.R")
usethis::use_r(name = "plan.R") # Do not forget to create the "plan" of your Drake pipeline

# Then I have to create the make.R file (kind of master script) and a "_drake.R" file (but I'm not
# sure why yet, possibly to contain the logs):
file.create("make.R")
file.create("_drake.R")

# To create new folders (directories) to store the various kind of results my analysis will produce:
dir.create("output")
dir.create("output/plots")
dir.create("output/tables")
dir.create("output/texts")

# I also need to ignore most of these files:
usethis::use_build_ignore(".drake")
usethis::use_build_ignore("_drake.R")
usethis::use_build_ignore("make.R")
usethis::use_git_ignore(".drake")
usethis::use_build_ignore("output/")
usethis::use_build_ignore("tables/")
usethis::use_build_ignore("texts/")
usethis::use_build_ignore("plots/")





# ________________________________________________________
##### STEP 2 - Populating my data processing scripts #####

### To analyse my data, I will use various packages (dependencies), so I need to fill-in the DESCRIPTION
# file about them in order to load them afterwards with devtools (in case I don't carry on with the
# package creation):
usethis::use_package("drake")
usethis::use_package("dplyr")
usethis::use_package("ggplot2")
usethis::use_package("forcats")
usethis::use_package("fishualize")
usethis::use_package("Hmisc")
usethis::use_package("ggcorrplot")
usethis::use_package("ade4")
usethis::use_package("FactoMineR")
usethis::use_package("factoextra")
usethis::use_package("gridExtra")
#devtools::load_all() # Not sure if I truly need to do this...



### Up to this point, I have created the "cleaning" functions for the first part of my data processing.
# So my 01_data_cleaning.R file is populated with some functions and the package works. So NOW, I need
# to start addressing my data preparation step and start keeping tracks of my data wrangling and future data
# analyses in a .Rmd report:
file.create(... = "output/texts/jk.d.tarp.data_preparation_report.Rmd") # Using this command, the .Rmd file
# will be created but will lack the YAML header skeleton that should thus be manually placed at the top of
# the document.
usethis::use_git(message = ":metal: First results for lreg_edges")
system("git push")



##### FOR DATA WRANGLING: do not forget to EMPTY your environment and DELETE files after wrangling!

