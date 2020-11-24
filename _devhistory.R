####################################################
##### Developement history for jk.dusz.tarping #####
####################################################

# ________________________
##### Project set up #####
### First, I created a package named "jk.dusz.tarping" (.Rproj) with the following command:
# usethis::create_package("/Users/francois/Desktop/Recherche & Environnement/Projets Mac/2020_Dusz_Tarping/Data analyses/jk.dusz.tarping")
# in which I specified the ABSOLUTE path to my package's folder (alternatively, I could have clicked the
# buttons in RStudio to create a package with devtools).
# IMPORTANT NOTE: to recreate this project (or reuse the code in another project), you do not need to re-run
# this line of code and recreate the package's folder.

### Then, I created a "dev history" file to keep track of everything I do in the project:
# usethis::edit_file("_devhistory.R")
# Here, I do not need to run it because the present script file (the one you are reading) ALREADY IS my
# development history file!
# From now, all the command lines will be written in this file.
# ________________________


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

