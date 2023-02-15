


# check data/folder package structure
## keep project data external and explain well how to load it (file.path())
### TODO use system.file() and file.path()

# devtools::document() devtools::check() BiocCheck::BiocCheck() as well
## Bioinformatics focused, CRAN more statistics

# get project data folder
## extract sample_names and validate with sample names given by user





# TODO fill in readme.rmd
## fill liberally with info from description

### TIPS
# setwd("package_home/")
## use load_all() (cntrl+shift+L) to simulate building->installing->attaching package

# use_R("file_name") to create new R/file
# (R CMD) check() (cntrl+shift+E)
# cntrol+. to open Rstudio file/function searcher
# Insert Roxygen2 documentation template: cntrl+alt_shift+R
## then document() (cntrl+shift+D)
# install() (& restart) (cntrl+shift+B)
# use_test("file_name") to create new or edit test-file_name
## to execute tests: test() (cntrl+shift+T)
# rename_files("old_file_name", "new_file_name")
# render readme with build_readme() (example: https://github.com/jennybc/regexcite#readme)

## QUESTIONS
# clean and rebuild option? (build() exists?)
# are tests performed automatically on check() / install() / (build()?)

# writing R extensions: https://cran.r-project.org/doc/manuals/R-exts.html#The-DESCRIPTION-file
# R packages: https://r-pkgs.org/Metadata.html
# Git R: http://happygitwithr.com/
# compiled code: https://adv-r.hadley.nz/rcpp.html & https://cpp11.r-lib.org/
# learn Markdown and RMarkdown
# HW advanced R: https://adv-r.hadley.nz/functions.html
# R for Data Science: https://r4ds.had.co.nz/functions.html
