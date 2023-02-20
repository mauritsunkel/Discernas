output_dir <- file.path("EMC-SKlab-scRNAseq", "results", "test_pipe")
# files and sample names
sample_files <- c("C:/Users/mauri/Desktop/Single Cell RNA Sequencing/EMC-SKlab-scRNAseq/results/Pipe_SCTv2_23-06 cluster_level_selection/BL_A/BL_A.rds",
                  "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/EMC-SKlab-scRNAseq/results/Pipe_SCTv2_23-06 cluster_level_selection/BL_C/BL_C.rds")
sample_names <- c("BL_A", "BL_C")
selection_panel <- c("VIM", "S100B", "SOX9")
# reticulate::use_python("C:/Users/mauri/AppData/Local/Programs/Python/Python38/python.exe")
EMC_SKlab_scRNAseq::samples_integration(sample_files, sample_names, output_dir,
                    selection_panel = selection_panel)

# could not find function "grid.arrange"
#  Called from: do.call("grid.arrange", c(plot_list, ncol = 2))






integrated <- readRDS("C:/Users/mauri/Desktop/test.rds")
# run FindClusters()

reticulate::use_python("C:/Users/mauri/AppData/Local/Programs/Python/Python38/python.exe")

G <- igraph::graph.famous('Zachary')
leiden::leiden(G)

library(leiden) # TODO test if this changes python path used

Cluster_assignment <- leiden(integrated@graphs$integrated_snn)

# reinstall leiden package from CRAN instead of Github version
## check to reinstall igrpah too (CRAN or not?)

# pip install leidenalg

# last time I got it working I think it was someting about the order of installation/running lol

# try method = 'matrix' instead of 'igrpah' now that you have more RAM?

# see leidenbase package for running leiden in C instead of Python (reticulate) dependency

# https://github.com/satijalab/seurat/issues/3168


# test  py_module_available(module = "leidenalg")
## still seems than once running FindClusters it seems to call it differently somehow

# put reticulate in package description imports?

# install leidenalg via pip to get in local Python library
## NOT FROM CONDA ENVIRONMENT
### uninstall conda!

# try leidenalg==0.8.8

usethis::edit_r_environ()
usethis::edit_r_profile()
# save new line in Renviron with the following
RETICULATE_PYTHON = "/path_here"
Sys.setenv(RETICULATE_PYTHON = "/usr/bin/python3")


# restart R and check to see if reticulate is pointed to the right place
reticulate::py_discover_config() # USE INSTEAD OF reticulate::py_config()

# see https://github.com/vtraag/leidenalg#troubleshooting
