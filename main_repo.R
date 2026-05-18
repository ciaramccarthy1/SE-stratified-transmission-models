################################################################################

library(here)

### folders
input_dir  <- paste0(here("data"))
source_dir <- paste0(here("codes"))
output_dir <- paste0(here("output"))
#fs::dir_create(output_dir)

### Basic setting
source(file = paste0(source_dir,"/setup.r")) #repo

### Diseases cycle
pset$Disease <- "Influenza"
#   run model, plot figures, write summaries
source(file = paste0(source_dir,"/modelrun.r")) #repo

pset$Disease <- "COVID-19"
#   run model, plot figures, write summaries
#pset$COMPILE <- 0
source(file = paste0(source_dir,"/modelrun.r")) #repo

pset$Disease <- "RSV-illness"
#   run model, plot figures, write summaries
#pset$COMPILE <- 0
source(file = paste0(source_dir,"/modelrun.r")) #repo

################################################################################
