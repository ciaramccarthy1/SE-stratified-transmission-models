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
### Currently RSV-focused; uncomment Influenza/COVID-19 blocks (and the
### all-diseases plot section in modelrun.r) to bring them back.
# pset$Disease <- "Influenza"
# source(file = paste0(source_dir,"/modelrun.r")) #repo
#
# pset$Disease <- "COVID-19"
# source(file = paste0(source_dir,"/modelrun.r")) #repo

pset$Disease <- "RSV-illness"
#   run model, plot figures, write summaries
source(file = paste0(source_dir,"/modelrun.r")) #repo

################################################################################
