# Read, label and concatenate data on US births, infant- and fetal deaths
#
# (1) Read NCHS data into R applying variable specifications stored
#     in custom codebook
# (2) Concatenate NCHS data across multiple years

# Init --------------------------------------------------------------------

here::i_am('src/20-prepare_us_fetoinfants.R'); setwd(here::here())

library(yaml); library(dplyr)
library(qs2)

#memory.limit(64000)

# Constants ---------------------------------------------------------------

paths <- list()
paths$input <- list(
  codebook_functions.R = 'src/00-codebook.R',
  codebook_infant.yaml = 'cfg/codebook-us_cohort_infant_births_deaths_minimal.yaml',
  codebook_fetus.yaml = 'cfg/codebook-us_fetal_deaths_minimal.yaml',
  zip_infant = 'dat/10-nchs-us_cohort_linked_infant_deaths_births/',
  zip_fetus = 'dat/10-nchs-us_fetal_deaths/'
)
paths$output <- list(
  fetoinfant.qs = 'tmp/20-fetoinfant.qs'
)

# codebook function
source(paths$input$codebook_functions.R)

# Read codebook -----------------------------------------------------------

codebook <- list()

codebook$fetus <- ReadCodebook(paths$input$codebook_fetus.yaml)
codebook$infant <- ReadCodebook(paths$input$codebook_infant.yaml)

# Read data into R and apply varspecs -------------------------------------

infant <- ReadFromZip(
  codebook$infant, paths$input$zip_infant, subset = c(
    'Cohort1989','Cohort1990',
    'Cohort1999','Cohort2000',
    'Cohort2009','Cohort2010',
    'Cohort2014','Cohort2015'
  ))
fetus <- ReadFromZip(
  codebook$fetus, paths$input$zip_fetus, subset = c(
    'Period1989','Period1990',
    'Period1999','Period2000',
    'Period2009','Period2010',
    'Period2014','Period2015'
  ))

# Concatenate data --------------------------------------------------------

# merge data on births, fetal- and infant deaths cross years
fetoinfant <-
  bind_rows(
    infant = bind_rows(infant),
    fetus = bind_rows(fetus),
    .id = 'type'
  )

# Export ------------------------------------------------------------------

qs_save(
  fetoinfant,
  file = paths$output$fetoinfant.qs
)
