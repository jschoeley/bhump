here::i_am('src/_install_dependencies.R'); setwd(here::here())

# specify packages used by project and the repositories from which
# they should be installed
# to find the project dependencies run
# dput(unique(renv::dependencies()$Package))
packages <- list(
  posit = list(
    url = "https://packagemanager.posit.co/cran/2026-04-01/",
    packages = c(
	"here", "readr", "yaml", "svglite", "tidyverse", "dplyr", "rlang",
	"cowplot", "DEoptim", "MASS", "maxLik", "numDeriv", "pracma",
	"qs2", "data.table", "lubridate", "magrittr", "forcats", "scales",
	"DemoDecomp", "gt"
    )
  )
)

# install/update packages
lapply(packages, function(X) install.packages(X$packages, repos = X$url))

# write out list of dependencies and currently used versions
versions <- installed.packages()[
  unlist(lapply(packages, function(X) X$packages)),
  c("Package", "Version")
]
write.csv(versions, "out/_package_versions.csv", row.names = FALSE)
