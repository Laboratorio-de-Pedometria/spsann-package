# Build package

# Automatically build functions
require(autofun)
fun.name <- c(".check_spsann_arguments", ".plotting_options", 
              ".prepare_jittering", ".prepare_points", ".plot_and_jitter",
              ".prepare_output", ".prepare_acdc_covars")
read.file <- c("R-autofun/check-spsann-arguments.R", 
               "R-autofun/plotting-options.R",
               "R-autofun/prepare-jittering.R",
               "R-autofun/prepare-points.R",
               "R-autofun/plot-and-jitter.R",
               "R-autofun/prepare-output.R",
               "R-autofun/prepare-acdc-covars.R")
write.file <- c("R/check-spsann-arguments.R", "R/plotting-options.R",
                "R/prepare-jittering.R", "R/prepare-points.R",
                "R/plot-and-jitter.R", "R/prepare-output.R",
                "R/prepare-acdc-covars.R")
lapply(1:length(fun.name), function (i) 
  autofun::autofun(fun.name[i], read.file[i], write.file[i]))

# turn on/off development mode
devtools::dev_mode()

# check examples and documentation
devtools::run_examples()
devtools::check_doc()

# check the package for Linux and Windows
devtools::check()
devtools::build_win()

# build package
devtools::build()

# upload to CRAN
devtools::release()




