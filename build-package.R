# Build package

#use_travis()
#use_cran_comments()
#use_cran_badge()

# Automatically build functions
fun.name <- c(".check_spsann_arguments", ".plotting_options", 
              ".prepare_jittering", ".prepare_points", ".plot_and_jitter",
              ".prepare_output")
read.file <- c("R-autoFunction/check-spsann-arguments.R", 
               "R-autoFunction/plotting-options.R",
               "R-autoFunction/prepare-jittering.R",
               "R-autoFunction/prepare-points.R",
               "R-autoFunction/plot-and-jitter.R",
               "R-autoFunction/prepare-output.R")
write.file <- c("R/check-spsann-arguments.R", "R/plotting-options.R",
                "R/prepare-jittering.R", "R/prepare-points.R",
                "R/plot-and-jitter.R", "R/prepare-output.R")
lapply(1:length(fun.name), function (i) 
  ASRtools::autoFunction(fun.name = fun.name[i], read.file = read.file[i], 
                         write.file = write.file[i]))

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




