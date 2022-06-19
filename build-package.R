# Optimization of Spatial Samples via Simulated Annealing
# Alessandro Samuel-Rosa

# Test package #####################################################################################
# Dependencies
update(remotes::package_deps("spsann"))
update(remotes::package_deps("devtools"))
if (!require(autofun)) {
  remotes::install_github(repo = "samuel-rosa/autofun")
}
# Reverse dependency tools
devtools::revdep("spsann-package")
# Render README
rmarkdown::render("spsann-package/README.Rmd")
# Automatically build functions
fun_name <- c(
  ".check_spsann_arguments",
  ".plotting_options",
  ".prepare_jittering",
  ".prepare_points",
  ".plot_and_jitter",
  ".prepare_output",
  ".prepare_acdc_covars",
  ".check_suggests",
  ".prepare_clhs_covars",
  ".set_progress",
  ".update_progress",
  ".check_first_chain")
read_file <- c(
  "check-spsann-arguments.R",
  "plotting-options.R",
  "prepare-jittering.R",
  "prepare-points.R",
  "plot-and-jitter.R",
  "prepare-output.R",
  "prepare-acdc-covars.R",
  "check-suggests.R",
  "prepare-clhs-covars.R",
  "set-progress.R",
  "update-progress.R",
  "check-first-chain.R")
read_file <- paste0("spsann-package/R-autofun/", read_file)
write_file <- c(
  "check-spsann-arguments.R",
  "plotting-options.R",
  "prepare-jittering.R",
  "prepare-points.R",
  "plot-and-jitter.R",
  "prepare-output.R",
  "prepare-acdc-covars.R",
  "check-suggests.R",
  "prepare-clhs-covars.R",
  "set-progress.R",
  "update-progress.R",
  "check-first-chain.R")
write_file <- paste0("spsann-package/R/", write_file)
lapply(seq_along(fun_name), function(i) {
  autofun::autofun(fun_name[i], read_file[i], write_file[i])
})
rm(fun_name, read_file, write_file)

# Rcpp::compileAttributes()

# check documentation ----
roxygen2::roxygenise("spsann-package")
devtools::check_man("spsann-package")
devtools::spell_check("spsann-package")
# spelling::update_wordlist("spsann-package")

# check examples ----
devtools::run_examples("spsann-package/")

# check for Linux (local) ----
devtools::check("spsann-package/",
  env_vars = c(`_R_CHECK_DEPENDS_ONLY_` = TRUE))
devtools::check("spsann-package/",
  document = TRUE, manual = TRUE, vignettes = TRUE, force_suggests = TRUE, incoming = TRUE,
  remote = TRUE)

# check for Windows (remote) ----
devtools::check_win_oldrelease()
devtools::check_win_release()
devtools::check_win_devel()

# check in R-hub ----
# rhub::validate_email(email = "alessandrosamuelrosa@gmail.com")
# rhub::check_on_windows()
# rhub::platforms()
platforms <- c("fedora-clang-devel",
  "ubuntu-gcc-release", "debian-clang-devel", "windows-x86_64-devel")
devtools::check_rhub(platforms = platforms)

devtools::build()

# Load package
devtools::load_all("spsann-package")

# upload to CRAN
devtools::release(check = FALSE)
