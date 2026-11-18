# Run from any directory: Rscript --vanilla /path/to/fitrec-heat-running/run.R
local({
  if (.Platform$OS.type=="windows" && !l10n_info()[["UTF-8"]]) Sys.setlocale("LC_CTYPE",".UTF-8")
  args <- commandArgs(trailingOnly=TRUE)
  if (any(!args %in% "--install")) stop("Usage: Rscript --vanilla run.R [--install]")
  script_arg <- grep("^--file=",commandArgs(trailingOnly=FALSE),value=TRUE)
  repo <- if (length(script_arg)) dirname(normalizePath(sub("^--file=","",script_arg[1]),winslash="/")) else getwd()
  library_dir <- file.path(repo,".R-library")
  if (dir.exists(library_dir)) .libPaths(c(library_dir,.libPaths()))
  if (!requireNamespace("data.table",quietly=TRUE) && "--install" %in% args) {
    dir.create(library_dir,showWarnings=FALSE)
    .libPaths(c(library_dir,.libPaths()))
    install.packages("data.table",lib=library_dir,repos="https://cloud.r-project.org")
  }
  if (!requireNamespace("data.table",quietly=TRUE)) {
    stop("The analysis needs data.table. Run Rscript --vanilla run.R --install to install it and run the analysis.")
  }
  library(data.table)
  public_out <- Sys.getenv("FITREC_OUTPUT",unset=file.path(repo,"outputs","reproduced"))
  if (!nzchar(public_out)) public_out <- file.path(repo,"outputs","reproduced")
  dir.create(public_out,recursive=TRUE,showWarnings=FALSE)
  scripts <- c("temperature_models.R","paired_heart_rate.R","heart_rate_decomposition.R","activity_models.R",
               "tables.R","figures.R","check_core_results.R")
  for (script in scripts) {
    message("Running ",script)
    source(file.path(repo,"scripts",script),local=TRUE,encoding="UTF-8")
  }
  writeLines(capture.output(sessionInfo()),file.path(public_out,"sessionInfo.txt"))
  message("Completed. Tables, figures, fitted contrasts and checks are in ",
          normalizePath(public_out,winslash="/"))
})
