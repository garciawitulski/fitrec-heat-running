# Run from the repository root with Rscript --vanilla record_level/run.R.
local({
  if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) {
    Sys.setlocale("LC_CTYPE", ".UTF-8")
  }
  repo <- normalizePath(Sys.getenv("FITREC_REPO", unset = getwd()), winslash = "/")
  script_dir <- file.path(repo, "record_level", "scripts")
  if (!file.exists(file.path(script_dir, "00_setup.R"))) {
    stop("Run this script from the repository root or set FITREC_REPO.")
  }
  source(file.path(script_dir, "00_setup.R"), local = TRUE, encoding = "UTF-8")
  tryCatch({
    source(file.path(repo, "scripts", "model_tools.R"), local = TRUE, encoding = "UTF-8")
    source(file.path(script_dir, "panel_checks.R"), local = TRUE, encoding = "UTF-8")
    scripts <- sort(list.files(script_dir, pattern = "^[0-9].*\\.R$"))
    for (script in setdiff(scripts, "00_setup.R")) {
      message("\n", script)
      source(file.path(script_dir, script), local = TRUE, encoding = "UTF-8")
    }
    writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
    if (!nzchar(Sys.getenv("FITREC_OUTPUT"))) {
      source(file.path(script_dir, "save_outputs.R"), local = TRUE, encoding = "UTF-8")
    }
  }, finally = DBI::dbDisconnect(db, shutdown = TRUE))
})
