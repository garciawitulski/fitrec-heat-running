# The analysis scripts share one working directory. Once all checks finish,
# copy their outputs into the folders used to navigate the repository.
files <- list.files(out)
is_figure <- tools::file_ext(files) %in% c("pdf", "png", "svg")
is_check <- grepl("^check_|_verification\\.csv$|_audit\\.csv$|_comparison\\.csv$", files)
folders <- rep("derived", length(files))
folders[grepl("curves?|coefficients", files)] <- "derived/coefficients"
folders[files %in% c("heart_bootstrap_numpy.csv", "heart_lee_numpy.csv")] <- "derived/original_draws"
folders[files %in% c("table1.csv", "table2.csv")] <- "tables/manuscript"
folders[files %in% c("cohort_flow.csv", "supp_session.csv", "supp_registered.csv")] <- "tables/supplementary"
folders[is_figure] <- "figures/manuscript"
folders[is_figure & startsWith(files, "figure_S")] <- "figures/supplementary"
folders[is_check | files == "verification_summary.csv"] <- "validation"
folders[files == "sessionInfo.txt"] <- "logs"
for (folder in unique(folders)) {
  dest <- file.path(repo, "outputs", "record_level", folder)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  selected <- files[folders == folder]
  if (!all(file.copy(file.path(out, selected), dest, overwrite = TRUE))) {
    stop("Could not save all files to ", dest)
  }
}
message("Saved the tables, figures and verification files in outputs/record_level/.")
