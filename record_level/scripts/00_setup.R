# Paths in the inventory are relative to the directory of licensed inputs.
# A two-column CSV can point individual input names to other local files.
library(data.table)

input_dir <- Sys.getenv("FITREC_INPUT_DIR", unset = file.path(repo, "data", "record_level"))
map_file <- Sys.getenv("FITREC_INPUT_MAP")
input_locations <- character()
if (nzchar(map_file)) {
  locations <- fread(map_file, colClasses = "character")
  stopifnot(identical(names(locations), c("input", "path")),
            !anyNA(locations), all(nzchar(locations$input)), all(nzchar(locations$path)),
            !anyDuplicated(locations$input))
  paths <- path.expand(locations$path)
  relative <- !grepl("^([A-Za-z]:[/\\\\]|[/\\\\])", paths)
  paths[relative] <- file.path(dirname(normalizePath(map_file)), paths[relative])
  input_locations <- setNames(paths, locations$input)
}

input_file <- function(...) {
  keys <- file.path(...)
  paths <- file.path(input_dir, keys)
  mapped <- keys %in% names(input_locations)
  paths[mapped] <- unname(input_locations[keys[mapped]])
  missing <- !file.exists(paths)
  if (any(missing)) {
    stop("Missing record-level input: ", paste(keys[missing], collapse = ", "),
         ". See record_level/INPUTS.md for the file layout and mapping option.")
  }
  paths
}

out <- Sys.getenv("FITREC_OUTPUT", unset = "")
if (!nzchar(out)) out <- tempfile("fitrec-")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
db <- DBI::dbConnect(duckdb::duckdb(home = tempdir()), dbdir = ":memory:")
setDTthreads(2L)

read_fitrec <- function(input, columns = NULL) {
  path <- input_file(input)
  fields <- if (is.null(columns)) "*" else {
    paste(DBI::dbQuoteIdentifier(db, columns), collapse = ", ")
  }
  sql <- paste("SELECT", fields, "FROM read_parquet(",
               DBI::dbQuoteString(db, normalizePath(path, winslash = "/")), ")")
  as.data.table(DBI::dbGetQuery(db, sql))
}
