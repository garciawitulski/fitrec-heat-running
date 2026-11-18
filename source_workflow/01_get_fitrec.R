# Run from the root of the replication repository. An optional argument points
# to a directory that already contains the two original FitRec downloads.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) <= 1L)
raw_dir <- if (length(args)) args[1L] else file.path("data", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
manifest <- data.frame(
  name = c("endomondoHR.json.gz", "endomondoMeta.json.gz"),
  bytes = c(1995378519, 1088266461),
  sha256 = c("31a49f8ec2fa5d038983e7b3bfe35f5f7eabee581180c803023f011f046681b9",
             "fb6956c506aa068f2f7d2a64478a286b5d663e45052caca1020988cffc475b6f")
)
origin <- "https://mcauleylab.ucsd.edu/public_datasets/gdrive/fitrec/"
options(timeout = max(14400, getOption("timeout")))
for (i in seq_len(nrow(manifest))) {
  dest <- file.path(raw_dir, manifest$name[i])
  new_download <- !file.exists(dest)
  path <- if (new_download) tempfile("fitrec-", tmpdir = raw_dir) else dest
  if (new_download) {
    message("Downloading ", manifest$name[i], " from UC San Diego")
    tryCatch(
      utils::download.file(paste0(origin, manifest$name[i]), path, mode = "wb", method = "libcurl"),
      error = function(e) {
        unlink(path)
        stop(conditionMessage(e), call. = FALSE)
      }
    )
  }
  if (file.info(path)$size != manifest$bytes[i] ||
      tolower(unname(tools::sha256sum(path))) != manifest$sha256[i]) {
    if (new_download) unlink(path)
    stop("The size or SHA-256 hash differs for ", manifest$name[i],
         ". Check the source file against the FitRec release used for the paper.")
  }
  if (new_download) {
    if (file.exists(dest) || !file.rename(path, dest)) {
      unlink(path)
      stop("Could not install the verified file. Check the destination and rerun.")
    }
  }
  message("Verified ", manifest$name[i], " (", manifest$bytes[i], " bytes)")
}
