# Compare each estimate, confidence limit and p-value with its reference value.
# Read the reference after fitting and write only aggregate comparisons.
compare_aggregate <- function(current, reference, keys, fields, label, tolerance = 1e-6) {
  current <- copy(current)
  reference <- copy(reference)
  for (key in keys) {
    # Temperature-band contrasts and joint tests can have an empty label.
    if (is.character(current[[key]]) || is.character(reference[[key]])) {
      set(current, j = key, value = fcoalesce(as.character(current[[key]]), ""))
      set(reference, j = key, value = fcoalesce(as.character(reference[[key]]), ""))
    }
  }
  stopifnot(!anyDuplicated(current, by = keys), !anyDuplicated(reference, by = keys))
  both <- merge(current, reference[, c(keys, fields), with = FALSE], by = keys,
                all = TRUE, suffixes = c("_R", "_paper"))
  answers <- lapply(fields, function(field) {
    a <- both[[paste0(field, "_R")]]
    b <- both[[paste0(field, "_paper")]]
    ans <- copy(both[, ..keys])
    ans[, `:=`(field = field, R = a, paper = b, absolute_difference = abs(a - b))]
    ans[, status := fcase(is.na(a) & is.na(b), "NOT_APPLICABLE",
                          is.na(a) | is.na(b), "MISSING",
                          abs(a - b) <= tolerance, "MATCH", default = "MISMATCH")]
    ans
  })
  checked <- rbindlist(answers)
  fwrite(checked, file.path(out, paste0("check_", label, ".csv")))
  print(checked[, .N, by = status])
  invisible(checked)
}

# Calculate a contrast on the model scale. The named vector terms holds its weights.
contrast_panel <- function(fit, terms) {
  stopifnot(all(names(terms) %in% names(fit$beta)))
  w <- setNames(numeric(length(fit$beta)), names(fit$beta))
  w[names(terms)] <- terms
  z <- lin_contrast(fit, w)
  data.table(estimate = z$estimate, se = sqrt(max(0, drop(crossprod(w, fit$vcov %*% w)))),
             ci_low = z$lower, ci_high = z$upper, p_value = z$p)
}
