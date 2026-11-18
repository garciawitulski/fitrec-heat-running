# Shared functions for the models in Supplement B.
# Design matrices have one row per observation and one column per covariate.

rcs_basis <- function(x, knots, center = 20, scale = 10) {
  stopifnot(length(knots) >= 3L, all(diff(knots) > 0), scale > 0)
  x <- (x - center) / scale
  k <- (knots - center) / scale
  first <- k[1L]
  last <- tail(k, 1L)
  penultimate <- k[length(k) - 1L]
  basis <- matrix(x, ncol = 1L)
  for (j in seq_len(length(k) - 2L)) {
    # Subtract the value at zero to centre the cubic term as well.
    cubic <- function(z) {
      (pmax(z - k[j], 0)^3 -
         (last - k[j]) / (last - penultimate) * pmax(z - penultimate, 0)^3 +
         (penultimate - k[j]) / (last - penultimate) * pmax(z - last, 0)^3) /
        (last - first)^2
    }
    basis <- cbind(basis, cubic(x) - cubic(0))
  }
  colnames(basis) <- c("linear", paste0("cubic",seq_len(ncol(basis)-1L)))
  basis
}

season_terms <- function(time_utc, south) {
  date <- as.POSIXct(time_utc, origin = "1970-01-01", tz = "UTC")
  year <- as.integer(format(date, "%Y", tz = "UTC"))
  years <- sort(unique(year[!is.na(year)]))
  starts <- as.numeric(as.POSIXct(paste0(years, "-01-01"), tz = "UTC"))
  ends <- as.numeric(as.POSIXct(paste0(years + 1L, "-01-01"), tz = "UTC"))
  i <- match(year, years)
  phase <- (as.numeric(date) - starts[i]) / (ends[i] - starts[i]) + 0.5 * south
  angle <- 2 * pi * phase
  cbind(sin1 = sin(angle), cos1 = cos(angle), sin2 = sin(2*angle), cos2 = cos(2*angle))
}

fit_fe <- function(y, X, fe, cluster, weights = NULL, cluster2 = NULL) {
  X <- as.matrix(X)
  stopifnot(nrow(X) == length(y), length(cluster) == length(y))
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
  if (anyDuplicated(colnames(X))) stop("The model columns must have unique names.")
  if (is.null(weights)) weights <- rep(1, length(y))
  keep <- is.finite(y) & complete.cases(X) & is.finite(weights) & weights > 0 & !is.na(cluster)
  if (!is.null(fe)) keep <- keep & !is.na(fe)
  if (!is.null(cluster2)) keep <- keep & !is.na(cluster2)
  x <- X[keep, , drop = FALSE]
  outcome <- y[keep]
  w <- weights[keep]
  runner <- cluster[keep]
  n <- length(outcome)
  if (!is.null(fe)) {
    group <- match(fe[keep], unique(fe[keep]))
    total_weight <- as.numeric(rowsum(w, group, reorder = FALSE))
    x_mean <- rowsum(x * w, group, reorder = FALSE) / total_weight
    y_mean <- as.numeric(rowsum(outcome * w, group, reorder = FALSE)) / total_weight
    x <- x - x_mean[group, , drop = FALSE]
    outcome <- outcome - y_mean[group]
  }
  x <- x * sqrt(w)
  outcome <- outcome * sqrt(w)
  # Scale the predictors before fitting, then return coefficients in their original units.
  norm <- sqrt(colSums(x^2))
  active <- which(norm > sqrt(.Machine$double.eps))
  z <- sweep(x[,active,drop = FALSE], 2L, norm[active], "/")
  q <- qr(z, tol = 1e-10)
  active <- active[q$pivot[seq_len(q$rank)]]
  z <- sweep(x[,active,drop = FALSE], 2L, norm[active], "/")
  q <- qr(z)
  coef_z <- qr.coef(q, outcome)
  residual <- as.numeric(outcome - z %*% coef_z)
  bread <- chol2inv(chol(crossprod(z)))
  scores <- z * residual
  k <- length(active)
  beta <- setNames(numeric(ncol(X)), colnames(X))
  beta[active] <- coef_z / norm[active]

  influence <- function(id) {
    score <- rowsum(scores, as.character(id), reorder = FALSE) %*% bread
    g <- nrow(score)
    stopifnot(g > 1L, n > k)
    score <- sweep(score, 2L, norm[active], "/") * sqrt((g/(g-1)) * ((n-1)/(n-k)))
    full <- matrix(0, g, length(beta), dimnames = list(rownames(score), names(beta)))
    full[, active] <- score
    full
  }
  score <- influence(runner)
  V <- crossprod(score)
  df <- nrow(score) - 1L
  if (!is.null(cluster2)) {
    other <- cluster2[keep]
    both <- paste(runner, other, sep = "|")
    V <- V + crossprod(influence(other)) - crossprod(influence(both))
    df <- min(df, length(unique(other))-1L)
  }
  list(beta = beta, vcov = V, df = df, keep = keep, score = score,
       n = n, groups = nrow(score))
}

lin_contrast <- function(fit, c, transform = "identity") {
  if (!is.null(names(c))) {
    if (!all(names(c) %in% names(fit$beta))) stop("The contrast contains an unknown column.")
    contrast <- setNames(numeric(length(fit$beta)), names(fit$beta))
    contrast[names(c)] <- c
  } else {
    stopifnot(length(c) == length(fit$beta))
    contrast <- c
  }
  effect <- sum(contrast * fit$beta)
  se <- sqrt(max(0, drop(crossprod(contrast, fit$vcov %*% contrast))))
  delta <- qt(.975, df = fit$df) * se
  values <- c(effect, effect - delta, effect + delta)
  if (transform == "percent") values <- 100 * expm1(values)
  if (transform == "points") values <- 100 * values
  p <- if (se > 0) 2 * pt(-abs(effect/se), fit$df) else NA_real_
  data.table::data.table(estimate = values[1L], lower = values[2L], upper = values[3L], p = p)
}
