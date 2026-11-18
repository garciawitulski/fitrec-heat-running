# Paired analysis of the same runs, using whole-cohort cross-products.
# Each block has two temperature spline terms after adjustment for the other
# predictors. The joint score matrix accounts for the correlation between
# the segment and whole-session estimates from the same runners.
local({
  base <- file.path(repo,"data","processed","paired_model")
  design <- fread(file.path(base,"design.csv"))
  outcomes <- fread(file.path(base,"outcomes.csv"))
  products <- fread(file.path(base,"joint_score_products.csv"))
  contrasts <- fread(file.path(base,"temperature_contrasts.csv"))
  fits <- c("fixed","free","whole")
  beta <- numeric(6L)
  bread <- matrix(0,6L,6L)
  for (j in seq_along(fits)) {
    ix <- (2L*j-1L):(2L*j)
    A <- matrix(design[model==fits[j]][order(column,row),value],2,2)
    b <- outcomes[model==fits[j]][order(row),value]
    beta[ix] <- solve(A,b)
    bread[ix,ix] <- solve(A)
  }
  meat <- matrix(products[order(column,row),value],6L,6L)
  V <- bread %*% meat %*% bread
  weights <- rbind(steady_fixed_speed=c(1,0,0),steady_speed_free=c(0,1,0),
    whole_session_mean=c(0,0,1),steady_speed_free_minus_whole_session=c(0,1,-1),
    steady_fixed_speed_minus_whole_session=c(1,0,-1),fixed_speed_minus_speed_free=c(1,-1,0))
  rows <- list()
  for (i in seq_len(nrow(contrasts))) for (j in seq_len(nrow(weights))) {
    v <- kronecker(weights[j,],c(contrasts$linear[i],contrasts$cubic[i]))
    estimate <- as.numeric(crossprod(v,beta))
    se <- sqrt(as.numeric(crossprod(v,V %*% v)))
    # The paper gives the runner count rounded down to 660. Use the two
    # endpoints of the possible t degrees of freedom (659 to 668).
    wide <- qt(.975,659)*se
    narrow <- qt(.975,668)*se
    rows[[length(rows)+1L]] <- data.table(quantity=rownames(weights)[j],target=contrasts$target[i],
      estimate,se,lower_conservative=estimate-wide,upper_conservative=estimate+wide,
      lower_at_max_df=estimate-narrow,upper_at_max_df=estimate+narrow,
      p_at_max_df=2*pt(-abs(estimate/se),668),p_conservative=2*pt(-abs(estimate/se),659))
  }
  dest <- file.path(public_out,"models")
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  fwrite(rbindlist(rows),file.path(dest,"paired_heart_rate.csv"))
})
