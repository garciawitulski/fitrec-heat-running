# Figure 3, Tables S5 and S18. Resample runners jointly across all outcome models.
# Duplicating a runner duplicates the same within-runner cross-products. Computing
# them once gives the refitted fixed-effects model for each bootstrap draw.
# Native R samples have the same design and replicate count as the Python analysis.
set.seed(20261002)
ids <- sort(unique(c(hr_main$user_hash,hr_phys$user_hash)))
G <- length(ids)
B <- 999L
weights_boot <- cbind(rep(1L,G),replicate(B,tabulate(sample.int(G,G,replace=TRUE),G)))
python_rng <- Sys.getenv("FITREC_PYTHON")
plan_names <- "r"
# To replay the original draws, pass Python the sample sizes and seeds. Python
# returns selection frequencies. R reads the temporary file and then deletes it.
read_numpy_plan <- function(design,sizes,replicates) {
    tmp <- tempfile()
    on.exit(unlink(tmp))
    status <- system2(python_rng,c(shQuote(file.path(repo,"record_level","scripts","bootstrap_indices.py")),
                                   shQuote(tmp),design,as.character(sizes)))
    if (status!=0) stop("FITREC_PYTHON must point to Python with NumPy installed.")
    count <- as.integer(sum(sizes)*replicates)
    values <- readBin(tmp,integer(),n=count,size=4,endian="little")
    stopifnot(length(values)==count,file.size(tmp)==count*4)
    values
}
if (nzchar(python_rng)) {
    numeric_plan <- matrix(read_numpy_plan("joint",G,B),nrow=G)
    stopifnot(all(colSums(numeric_plan)==G))
    weights_boot <- cbind(weights_boot,numeric_plan)
    plan_names <- c("r","numpy")
}

hr_boot_fit <- function(X,y,runner,model) {
    ok <- is.finite(y) & complete.cases(X)
    runner <- runner[ok]
    X <- X[ok,,drop=FALSE]
    y <- y[ok]
    p <- ncol(X)
    ix <- split(seq_along(runner),factor(runner,levels=ids),drop=TRUE)
    product <- matrix(0,G,p*p+p)
    for (id in names(ix)) {
        rows <- ix[[id]]
        z <- sweep(X[rows,,drop=FALSE],2L,colMeans(X[rows,,drop=FALSE]))
        product[match(id,ids),] <- c(crossprod(z),crossprod(z,y[rows]-mean(y[rows])))
    }
    combined <- crossprod(product,weights_boot)
    keep <- intersect(c("T1","T2","V1","V2"),colnames(X))
    beta <- matrix(0,nrow=length(keep),ncol=ncol(weights_boot),dimnames=list(keep,NULL))
    for (j in seq_len(ncol(weights_boot))) {
        xx <- matrix(combined[seq_len(p*p),j],p,p)
        xy <- combined[p*p+seq_len(p),j]
        active <- which(diag(xx)>1e-18)
        norm <- sqrt(diag(xx)[active])
        gram <- xx[active,active,drop=FALSE]/outer(norm,norm)
        rhs <- xy[active]/norm
        # A start hour can disappear in a draw. Keep the estimable column space.
        e <- eigen(gram,symmetric=TRUE)
        estimable <- e$values > max(e$values)*1e-12
        Q <- e$vectors[,estimable,drop=FALSE]
        solution <- setNames(numeric(p),colnames(X))
        solution[active] <- as.vector(Q %*% (crossprod(Q,rhs)/e$values[estimable]))/norm
        beta[,j] <- solution[keep]
    }
    stopifnot(max(abs(beta[,1]-model$beta[keep]))<1e-7)
    beta
}
boot_fits <- list(
    gap=hr_boot_fit(hr_xmain,hr_main$heat_avoided_c,hr_main$user_hash,hr_gap),
    hour=hr_boot_fit(hr_xhour,hr_main$heat_avoided_c,hr_main$user_hash,hr_hour),
    speed=hr_boot_fit(hr_xmain,hr_main$gps_speed_median_m_s,hr_main$user_hash,hr_speed),
    net=hr_boot_fit(hr_xmain,hr_main$hr_mean,hr_main$user_hash,hr_net),
    phys=hr_boot_fit(hr_xphys,hr_phys$hr_lag60,hr_phys$user_hash,hr_phfit),
    fixed=hr_boot_fit(hr_xfixed,hr_day$hr_lag60,hr_day$user_hash,hr_fixed),
    free=hr_boot_fit(hr_xfree,hr_day$hr_lag60,hr_day$user_hash,hr_free),
    whole=hr_boot_fit(hr_xfree,hr_day$session_hr_mean,hr_day$user_hash,hr_whole))
for (id in names(boot_fits)) rownames(boot_fits[[id]]) <- paste(id,rownames(boot_fits[[id]]),sep=".")
betas <- do.call(rbind,boot_fits)

# Each draw also recalculates the mean thermal gap and the mean temperature basis.
sums <- rowsum(cbind(n=1,gap=hr_main$heat_avoided_c,hr_xmain[,c("T1","T2")]),hr_main$user_hash)
averages <- crossprod(sums,weights_boot[match(rownames(sums),ids),,drop=FALSE])
averages <- sweep(averages,2L,averages["n",],"/")
basis_ref <- as.numeric(rcs_basis(hr_targets[["P50"]],hr_kday))
hr_boot <- list()
for (label in c("P90","P95")) {
    target <- hr_targets[[label]]
    v <- hr_spline_change(target,hr_targets[["P50"]],hr_kday)
    draws <- lapply(seq_len(ncol(weights_boot)),function(j) {
        b <- betas[,j]
        parts <- hr_parts(b,target,unname(averages["gap",j]),basis_ref-averages[c("T1","T2"),j])
        parts <- parts[!names(parts)%in%c("session_temperature_on_reference_day_c","speed_change_m_s")]
        fixed <- sum(v*b[c("fixed.T1","fixed.T2")])
        free <- sum(v*b[c("free.T1","free.T2")])
        whole <- sum(v*b[c("whole.T1","whole.T2")])
        c(parts,steady_fixed_speed=fixed,steady_speed_free=free,whole_session_common_sessions=whole,
          link_transfer=parts[["rise_with_observed_timing"]]-fixed,link_speed=fixed-free,
          link_measure=free-whole,link_sample_and_covariates=whole-parts[["net_heart_rate_change"]],
          steady_fixed_speed_minus_whole_session=fixed-whole)
    })
    draws <- do.call(cbind,draws)
    stopifnot(all(is.finite(draws)))
    for (k in seq_along(plan_names)) {
        simulations <- draws[,1L+(k-1L)*B+seq_len(B),drop=FALSE]
        hr_boot[[paste(label,k)]] <- data.table(quantity=paste(rownames(draws),label,sep="@"),estimate=draws[,1],
            bootstrap_se=apply(simulations,1,sd),ci_low=apply(simulations,1,quantile,probs=.025),
            ci_high=apply(simulations,1,quantile,probs=.975),replicates=B,plan=plan_names[k])
    }
}
bootstrap_tables <- rbindlist(hr_boot)
for (engine in plan_names) {
    name <- if (engine=="r") "heart_bootstrap.csv" else "heart_bootstrap_numpy.csv"
    fwrite(bootstrap_tables[plan==engine,!"plan"],file.path(out,name))
}
message("Runner cluster bootstrap: ",B," complete draws at both temperature contrasts")

# Table S26. On the P50 window recover the participation rate and fit the registered
# participation model with runner-by-month effects, weekday controls and CR1.
panel <- read_fitrec("panels/runner_days.parquet",
                    c("user_hash","run","day_max_c","local_weekday","year_month_local"))
panel <- panel[is.finite(day_max_c)]
q_ref <- panel[abs(day_max_c-hr_targets[["P50"]])<=2,mean(as.numeric(run))]
panel_x <- hr_design(panel,season=FALSE)
panel_fit <- fit_fe(as.numeric(panel$run),panel_x,paste(panel$user_hash,panel$year_month_local),panel$user_hash)
selection <- hr_contrasts(panel_fit,"p_run")[target_label=="P95"]
trim_shares <- c(point_estimate=max(0,-selection$estimate)/q_ref,
                 lower_95_limit=max(0,-selection$lower)/q_ref)
stopifnot(all(trim_shares>=0 & trim_shares<1))
set.seed(20261104)
hr_bounds <- list()
if (nzchar(python_rng)) {
    sizes <- c(uniqueN(hr_main$user_hash[hr_valid]),
               uniqueN(hr_main$user_hash[is.finite(hr_main$gps_speed_median_m_s)]))
    lee_plan <- read_numpy_plan("lee",sizes,2L*199L)
    offset <- 0L
}
for (id in c("h1_hr_mean","b2_speed")) {
    is_hr <- id=="h1_hr_mean"
    fit <- if (is_hr) hr_net else hr_speed
    ok <- if (is_hr) hr_valid else is.finite(hr_main$gps_speed_median_m_s)
    X <- hr_xmain[ok,,drop=FALSE]
    y <- if (is_hr) hr_main$hr_mean[ok] else hr_main$gps_speed_median_m_s[ok]
    d <- data.table(user=hr_main$user_hash[ok],exposure=hr_main$day_max_c[ok],
                    adjusted=y-as.vector(X%*%fit$beta))
    d[, adjusted:=adjusted-mean(adjusted),by=user]
    d[, adjusted:=adjusted+as.vector(X[,c("T1","T2")]%*%fit$beta[c("T1","T2")])]
    roster <- sort(unique(d$user))
    # Sort the reference outcome once. Repeated draws are frequency weights.
    ref <- d[abs(exposure-hr_targets[["P50"]])<=2][order(adjusted)]
    hot <- d[abs(exposure-hr_targets[["P95"]])<=2]
    jref <- match(ref$user,roster)
    jhot <- match(hot$user,roster)
    spline <- hr_contrasts(fit,id)[target_label=="P95"]
    for (method in names(trim_shares)) {
        share <- trim_shares[[method]]
        bound <- function(w) {
            weight <- w[jref]
            n <- sum(weight)
            k <- floor(share*n)
            # Remove exactly k sessions. A duplicated outcome can straddle the trim.
            left <- pmin(weight,pmax(0,cumsum(weight)-k))
            right <- pmin(weight,pmax(0,n-k-(cumsum(weight)-weight)))
            target_mean <- weighted.mean(hot$adjusted,w[jhot])
            c(lower=target_mean-weighted.mean(ref$adjusted,left),
              upper=target_mean-weighted.mean(ref$adjusted,right))
        }
        n <- length(roster)
        point <- bound(rep(1,n))
        draw <- replicate(199L,bound(tabulate(sample.int(n,n,replace=TRUE),n)))
        row <- data.table(model_id=id,trimming_basis=method,
            q_reference=q_ref,selection_change_p_run_p95=selection$estimate,selection_change_lower_95=selection$lower,
            spline_contrast_p95=spline$estimate,spline_ci_low=spline$lower,spline_ci_high=spline$upper,
            windowed_difference=mean(hot$adjusted)-mean(ref$adjusted),trimming_share=share,
            lower_bound=point[1],upper_bound=point[2],
            lower_bound_2_5_percentile=unname(quantile(draw[1,],.025)),
            upper_bound_97_5_percentile=unname(quantile(draw[2,],.975)),bootstrap_replicates=199L,plan="r")
        hr_bounds[[length(hr_bounds)+1L]] <- row
        if (nzchar(python_rng)) {
            draw_weights <- matrix(lee_plan[offset+seq_len(199L*n)],nrow=n)
            offset <- offset+199L*n
            stopifnot(all(colSums(draw_weights)==n))
            exact <- apply(draw_weights,2L,bound)
            hr_bounds[[length(hr_bounds)+1L]] <- copy(row)[, ':='(plan="numpy",
                lower_bound_2_5_percentile=unname(quantile(exact[1,],.025)),
                upper_bound_97_5_percentile=unname(quantile(exact[2,],.975)))]
        }
    }
}
bounds_tables <- rbindlist(hr_bounds)
for (engine in plan_names) {
    name <- if (engine=="r") "heart_lee.csv" else "heart_lee_numpy.csv"
    fwrite(bounds_tables[plan==engine,!"plan"],file.path(out,name))
}
