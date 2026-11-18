# Compare published aggregates only after estimating the models from session records.
# Source-level checks are separate from the manuscript's rounded display values.
hr_equal <- function(actual, published, keys, fields, source, tolerance = 1e-6,
                     simulation = character()) {
    actual <- as.data.table(actual)
    published <- as.data.table(published)
    stopifnot(!anyDuplicated(actual[, ..keys]), !anyDuplicated(published[, ..keys]))
    index <- published[actual, on = keys, which = TRUE]
    stopifnot(length(index) == nrow(actual), !anyNA(index))
    rows <- lapply(names(fields), function(field_name) {
        simulated <- field_name %in% simulation
        answer <- copy(actual[, ..keys])
        answer[, quantity := do.call(paste, c(.SD, sep = " / ")), .SDcols = keys]
        answer <- answer[, .(quantity)]
        answer[, field := field_name]
        answer[, r_value := as.numeric(actual[[field_name]])]
        answer[, reference := as.numeric(published[[fields[[field_name]]]][index])]
        answer[, difference := r_value - reference]
        answer[, tolerance := if (simulated) NA_real_ else tolerance]
        answer[, status := fifelse(field %in% simulation, "different_rng",
            fifelse(is.finite(difference) & abs(difference) <= tolerance, "match", "different"))]
        answer[, source_file := source]
        answer
    })
    rbindlist(rows)
}
paired_ref <- "references/heart_rate_and_history"
sampling_ref <- "references/additional"
model_ref <- "references"
hr_fields <- c(estimate="estimate", lower="ci_low", upper="ci_high", p="p_value")
hr_check <- list()
for (item in list(c("heart_paired.csv", paired_ref, "paired_differences.csv"),
                  c("heart_chain.csv", paired_ref, "chain_delta.csv"),
                  c("heart_time_weighted_pairs.csv", paired_ref, "time_weighted_paired.csv"))) {
    file <- file.path(out, item[1])
    src <- file.path(item[2], item[3])
    hr_check[[src]] <- hr_equal(fread(file), fread(input_file(src)),
        c("quantity","target_label"), c(hr_fields,se="se"), src)
}
for (item in list(c("heart_decomposition.csv","main/heart_rate_decomposition.json"),
                  c("heart_hour_decomposition.csv","sensitivity/benchmark_heart_rate_decomposition.json"),
                  c("heart_month_decomposition.csv","sensitivity/runner_month_heart_rate_decomposition.json"),
                  c("heart_edge_decomposition.csv","edge_heart_rate_decomposition.json"))) {
    src <- file.path(if (item[1]=="heart_edge_decomposition.csv") paired_ref else model_ref,item[2])
    published <- jsonlite::fromJSON(input_file(src),simplifyVector=FALSE)
    ref <- rbindlist(lapply(c("P90","P95"), function(lab) {
        items <- Filter(is.list,published[[lab]])
        rbindlist(lapply(names(items), function(q)
            as.data.table(items[[q]])[, ':='(quantity=q,target_label=lab)]))
    }))
    ans <- fread(file.path(out,item[1]))
    # Two additional quantities from the original benchmark are saved in another file.
    ans <- ans[quantity %in% ref$quantity]
    stopifnot(nrow(ans)==nrow(ref))
    hr_check[[src]] <- hr_equal(ans,ref,c("quantity","target_label"),c(hr_fields[1:3],se="se"),src)
}
for (item in list(c("heart_time_weighted.csv",paired_ref,"time_weighted_contrasts.csv"),
                  c("heart_recording_resolution.csv",sampling_ref,"heart_rate_sampling_resolution.csv"),
                  c("heart_resolution_groups.csv",sampling_ref,"sampling_resolution_groups.csv"),
                  c("heart_place_checks.csv",model_ref,"sensitivity/place_contrasts.csv"),
                  c("heart_month_checks.csv",model_ref,"sensitivity/runner_month_contrasts.csv"),
                  c("heart_edge_checks.csv",paired_ref,"edge_contrasts.csv"),
                  c("heart_benchmark.csv",model_ref,"sensitivity/benchmark_contrasts.csv"))) {
    src <- file.path(item[2],item[3])
    keys <- c("model_id","target_label", if (item[1]=="heart_resolution_groups.csv") "group_level")
    hr_check[[src]] <- hr_equal(fread(file.path(out,item[1])),fread(input_file(src)),
        keys,hr_fields,src)
}
src <- file.path(model_ref,"main/heart_rate_decomposition_curves.csv")
ref <- fread(input_file(src))
ans <- fread(file.path(out,"heart_decomposition_curves.csv"))[quantity %in% ref$quantity]
ref[, exposure_c:=round(exposure_c,8)]
ans[, exposure_c:=round(exposure_c,8)]
stopifnot(nrow(ref)==nrow(ans))
hr_check[[src]] <- hr_equal(ans,ref,c("exposure_c","quantity"),c(hr_fields[1:3],se="se"),src)
if (file.exists(file.path(out,"heart_resolution_tests.csv"))) {
    src <- file.path(sampling_ref,"sampling_resolution_tests.csv")
    hr_check[[src]] <- hr_equal(fread(file.path(out,"heart_resolution_tests.csv")),
        fread(input_file(src)),"model_id",c(F_statistic="F_statistic",p="p_value"),src)
}
if (file.exists(file.path(out,"heart_lee.csv"))) {
    src <- file.path(paired_ref,"lee_bounds.csv")
    ans <- fread(file.path(out,"heart_lee.csv"))
    keys <- c("model_id","trimming_basis")
    vars <- setdiff(names(ans),keys)
    hr_check[[src]] <- hr_equal(ans, fread(input_file(src)),keys,setNames(vars,vars),src,
        simulation=c("lower_bound_2_5_percentile","upper_bound_97_5_percentile"))
    if (nzchar(Sys.getenv("FITREC_PYTHON"))) {
        exact_fields <- c("lower_bound_2_5_percentile","upper_bound_97_5_percentile")
        hr_check[[paste(src,"numpy")]] <- hr_equal(fread(file.path(out,"heart_lee_numpy.csv")),
            fread(input_file(src)),keys,setNames(exact_fields,exact_fields),src)[,rng:="paper"]
    }
}
if (file.exists(file.path(out,"heart_bootstrap.csv"))) {
    src <- file.path(paired_ref,"paired_model_bootstrap.json")
    paper_boot <- jsonlite::fromJSON(input_file(src),simplifyVector=FALSE)$estimates
    ref <- rbindlist(lapply(names(paper_boot),function(q) as.data.table(paper_boot[[q]])[,quantity:=q]))
    ans <- fread(file.path(out,"heart_bootstrap.csv"))
    stopifnot(setequal(ans$quantity,ref$quantity))
    vars <- c("estimate","bootstrap_se","ci_low","ci_high","replicates")
    hr_check[[src]] <- hr_equal(ans,ref,"quantity",setNames(vars,vars),src,
        simulation=c("bootstrap_se","ci_low","ci_high"))
    if (nzchar(Sys.getenv("FITREC_PYTHON"))) {
        exact_fields <- c("bootstrap_se","ci_low","ci_high")
        hr_check[[paste(src,"numpy")]] <- hr_equal(fread(file.path(out,"heart_bootstrap_numpy.csv")),
            ref,"quantity",setNames(exact_fields,exact_fields),src)[,rng:="paper"]
    }
}
hr_verification <- rbindlist(hr_check,fill=TRUE)
if (!"rng" %in% names(hr_verification)) hr_verification[,rng:=NA_character_]
hr_verification[status=="different_rng",rng:="independent_R"]
fwrite(hr_verification,file.path(out,"heart_verification.csv"))
message("Heart-rate comparisons: ", sum(hr_verification$status=="match"), " match, ",
    sum(hr_verification$status=="different"), " different, ",
    sum(hr_verification$status=="different_rng"), " bootstrap quantities from independent draws")
