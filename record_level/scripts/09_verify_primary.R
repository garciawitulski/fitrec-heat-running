# Compare values at the precision reported in the paper.
ref <- fread(file.path(repo,"data","processed","primary_printed.csv"), colClasses=c(printed="character",model="character",outcome="character"))
t1 <- fread(file.path(out,"table1.csv"))[,.(table="T1",model="main",outcome=as.character(row),target="sample",statistic,value)]
t2 <- fread(file.path(out,"table2.csv"))[,':='(table="T2",model="main")]
blocks <- list(t2)
for (fname in c("supp_session.csv","supp_registered.csv")) {
  if (file.exists(file.path(out,fname))) {
    block <- fread(file.path(out,fname),colClasses=c(model="character"))
    blocks[[length(blocks)+1L]] <- block
  }
}
wide <- rbindlist(blocks,fill=TRUE)
keys <- c("table","model","outcome","target")
stats <- intersect(unique(ref$statistic),names(wide))
wide[,(stats):=lapply(.SD,as.numeric),.SDcols=stats]
long <- melt(wide,id.vars=keys,measure.vars=stats,
             variable.name="statistic",variable.factor=FALSE)
values <- rbind(t1,long)
checked <- merge(ref,values,by=c(keys,"statistic"),all.x=TRUE,sort=FALSE)
stopifnot(nrow(checked)==nrow(ref))
checked[statistic %in% c("n","n_sessions","n_runners"),tolerance:=0]
checked[,abs_difference:=abs(value-reference)]
checked[,rel_difference:=fifelse(reference!=0,abs_difference/abs(reference),NA_real_)]
checked[,status:=fcase(is.na(value),"NOT REPLICABLE",
                       upper_bound & value<reference,"MATCH",
                       !upper_bound & abs_difference<=tolerance+1e-12,"MATCH",
                       default="MISMATCH")]
conventions <- fread(file.path(repo,"data","processed","primary_conventions.csv"))
checked <- merge(checked,conventions,by=c(keys,"statistic"),all.x=TRUE,sort=FALSE)
checked[status=="MATCH" & !is.na(note),status:="EXPLAINED"]
# Four values in the supplement were printed from a Python CSV rounded to four decimal places.
rounding <- fread(file.path(repo,"data","processed","primary_rounding.csv"))
checked <- merge(checked,rounding,by=c(keys,"statistic"),all.x=TRUE,sort=FALSE)
checked[!is.na(intermediate) & abs(value/unit_multiplier-intermediate)<=0.00005+1e-12 &
        python_display==reference,":="(status="EXPLAINED",
        note="The R estimate matches the CSV value to four decimal places. Formatting that value as in the original Python code reproduces the table entry.")]
fwrite(checked,file.path(out,"primary_comparison.csv"))
print(checked[,.N,by=.(table,status)][order(table,status)])




# Table 2. Compare the results before rounding with the original Python estimates,
# including standard errors on the model scale.
original <- fread(input_file("references/main/table2_contrasts.csv"))
original[outcome=="speed_lag60",outcome:="steady_speed"]
original <- original[,.(outcome,target=target_label,estimate=display_estimate,lower=display_ci_low,
                        upper=display_ci_high,se_model=se,p=p_value,p_holm=p_holm_primary_family)]
contrast_keys <- c("outcome","target")
fields <- c("estimate","lower","upper","se_model","p","p_holm")
archived_values <- melt(original,id.vars=contrast_keys,measure.vars=fields,
                  variable.name="statistic",value.name="reference")[!is.na(reference)]
replicated_values <- melt(t2,id.vars=contrast_keys,measure.vars=fields,
                    variable.name="statistic",value.name="value")
full <- merge(archived_values,replicated_values,by=c(contrast_keys,"statistic"),all.x=TRUE)
stopifnot(nrow(full)==nrow(archived_values),!anyNA(full$value))
full[,abs_difference:=abs(value-reference)]
full[,rel_difference:=abs_difference/pmax(abs(reference),1e-300)]
# The public file stores some p-values to 16 decimal places. Allow for that
# rounding and use a relative tolerance of 1e-7 for the other values.
full[,tolerance:=pmax(ifelse(statistic %in% c("p","p_holm"),1e-16,1e-12),1e-7*abs(reference))]
full[,status:=fifelse(abs_difference<=tolerance,"MATCH","MISMATCH")]
fwrite(full,file.path(out,"check_02_table2_full.csv"))
print(full[,.N,by=status])
