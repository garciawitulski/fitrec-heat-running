# Rebuild the manuscript and supplementary tables from the printed cells.
# Values shown as lower bounds retain that notation. The temperature models
# separately reproduce the Table 2 contrasts.
local({
  published <- file.path(repo, "data", "published")
  catalog <- fread(file.path(published, "table_catalog.csv"), colClasses = "character", encoding = "UTF-8")
  cells <- fread(file.path(published, "table_cells.csv"), colClasses = c(table_id="character",text="character"),
                 na.strings = NULL, encoding = "UTF-8")
  stopifnot(nrow(catalog)==31L, !anyDuplicated(catalog$table_id),
            setequal(catalog$table_id, c("1","2",paste0("S",1:29))),
            setequal(catalog$table_id,cells$table_id))

  escape <- function(s) {
    s <- gsub("&", "&amp;", s, fixed=TRUE)
    s <- gsub("<", "&lt;", s, fixed=TRUE)
    gsub(">", "&gt;", s, fixed=TRUE)
  }
  style <- paste(
    "body{font:16px/1.45 system-ui,sans-serif;color:#192a30;max-width:1050px;margin:32px auto;padding:0 16px}",
    "h1{font-size:1.35rem}h2{font-size:1.05rem}table{border-collapse:collapse;width:100%;margin:.8em 0}",
    "td,th{padding:5px 8px;vertical-align:top;text-align:right}td:first-child,th:first-child{text-align:left}",
    "thead{border-top:2px solid #293d43;border-bottom:1px solid #758184}",
    "tbody{border-bottom:2px solid #293d43}.section th{text-align:left;background:#edf2f2}",
    ".note{font-size:.9rem}.scroll{overflow-x:auto}@media print{body{font-size:11px;margin:0}}")
  for (i in seq_len(nrow(catalog))) {
    id <- catalog$table_id[i]
    destination <- file.path(public_out,"tables",if(startsWith(id,"S")) "supplementary" else "manuscript")
    dir.create(destination,recursive=TRUE,showWarnings=FALSE)
    tab <- cells[table_id==id]
    stopifnot(!anyDuplicated(tab,by=c("row","column")),all(tab$span>=1L))
    title <- paste0("Table ",id,". ",catalog$title[i])
    page <- c("<!doctype html><html lang='en'><head><meta charset='utf-8'>",
              "<meta name='viewport' content='width=device-width, initial-scale=1'>",
              paste0("<title>",escape(title),"</title><style>",style,"</style></head><body>"),
              paste0("<h1>",escape(title),"</h1>"))
    for (panel_number in unique(tab$panel_id)) {
      p <- tab[panel_id==panel_number]
      stopifnot(nrow(p)>0L, p$role[1]=="header", any(p$role=="data"))
      if (nzchar(p$panel[1])) page <- c(page,paste0("<h2>",escape(p$panel[1]),"</h2>"))
      page <- c(page,"<div class='scroll'><table><thead>")
      ended_header <- FALSE
      for (j in unique(p$row)) {
        row <- p[row==j][order(column)]
        kind <- row$role[1]
        if (!ended_header && kind!="header") {
          page <- c(page,"</thead><tbody>")
          ended_header <- TRUE
        }
        tag <- if (kind=="data") "td" else "th"
        scope <- if (kind=="header") ifelse(row$span>1," scope='colgroup'"," scope='col'") else ""
        html_cells <- paste0("<",tag,scope," colspan='",row$span,"'>",escape(row$text),"</",tag,">")
        page <- c(page,paste0("<tr class='",kind,"'>",paste(html_cells,collapse=""),"</tr>"))
      }
      page <- c(page,"</tbody></table></div>")
    }
    page <- c(page,paste0("<p class='note'><strong>Note.</strong> ",escape(catalog$note[i]),
                          "</p></body></html>"))
    writeLines(enc2utf8(page),file.path(destination,paste0("table_",id,".html")),useBytes=TRUE)
    # The CSV retains each panel's original column headings as rows.
    wide <- dcast(tab,panel_id+panel+row+role~column,value.var="text",fill="")
    setorder(wide,row)
    fwrite(wide,file.path(destination,paste0("table_",id,".csv")))
  }
})
