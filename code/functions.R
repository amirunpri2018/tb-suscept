
#' Save ExpressionSet as rds file
download_thuong <- function(fname = "GSE11199.rds", overwrite = FALSE) {
  if (!file.exists(geo_fname) | overwrite) {
    gse <- GEOquery::getGEO(GEO = "GSE11199")
    gse <- gse[[1]]
    saveRDS(gse, fname)
  }
}
