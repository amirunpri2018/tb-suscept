#' @importFrom magrittr %>%
#' @export
download_tss <- function(gene_names,
                         archive = "dec2015.archive.ensembl.org",
                         tss_all_fname = NULL,
                         overwrite = FALSE
                         ) {
  # Create connection to Ensembl database
  ensembl <- biomaRt::useMart(host = archive,
                              biomart = "ENSEMBL_MART_ENSEMBL",
                              dataset = "hsapiens_gene_ensembl")
  # Download TSS
  if (file.exists(tss_all_fname) & !overwrite) {
    tss_all <- readRDS(tss_all_fname)
  } else {
    tss_all <- biomaRt::getBM(attributes = c("ensembl_gene_id", "chromosome_name",
                                             "transcription_start_site", "strand"),
                              filters = "ensembl_gene_id",
                              values = gene_names,
                              mart = ensembl)
    if (!is.null(tss_all_fname)) {
      saveRDS(tss_all, file = tss_all_fname)
    }
  }
  # One tss per gene
  tss <- tss_all %>%
    dplyr::group_by_(~ensembl_gene_id) %>%
    dplyr::summarize_(chr = ~chromosome_name[1],
               strand = ~strand[1],
               tss = ~if (strand == 1) min(transcription_start_site) else max(transcription_start_site))
  stopifnot(nrow(tss) == length(gene_names))
  return(tss)
}

#' @export
define_gene_windows <- function(x,
                                window_size
                                ) {
  stopifnot(is.data.frame(x),
            c("ensembl_gene_id", "chr", "strand", "tss") %in% colnames(x),
            is.numeric(window_size),
            length(window_size) == 1)
  x$start <- x$tss - window_size
  x$end <- x$tss + window_size
  x$strand <- ifelse(x$strand == 1, "+", "-")
  result <- GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
  GenomeInfoDb::seqlevels(result) <- paste0("ch", GenomeInfoDb::seqlevels(result))
  stopifnot(class(result) == "GRanges")
  return(result)
}

#' @import SNPlocs.Hsapiens.dbSNP144.GRCh38
#' @export
obtain_snp_coords <- function(rsid,
                              snp_coords_fname,
                              overwrite = FALSE
                              ) {
  if (file.exists(snp_coords_fname)) {
    snp_coords <- readRDS(snp_coords_fname)
  } else {
    snp_coords <- BSgenome::snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh38,
                                     rsid,
                                     ifnotfound = "drop")
    if (!is.null(snp_coords_fname)) {
      saveRDS(snp_coords, file = snp_coords_fname)
    }
  }
  stopifnot(S4Vectors::mcols(snp_coords)$RefSNP_id %in% rsid)
  return(snp_coords)
}

#' @export
combine_snps_and_genes <- function(snp_coords,
                                   gene_windows,
                                   pval
                                   ) {
  # Overlap snps and genes
  overlaps <- GenomicRanges::findOverlaps(snp_coords, gene_windows, ignore.strand = TRUE)

  # Add back names
  results <- data.frame(S4Vectors::as.matrix(overlaps))
  colnames(results) <- c("rsID", "gene")
  results$rsID <- S4Vectors::mcols(snp_coords)$RefSNP_id[results$rsID]
  results$gene <- S4Vectors::mcols(gene_windows)$ensembl_gene_id[results$gene]

  return(results)
}

#' @export
add_gwas_pval <- function(x,
                          rsid,
                          pval
                          ) {
  stopifnot(length(rsid) == length(pval))
  names(pval) <- rsid
  x$gwas_p <- pval[x$rsID]
  stopifnot(!is.na(x$gwas_p))
  return(x)
}

#' @importFrom magrittr %>%
#' @export
assign_min_pval <- function(x
                            ) {
  stopifnot(c("rsID", "gene", "gwas_p") %in% colnames(x))

  # Assign minimum
  result <- x %>%
    dplyr::group_by_(~gene) %>%
    dplyr::summarize_(gwas_p = ~min(gwas_p),
                     n_snps = ~n())
  return(result)
}

#' @export
run_gwas_enrich <- function(gene_names,
                            archive = "dec2015.archive.ensembl.org",
                            tss_all_fname = NULL,
                            overwrite = FALSE,
                            window_size,
                            rsid,
                            snp_coords_fname,
                            pval
                            ) {
  # Download the most upstream TSS for each gene
  tss <- download_tss(gene_names = gene_names,
                    archive = archive,
                    tss_all_fname = tss_all_fname,
                    overwrite = overwrite)
  tss_gr <- define_gene_windows(tss, window_size = window_size)
  # Obtain SNP coordinates
  snp_coords <- obtain_snp_coords(rsid = rsid,
                                  snp_coords_fname = snp_coords_fname,
                                  overwrite = overwrite)
  snp_genes <- combine_snps_and_genes(snp_coords, tss_gr)
  snp_genes_pval <- add_gwas_pval(snp_genes, rsid, pval)
  snp_genes_pval_unique <- assign_min_pval(snp_genes_pval)

  return(snp_genes_pval_unique)
}
