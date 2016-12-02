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
  x$pval <- pval[x$rsID]
  stopifnot(!is.na(x$pval))
  return(x)
}

#' @importFrom magrittr %>%
#' @export
assign_min_pval <- function(x
                            ) {
  stopifnot(c("rsID", "gene", "pval") %in% colnames(x))

  # Assign minimum
  result <- x %>%
    dplyr::group_by_(~gene) %>%
    dplyr::summarize_(pval = ~min(pval),
                     n_snps = ~n())
  return(result)
}

#' @export
add_effect_size <- function(x,
                            effect_size
                            ) {
  stopifnot(x$gene %in% names(effect_size))
  x$effect_size <- effect_size[x$gene]
  if (any(x$effect_size < 0)) {
    message("Effect size contained negative values. Converted to absolute value.")
    x$effect_size <- abs(x$effect_size)
  }
  return(x)
}

#' Main enrichment function.
#'
#' For each interval, calculate the enrichment compared to the background enrichment.
#'
#' x - the variable that is split into intervals
#' y - the variable that is being tested for enrichment
#' cutoff - the value of y at which to separate genes when cacluting enrichment
#' m - the average number of genes in each interval
#' x_direction - With increasing values of x, should genes be included if they are greater than or less than
#'               the metric that defines that interval. Default "greater", anything else with do "lesser".
#' cutoff_direction - When testing for enrichment of y at the given interval of x, are genes counted if they are
#'                    greater or lesser than the cuottf.
#'
#' @export
enrich <- function(x, y, cutoff, m = 50,
                   x_direction = "greater",
                   cutoff_direction = "greater") {
  intervals <- Hmisc::cut2(x, m = m, onlycuts = TRUE)
  enrichment <- numeric(length = length(intervals))
  sizes <- numeric(length = length(intervals))
  if (cutoff_direction == "greater") {
    background_enrich <- sum(y > cutoff) / length(y)
  } else {
    background_enrich <- sum(y < cutoff) / length(y)
  }
  for (i in seq_along(intervals)) {
    if (x_direction == "greater") {
      y_sub <- y[x > intervals[i]]
    } else {
      y_sub <- y[x < intervals[i]]
    }
    sizes[i] <- length(y_sub)
    # browser()
    if (cutoff_direction == "greater") {
      enrichment[i] <- sum(y_sub > cutoff) / sizes[i] / background_enrich
    } else {
      enrichment[i] <- sum(y_sub < cutoff) / sizes[i] / background_enrich
    }
  }
  return(data.frame(enrichment = enrichment, intervals = intervals, sizes = sizes))
}

#' Calculates the enrichment of the actual data (returned in column 1) and the
#' enrichment of perumutations.
#'
#' iterations - The number of iterations to calculate enrichment with permuted data.
#'
#' See function enrich for descriptions of other arguments.
#'
#' @export
enrich_full <- function(x, y, cutoff, m = 50,
                        x_direction = "greater",
                        cutoff_direction = "greater",
                        iterations = 100) {
  result <- list()
  class(result) <- "enrichment"
  result$main <- enrich(x = x, y = y, cutoff = cutoff,
                        m = m, x_direction = x_direction,
                        cutoff_direction = cutoff_direction)
  mat_enrichment <- matrix(nrow = nrow(result$main), ncol = iterations)
  mat_intervals <- matrix(nrow = nrow(result$main), ncol = iterations)
  mat_sizes <- matrix(nrow = nrow(result$main), ncol = iterations)
#   mat_enrichment[, 1] <- main$enrichment
#   mat_intervals[, 1] <- main$intervals
#   mat_sizes[, 1] <- main$sizes
  # browser()
  for (iter in 1:iterations) {
    permuted <- enrich(x = sample(x), y = y, cutoff = cutoff,
                       m = m, x_direction = x_direction,
                       cutoff_direction = cutoff_direction)
    mat_enrichment[, iter] <- permuted$enrichment
    mat_intervals[, iter] <- permuted$intervals
    mat_sizes[, iter] <- permuted$sizes
  }
  result$permutations <- list(enrichment = mat_enrichment, intervals = mat_intervals,
                              sizes = mat_sizes)
  return(result)
}

print.enrichment <- function(x, ...) {
  str(x)
}

#' Function to create enrichment plot
#'
#' Input is output from enrich_full.
#'
#' enrichment - fold enrichment of a particular metric at increasing stringency
#'              of another metric
#' intervals - the cutoff value for the metric
#' sizes - the number of genes considered at each interval
#' ... - arguments passed to `plot`
#' @export
plot.enrichment <- function(x, ...) {
  enrichment <- x$main$enrichment
  ymin <- min(cbind(enrichment, x$permutations$enrichment), na.rm = TRUE)
  ymax <- max(cbind(enrichment, x$permutations$enrichment), na.rm = TRUE)
  plot(enrichment, type = "l", col = "red",
       ylim = c(ymin, ymax),
       ylab = "Fold enrichment of GWAS P < 0.05",
       xlab = "Number of genes (effect size cutoff)",
       xaxt = "n", ...)
  n <- length(enrichment)
  spacing <- seq(1, n, by = 50)
  axis(side = 1, at = (1:n)[spacing], padj = 0.5,
       labels = paste0(x$main$sizes, "\n(", round(x$main$intervals, digits = 2),
                       ")")[spacing])
  apply(x$permutations$enrichment, 2, lines, col = "grey75")
  lines(enrichment, col = "red")
  abline(h = 1, col = "blue", lty = 2)
}

#' @export
calc_auc <- function(x) {
  stopifnot(class(x) == "enrichment")
  x$auc$main <- flux::auc(x = 1:length(x$main$enrichment),
                          y = x$main$enrichment)
  x$auc$permutations <- apply(x$permutations$enrichment, 2,
                              function(r) flux::auc(x = 1:length(r), y = r))
  x$auc$signif <- sum(x$auc$main < x$auc$permutations) / length(x$auc$permutations)
  return(x)
}

#' @export
run_gwas_enrich <- function(gene_names,
                            archive = "dec2015.archive.ensembl.org",
                            tss_all_fname = NULL,
                            overwrite = FALSE,
                            window_size,
                            rsid,
                            snp_coords_fname,
                            pval,
                            effect_size,
                            cutoff = 0.05,
                            m = 25,
                            x_direction = "greater",
                            cutoff_direction = "lesser"
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
  snp_genes_final <- add_effect_size(snp_genes_pval_unique, effect_size)
  enrich_result <- enrich_full(x = snp_genes_final$effect_size,
                               y = snp_genes_final$pval,
                               cutoff = cutoff,
                               m = m,
                               x_direction = x_direction,
                               cutoff_direction = cutoff_direction)
  enrich_result <- calc_auc(enrich_result)
  return(enrich_result)
}

boxplot_enrich <- function(auc, permutations) {
  stopifnot(length(auc) == ncol(permutations))
  ymin <- min(auc, permutations, na.rm = TRUE)
  ymax <- max(auc, permutations, na.rm = TRUE)
  boxplot(permutations, ylim = c(ymin, ymax))
  points(auc, col = "red", pch = 19)
}
