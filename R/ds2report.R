#!/usr/bin/R

#' DGEA HTML Report
#'
#' This function performs semi-automated differential gene expression analysis
#' HTML Report
#'
#' @param ds2object DESeq2 object or data.frame with results.
#' @param padjthr Adjusted P-value threshold.
#' @param reportdir Report directory name.
#' @param rewrite Re-write results, otherwise append.
#' @param verbose Verbose.
#' @keywords DGEA, OrgDb
#' @references \url{https://rdrr.io/bioc/ReportingTools/man/HTMLReport.html}
#'
#' @return NULL
#'
#' @importFrom org.Hs.eg.db org.Mm.eg.db ReportingTools DESeq2
#'
#' @export
#'
#' @examples
# ' void <- ds2report(ds2object = myobject)
#'

ds2report <- function(
  ds2object,
  padjthr = 0.05,
  sgenes = NULL,
  methd = "DESeq2",
  reportdir = "report",
  rewrite = FALSE,
  verbose = FALSE
){
  if(verbose){ cat("--- DGEA HTML Report ---\n"); timestamp() }
  suppressPackageStartupMessages(library(ReportingTools))
  suppressPackageStartupMessages(library(DESeq2))
  if(isTRUE(rewrite) && dir.exists(reportdir)) system(paste("rm -r", reportdir))
  keyt <- if(all(grepl("ENS", rownames(ds2object)[1:10]))) "ENSEMBL" else "SYMBOL"
  spec <- all(grepl("ENSG", rownames(ds2object)[1:10]))
  if(!spec) spec <- !all(grepl("ENSMUS", rownames(ds2object)[1:10]))

  if(verbose){
    cat("Using", keyt, "to find Entrez IDs (integers)\n")
    cat("Data base:", ifelse(spec, "Homo sapiens", "Mus musculus"), "\n")
  }

  sgenes <- sgenes[sgenes %in% rownames(ds2object)]
  if(!is.null(sgenes)){
    if(verbose) cat("Selected features:", length(sgenes), "\n")
    ds2object <- ds2object[sgenes, ]
  }
  ds2object <- ds2object[rowMeans(counts(ds2object)) > 1, ]
  rownames(ds2object) <- features_parse_ensembl(rownames(ds2object), keepens = TRUE)
  tvar <- gsub("\\..*", "", rownames(ds2object))
  if(any(duplicated(tvar))){
    if(verbose) cat("Duplicated gene names\n")
    tvar <- rownames(ds2object)
  }
  rownames(ds2object) <- tvar

  spec <- ifelse(spec, "Hs", "Mm")
  eval(parse(text = paste0("library(org.", spec, ".eg.db)")))
  eval(parse(text = paste0("uniKeys <- keys(x = org.", spec, ".eg.db, keytype = keyt)")))
  eval(parse(text = paste0(
    "gannot <- select(x = org.", spec, ".eg.db, keys = uniKeys,",
    "columns = c('SYMBOL', 'PATH', 'ENTREZID', 'ENSEMBL'), keytype = keyt)"
  )))

  gannot$gene_names <- gannot[, keyt]
  gannot <- gannot[!is.na(gannot$gene_names), ]
  gannot <- gannot[!duplicated(gannot$gene_names), ]
  rownames(gannot) <- gannot$gene_names
  if(verbose) str(gannot)
  tvar <- rownames(ds2object) %in% gannot$gene_names
  if(verbose) cat("Found features:", sum(tvar), "\n")
  if(sum(tvar) == 0){ return(NULL)}
  ds2object <- ds2object[tvar, ]
  rownames(ds2object) <- gannot[rownames(ds2object), 'ENTREZID']

  if(verbose) cat("Starting report\n")
  des2Report <- HTMLReport(
    shortName = paste0(casefold(methd), '_report'),
    title = paste0('RNA-seq analysis of differential gene expression using ', methd),
    reportDirectory = reportdir
  )

  if(casefold(class(ds2object)) == 'deseqdataset'){
    factores <- colData(ds2object)$condition
  }else{ factores = NULL }

  if(verbose) cat("Publishing\n")
  publish(
    object = ds2object,
    publicationType = des2Report,
    pvalueCutoff = padjthr,
    annotation.db = paste0("org.", spec, ".eg.db"),
    factor = factores,
    reportDir = reportdir
  ); if(verbose){ cat("----------\n"); ; timestamp()}

  finish(des2Report)
  return(NULL)
}
