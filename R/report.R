#!/usr/bin/R

###############
# DGEA report #
###############

suppressPackageStartupMessages({
  library(ggplot2)
})

source("/home/ciro/scripts/handy_functions/devel/utilities.R") # getDEGenes
source("/home/ciro/scripts/handy_functions/devel/volcano.R") # volplot
source("/home/ciro/scripts/handy_functions/devel/overlap.R") # overlap_list
source("/home/ciro/scripts/handy_functions/devel/plots.R") # make_breaks
source("/home/ciro/scripts/handy_functions/devel/filters.R") # getDEGenes

dgea_report = function(
  results,
  group1 = "group1",
  group2 = "group2",
  padjthr = 0.05,
  fcthr = 1,
  cols_filt = "minExp",
  output = "./",
  verbose = TRUE,
  return_report = FALSE,
  ##############
  pseuc = 1,
  dtype = "CST",
  # volcano parameters
  showgenes = NULL,
  # heatmap parameters
  mdata = NULL,
  hname = "condition",
  edata = NULL,
  ctrans = "log2",
  couls = NULL
) {
  if(verbose) cat("\n%%%%%%%%%%%%%%%%%%%%%% DGEA report %%%%%%%%%%%%%%%%%%%%%\n")
  means = NULL
  if(is.null(cols_filt)) cols_filt = "pattern123"
  the_report = list()
  if(is.null(names(group1))) names(group1) <- group1
  if(is.null(names(group2))) names(group2) <- group2

  if(verbose) cat("---------------------- Filtering results ---------------\n")
  resdf <- data.frame(
    results[which(results$padj <= 1), ],
    stringsAsFactors = FALSE, check.names = FALSE
  ) # order should ideally be by padj all the way through
  resdf <- resdf[order(resdf$padj), ]
  if(!"gene" %in% colnames(resdf)) resdf$genes = rownames(resdf)
  resdf[, "gene"] <- features_parse_ensembl(resdf[, "gene"])
  rownames(resdf) <- features_parse_ensembl(resdf[, "gene"])
  tvar <- list(
    which(resdf[, "log2FoldChange"] <= -fcthr),
    which(resdf[, "log2FoldChange"] >= fcthr)
  )
  if(!"group" %in% colnames(resdf)){
    if(verbose) cat("Addding 'group' column\n")
    resdf$group = NA; resdf$group[tvar[[1]]] <- group1
    resdf$group[tvar[[2]]] <- group2
  }else{
    tvar <- list(which(resdf$group == group1), which(resdf$group == group2))
  }
  group_m <- sapply(c(names(group1), names(group2)), function(x){
    grep(pattern = paste0("^", x, "_mean"),x = colnames(resdf), value = TRUE)
  })
  if(length(group_m) == 2){
    if(verbose){
      cat("Using means as colour\n"); str(tvar, vec.len = 10, no.list = 4)
      str(group_m, vec.len = 10, no.list = 4)
    }
    resdf$Mean <- 0; means = "Mean"
    resdf$Mean[tvar[[1]]] <- round(log2(resdf[tvar[[1]], group_m[1]] + pseuc), 1)
    resdf$Mean[tvar[[2]]] <- round(log2(resdf[tvar[[2]], group_m[2]] + pseuc), 1)
  }
  genes2plot <- mysignames <- getDEGenes(
    resdf, pv = padjthr, fc = fcthr,
    gene_name = "gene", further = NULL, verbose = verbose
  )
  if(any(grep(cols_filt, colnames(resdf)))){
    tvar <- grep(cols_filt, colnames(resdf), value = TRUE)
    genes2plot <- genes2plot[genes2plot %in% resdf$gene[resdf[, tvar]]]
    if(verbose){ cat("Filtered by", tvar, "\n"); str(genes2plot) }
  }
  resdf$degs <- "Not_significant"
  resdf$degs[resdf[, "gene"] %in% genes2plot] <- "DEG"
  if(length(genes2plot) == 0)
    genes2plot<-getDEGenes(resdf,pv=0.2,fc=fcthr,gene_name="gene",v=TRUE)
  if(length(genes2plot) == 0)
    genes2plot<-resdf[bordering(resdf,cnames="log2FoldChange",n=50),"gene"]
  if(length(group_m) == 2) resdf$Mean[!resdf[, "gene"] %in% genes2plot] <- NA
  resdf$group[!resdf[, "gene"] %in% mysignames] <- NA
  if("pct_diff" %in% colnames(resdf))
    resdf$pct_diff <- abs(resdf$pct_diff)
    resdf[!resdf[, "gene"] %in% genes2plot, ]$pct_diff <- 0
  if(is.null(showgenes)){
    showgenes <- unique(c(
      bordering(resdf[genes2plot, ], cnames = "log2FoldChange", n = 10),
      head(resdf[genes2plot, "gene"], 10)))
  }

  if(verbose) cat("---------------------- Volcano -------------------------\n")
  the_report$volcano <- try(volplot(
    resdf,
    pvalth = padjthr,
    lfcth = fcthr,
    pvaltype = "padj",
    lfctype = "log2FoldChange",
    col_feature = means,
    size_feature = if("pct_diff" %in% colnames(resdf)) "pct_diff",
    gene_name = "gene",
    group = "degs",
    check_genes = list(text = features_parse_ensembl(showgenes)),
    return_plot = TRUE,
    clipp = 4,
    verbose = verbose
  ) + labs(
    size = "Delta %", color = paste0("Mean (", dtype, ")"),
    title = paste(group2, "(-) vs ", group1, "(+)")
  ))
  if(grepl("volcano", return_report)) return(the_report)
  if(!isTRUE(return_report)){
    pdf(paste0(output, 'volcano.pdf'), width = 10, height = 10)
    print(the_report$volcano)
    graphics.off()
  }

  if(verbose) cat("---------------------- Histogram of p-values -----------\n")
  if(!isTRUE(return_report)){
    pdf(paste0(output, "pvalue_distribution.pdf"))
    the_report$histogram = hist(
      resdf$padj, col = seq(0, 1, length.out = 100) <= padjthr,
      breaks = 100, freq = FALSE,
      main = paste("Differentially Expressed Features:", length(mysignames))
    )
    graphics.off()
  }
  if(grepl("histogram", return_report)) return(the_report)

  if(verbose) cat("---------------------- Gene Ontology Analysis ----------\n")
  source('/home/ciro/scripts/gsea/topgo.R')
  tvar <- resdf$gene[1]
  tvar <- ifelse(
    grepl("ENSG", tvar) || tvar == casefold(tvar, upper = TRUE),
    yes = "Human",
    no = "Mouse"
  )
  the_report$go_analysis <- try(go_analysis_table(
    resdf,
    cgroup = "group",
    gene_name = "gene",
    path = if(!isTRUE(return_report)) paste0(output, "GO_Analysis"),
    verbose = verbose,
    organism = tvar
  ))
  if(grepl("go_analysis", return_report)) return(the_report)

  ddsfname <- try(readLines("a1_rdata/ddslocation"), silent = TRUE)
  if(file.exists(ddsfname)){
    load(ddsfname); dds <- dds[genes2plot, ]
  }else{
    dds <- resdf[genes2plot, ]
  }

  if(!is.null(mdata) & !is.null(edata)){
    if(verbose) cat("---------------------- Heatmap -------------------------\n")
    # Further normalisation, currently only log base 2
    if(grepl("log2", ctrans[1])){
      if(verbose)cat("Log2 normalisation\n"); edata <- log2(edata + pseuc);
      dtype <- paste0("log2(", dtype, " + ", pseuc, ")")
    }; if(verbose) cat('Data type for visualisation:', dtype, '\n')

    tmp <- grep(
      ignore.case = TRUE,
      pattern = "donor|type|issue|class|disease",
      value = TRUE, x = colnames(mdata))
    orignames <- filters_columns(mdata, onames = tmp, verbose = verbose)
    orignames <- c(hname, orignames[!orignames %in% hname])
    orignames <- orignames[orignames %in% colnames(mdata)]
    if(verbose) cat("Informative columns:", show_commas(orignames, Inf), "\n")
    symbols_to_plot <- bordering(
      resdf[genes2plot, ], cnames = "log2FoldChange", n = 50)
    mdata[, hname] <- factor(x = mdata[, hname], levels = c(group1, group2))
    data.table::setorderv(mdata, cols = orignames)
    annoc <- mdata[, rev(orignames), drop = FALSE]; ssamples <- rownames(annoc)
    anncolist <- lapply(annoc, function(x){
      v2cols(select = x, sour = couls, v = verbose)
    })

    mat_to_plot <- edata[symbols_to_plot, ssamples]
    tmp <- make.names(names = rownames(mat_to_plot), unique = TRUE)
    rownames(mat_to_plot) <- tmp
    tvar <- rownames(mat_to_plot) %in% mysignames
    rownames(mat_to_plot)[tvar] <- paste0(rownames(mat_to_plot)[tvar], '**')
    titles <- paste('Top', nrow(mat_to_plot), 'genes -',
      sum(symbols_to_plot %in% mysignames), 'significative'
    )

    mat_to_plot <- t(scale(t(mat_to_plot)))
    topz <- max(c(min(abs(c(range(mat_to_plot), 2))), 1))
    mat_to_plot[mat_to_plot > topz] <- topz;
    mat_to_plot[mat_to_plot < (-topz)] <- -topz;
    tvar <- ifelse(nrow(mat_to_plot) > 30, 8, 12)
    fname <- paste0(output, "heatmap.pdf")
    if(!isTRUE(return_report)) pdf(fname, width = tvar, height = 12, onefile = FALSE)
    the_report$heatmap = NMF::aheatmap(
      mat_to_plot, annCol = annoc, annColors = anncolist,
      main = paste(group1, "(-) vs ", group2, "(+)"), sub = titles,
      scale = 'none', Rowv = NA, Colv = NA, labCol = NA,
      col = rev(colorRampPalette(colors = c('yellow', 'black', 'blue'))(20)))
    if(!isTRUE(return_report)) graphics.off()
  }

  if(grepl("3.5|3.4", getRversion()) && !isTRUE(return_report)){
    if(verbose) cat("---------------------- Significative genes report ------\n")
    # deseq2 works with R 3.4, ... maybe 3.5
    source('/home/ciro/scripts/dgea/R/ds2report.R')
    void <- try(ds2report(
      ds2object = dds,
      sgenes = showgenes,
      padjthr = padjthr,
      rewrite = FALSE,
      v = TRUE
    ))
  }

  if(isTRUE(return_report)) return(the_report)
  if(verbose) cat("\n%%%%%%%%%%%%%%%%%%%%%% %%%% %%%%%% %%%%%%%%%%%%%%%%%%%%%\n")
}
