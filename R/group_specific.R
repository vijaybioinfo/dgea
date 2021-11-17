#!/usr/bin/R

#######################
# DGEA group specific #
#######################

source("/home/ciro/scripts/handy_functions/devel/utilities.R")
# show_commas, filters_thresholds
source("/home/ciro/scripts/handy_functions/devel/filters.R") # ident_combine
source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
# stats_summary_table moments
source("/home/ciro/scripts/handy_functions/devel/pheatmapCorrection.R")
# pheatmap

#' @title Specific features
#' @description Calculate the specific features for each group.
#' @param results List of the comparisons' data.frames. You need the columns
#' log2FoldChange and padj in each of them.
#' @param vs String used to state it's a comparison, Default: 'vs'.
#' @param sep String to indicate if which groups share features, Default: '&'.
#' @param padjthr Adjusted P-value threshold, Default: 0.05.
#' @param fcthr Fold change threshold, Default: 0.5.
#' @param redundant Add P-value and Fold change to summary.csv, Default: FALSE.
#' @param verbose Show progress, Default: TRUE.
#' @return data.frame with summary of group-specific features.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   gs_df_summary <- group_specific_features(
#'     list(G1vsG2 = comp1_df, G1vsG3 = comp2_df, , G2vsG3 = comp3_df)
#'   )
#' }
#' @seealso
#'  \code{\link[gtools]{mixedsort}}
#' @rdname group_specific_features
#' @export
#' @importFrom gtools mixedsort
group_specific_features <- function(
  results,
  vs = "vs", # symbol used to separate groups in comparisons
  sep = "&",
  padjthr = 0.05,
  fcthr = 0.5,
  redundant = FALSE,
  check = FALSE,
  verbose = TRUE
) {
  # 1. Per DGEA table, give the genes a group based on the thresholds
  # 2. Transform lists into a table
  # 3. Classify
  if(verbose){
    cat("\n-------------------- Group specific --------------------\n")
  }
  if(is.data.frame(results)) results = list(results)
  stat_df = NULL
  if(all(c("avg_logFC", "p_val_adj") %in% colnames(results[[1]]))){
    if(verbose) cat("Seurat results\n")
    if(!is.null(results[[1]]$gene_name)) # safer name
      results[[1]]$gene <- gsub("'", "", as.character(results[[1]]$gene_name))
    results[[1]]$group <- results[[1]]$cluster
    groups = unique(as.character(results[[1]]$group))
    res_de_df <- results[[1]][, c("gene", "group")]
    res_de_df = summarise_table(res_de_df, "gene", sep = sep, verbose = verbose > 1)
    res_de_df$group = as.character(res_de_df$group)
    tvar <- grepl("_mean|_percent", colnames(results[[1]]))
    if(any(tvar)){
      tmp <- !duplicated(as.character(results[[1]]$gene))
      stat_df <- results[[1]][tmp, tvar]
      rownames(stat_df) <- results[[1]]$gene[tmp]
      colnames(stat_df) <- gsub("Seurat.*|CPM.*", "", colnames(stat_df))
    }
    tvar <- c("avg_logFC", "p_val_adj", "gene")
    if(verbose) cat("Transforming to list\n")
    results = lapply(
      X = setNames(groups, paste0(groups, "vsREST")),
      FUN = function(x){
        tmp = results[[1]][which(results[[1]]$group == x), tvar]
        rownames(tmp) <- tmp$gene
        colnames(tmp) <- c("log2FoldChange", "padj", "gene"); tmp
      }
    )
  }else{
    res_de_df = group_specific_classify(
      results = results, vs = vs, sep = sep,
      padjthr = padjthr, fcthr = fcthr, verbose = verbose)
    groups = gtools::mixedsort(unique(unlist(strsplit(names(results), vs))))
  }
  if(verbose) cat("Groups:", show_commas(groups), "\n")
  if(check) return(list(results, res_de_df)) # check _summary step
  res_tests_summary = group_specific_summary(
    results = results,
    features_group = setNames(res_de_df$group, rownames(res_de_df)),
    vs = vs, redundant = redundant, verbose = verbose
  )
  res_de_df_summary <- joindf(res_de_df, res_tests_summary)
  tvar <- names(sort(table(res_de_df_summary$group), decreasing = TRUE))
  tvar <- c(tvar[tvar %in% groups], tvar[!tvar %in% c(groups, "NDE")], "NDE")
  res_de_df_summary$group <- factor(res_de_df_summary$group, tvar)
  res_de_df_summary <- res_de_df_summary[order(res_de_df_summary$group), ]

  if(verbose) cat("\n-------------------- ----- -------- --------------------\n")
  return(list(
    group_specific = res_de_df_summary,
    stats = stat_df
  ))
}

#' @title Classify features
#' @description Use tests' stats to assign the features to each group.
#' @param results List of the comparisons' data.frames. You need the columns
#' log2FoldChange and padj in each of them.
#' @param vs String used to state it's a comparison, Default: 'vs'.
#' @param sep String to indicate if which groups share features, Default: '&'.
#' @param padjthr Adjusted P-value threshold, Default: 0.05.
#' @param fcthr Fold change threshold, Default: 0.5.
#' @param verbose Show progress, Default: TRUE.
#' @return data.frame with classified features indicated in the colum 'group'.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   gs_df <- group_specific_classify(
#'     list(G1vsG2 = comp1_df, G1vsG3 = comp2_df, , G2vsG3 = comp3_df)
#'   )
#' }
#' @seealso
#'  \code{\link[gtools]{mixedsort}}
#'  \code{\link[reshape2]{cast}},\code{\link[reshape2]{melt}}
#'  \code{\link[stringr]{str_count}}
#' @rdname group_specific_classify
#' @export
#' @importFrom gtools mixedsort
#' @importFrom reshape2 dcast melt
#' @importFrom stringr str_count
group_specific_classify <- function(
  results,
  vs = "vs",
  sep = "&",
  padjthr = 0.05,
  fcthr = 0.5,
  verbose = TRUE
) {
  groups = gtools::mixedsort(unique(unlist(strsplit(names(results), vs))))
  if(verbose) cat("Getting differentially expressed\n")
  filters_thresholds_fun = function(x, ...)
    filters_thresholds(x, padjthr, fcthr, ...)
  res_de = lapply(
    X = setNames(names(results), names(results)),
    FUN = function(x){
      setNames(list(
        filters_thresholds_fun(results[[x]], upreg = TRUE),
        filters_thresholds_fun(results[[x]], upreg = FALSE)
      ), unlist(strsplit(x, vs)))
  })
  if(verbose > 1) str(unlist(x = res_de, recursive = FALSE), list.len = 12)
  # Creating data.frame were columns are comparison and rows features
  de_df = reshape2::dcast(
    data = reshape2::melt(res_de), formula = value ~ L1, value.var = "L2"
  )
  rownames(de_df) = as.character(de_df[, 1])
  de_df = de_df[, -1]; de_df[is.na(de_df)] <- "NDE"
  if(verbose > 1) str(de_df)

  if(verbose) cat('Classifying\n')
  de_df$combined = as.character(ident_combine(de_df, names(results), sep = sep))
  de_df$group <- "NDE"
  ntests_names = names(unlist(x = res_de, recursive = FALSE))
  groups_ntests = table(pattern = gsub(".*\\.", replacement = "", ntests_names))
  # must appear the same number of tests and will have at least one NDE
  ii = names(groups_ntests)
  for(i in ii){
    if(verbose) if(which(ii == i) < 10) cat("***", i, "\n") else cat(".")
    j <- stringr::str_count(de_df$combined, i) == groups_ntests[[i]] &
         grepl(pattern = "NDE", x = de_df$combined)
    de_df$group[j] <- i
  }
  group_specs = unique(de_df[de_df$group == "NDE", ]$combined)
  group_specs = sapply(
    X = setNames(strsplit(group_specs, "&"), group_specs),
    FUN = function(x){ unique(x[x %in% groups]) }
  )
  tvar <- sapply(group_specs, length) > 1 & grepl("NDE", names(group_specs))
  ii = names(group_specs[tvar])
  for(i in ii){
    if(verbose) if(which(ii == i) < 10) cat("***", i, "\n") else cat(".")
    nde_tests <- which(strsplit(i, "&")[[1]] == "NDE")
    group_specs_nde <- sapply(
      X = names(de_df)[nde_tests],
      function(x){ all(strsplit(x, vs)[[1]] %in% group_specs[[i]]) }
    )
    if(!all(group_specs_nde)) next
    j <- paste0(group_specs[[i]], collapse = "&")
    de_df$group[de_df$combined == i] <- j
  }; if(verbose) cat("\n")
  tvar <- table(de_df$combined, de_df$group)
  if(verbose > 1 & ncol(tvar) < 10) print(tvar, max = 200)
  de_df <- de_df[, c("group", names(results))]
  de_df <- cbind(gene_name = paste0("'", rownames(de_df)), de_df)
  return(de_df)
}

#' @title Tests summary of group-specific features
#' @description Calculate the summary of tests' stats for each group-specific
#' features.
#' @param results List of the comparisons' data.frames. You need the columns
#' log2FoldChange and padj in each of them.
#' @param features_group Named vector of assigned group with features as a names.
#' @param vs String used to state it's a comparison, Default: 'vs'.
#' @param redundant Add P-value and Fold change to summary.csv, Default: FALSE.
#' @param verbose Show progress, Default: TRUE.
#' @return data.frame with summarised tests' stats across comparisons.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   res_tests_summary = group_specific_summary(
#'     results = list(G1vsG2 = comp1_df, G1vsG3 = comp2_df, , G2vsG3 = comp3_df),
#'     features_group = setNames(gs_df$group, rownames(gs_df))
#'   )
#' }
#' @seealso
#'  \code{\link[data.table]{rbindlist}}
#'  \code{\link[metap]{wilkinsonp}}
#' @rdname group_specific_summary
#' @export
#' @importFrom data.table rbindlist
#' @importFrom metap minimump
group_specific_summary <- function(
  results,
  features_group,
  vs = "vs",
  redundant = FALSE,
  verbose = TRUE
) {
  if(verbose) cat("Merging test statistics\n")
  res_tests = data.frame(row.names = unique(unlist(lapply(results, rownames))))
  for(i in 1:length(results)){
    j <- results[[i]][, c("log2FoldChange", "padj"), drop = FALSE]
    colnames(j) <- paste0(colnames(j), "(", names(results)[i], ")")
    res_tests = joindf(res_tests, j)
  }
  if(isTRUE(redundant)) return(res_tests)
  groups_all = setdiff(unique(features_group), "NDE")
  groups = grep("&", groups_all, value = TRUE, invert = TRUE)
  if(verbose) cat("Summarising test statistics\n")
  tests_summary = data.frame(data.table::rbindlist(lapply(
    X = groups_all, # for each group, and group of groups (eg. G1&G2)
    FUN = function(x){
      if(verbose > 1) cat(x, "\n")
      groups_x = strsplit(x, "&")[[1]]
      genes_x = names(features_group)[features_group == x]
      z = data.table::rbindlist(lapply(
        X = setNames(groups_x, groups_x),
        FUN = function(group_i){
          tvar <- paste0(c("\\(", vs), group_i, c(vs, "\\)"))
          tmp = grepl(paste0(tvar, collapse = "|"), colnames(res_tests))
          y <- res_tests[genes_x, tmp]
          if(all(grepl("REST", colnames(y)))) group_i = "REST"
          tvar <- c(paste0(vs, group_i), paste0(group_i, vs))
          colnames(y) <- gsub(paste0(tvar, collapse = "|"), "", colnames(y))
          return(cbind(
            group_stat = as.character(colnames(y)),
            data.frame(t(y), check.names = FALSE)
          ))
      }));
      z$group_stat = as.character(z$group_stat)
      mydf = data.frame(features = colnames(z)[-1], row.names = colnames(z)[-1],
        check.names = FALSE)
      for(i in unique(z$group_stat)){
        mydf[, i] <- colMeans(z[z$group_stat == i, -1], na.rm = TRUE)
      }
      return(mydf)
    }
  ), fill = TRUE), check.names = FALSE)
  # tests_summary <- tests_summary[!grepl("^NA$|^NA\\.[0-9]", tests_summary$features), ]
  rownames(tests_summary) <- as.character(tests_summary$features)
  tvar <- unlist(lapply(
    X = paste0("\\(", groups, "\\)"),
    FUN = function(x) grep(x, colnames(tests_summary), value = TRUE)
  ))
  tests_summary <- tests_summary[, tvar]
  tests_summary$minimump <- apply(
    X = tests_summary[, grepl("^padj\\(", colnames(tests_summary))],
    MARGIN = 1,
    FUN = function(x) {
      y <- x[!is.na(x)]
      if(length(y) >= 2) metap::minimump(y)$p else y[1]
  }); if(verbose) cat("... done\n")
  return(tests_summary)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param results_list List of the comparisons' data.frames. You need the
#' columns log2FoldChange and padj in each of them.
#' @param path Path to store the results, Default: NULL.
#' @param colours Named vector of colours where the names are the groups,
#' Default: NULL
#' @param return_report Return tables and plots, Default: FALSE.
#' @param optimize Create PNG and not PDF, Default: FALSE.
#' @param verbose Show progress, Default: TRUE.
#' @param mdata Metadata, Default: NULL
#' @param edata Normalised matrix of expression, Default: NULL
#' @param hname Column where groups are in the metadata, Default: NULL
#' @param ... Extra parameters for the 'group_specific_features' function.
#' @return Nothing is returned unles indicated.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   fnames <- list.files(
#'     path = "/path/to/comparions", pattern = "results.csv", full.names = TRUE)
#'   names(fnames) <- basename(dirname(fnames)) # needs to be named
#'   results = lapply(X = fnames, read.csv, row.names = 1, check.names = FALSE)
#'   gs_report <- group_specific_report(
#'     results_list = results,
#'     path = "/path/to/comparions",
#'     mdata = mdata,
#'     edata = edata_norm,
#'     hname = "condition"
#'   )
#' }
#' @seealso
#'  \code{\link[grDevices]{colorRamp}}
#' @rdname group_specific_report
#' @export
#' @importFrom grDevices colorRampPalette
group_specific_report <- function(
  results_list,
  path = NULL,
  colours = NULL,
  return_report = FALSE,
  optimize = FALSE,
  verbose = TRUE,
  # to calculate the stats
  mdata = NULL,
  edata = NULL,
  hname = NULL,
  ...
) {
  if(verbose) cat("\n%%%%%%%%%%%%%%%%%%%% GSpec report %%%%%%%%%%%%%%%%%%%%\n")
  if(!is.null(path) && verbose) cat("Report at:\n", path, "\n")
  the_report = list()
  results_group_specific = group_specific_features(
    results = results_list, ...
  ) # check _summary
  if(isTRUE(list(...)$check)) return(results_group_specific)
  if(verbose > 1) str(results_group_specific)
  summ_df = results_group_specific$group_specific
  if(is.null(results_group_specific$stats)){
    stat_df = data.frame(
      row.names = unique(unlist(lapply(results_list, rownames))))
    stat_df$gene_name = rownames(stat_df)
    tmp <- paste0(levels(summ_df$group),
      rep(moments[c("mn", "p")], each = nlevels(summ_df$group)))
    for(i in tmp){
      for(j in 1:length(results_list)){
        y <- which(grepl(paste0("^", i), colnames(results_list[[j]])))
        if(any(grepl(paste0("^", i), colnames(stat_df)))) next
        if(length(y))
          stat_df = joindf(stat_df, results_list[[j]][, y, drop = FALSE])
      }
    }
    if(!is.null(mdata) & !is.null(edata)){
      results_group_specific$stats = stats_summary_table(
        mat = edata,
        groups = setNames(mdata[, hname], rownames(mdata)),
        rname = rownames(results_group_specific[[1]]),
        moments = c("mn", "p"),
        verbose = verbose
      )
    }else if(ncol(stat_df) > 1){
      if(verbose) cat("------- Stat summary -------\nTaken from results\n")
      if(verbose) warning("Some measurements may be missing from comparisons\n")
      colnames(stat_df) <- gsub(".eurat.*|CPM.*", "", colnames(stat_df))
      results_group_specific$stats = stat_df[, -1, drop = FALSE]; rm(stat_df)
    }
  }

  if(!is.null(results_group_specific$stats)){
    summ_df = joindf(summ_df, results_group_specific$stats)
  }
  # summ_df <- summ_df[order(summ_df$group), ]
  # groups = levels(summ_df$group)[levels(summ_df$group) %in% colnames(summ_df)]
  # tvar <- colnames(summ_df)[!colnames(summ_df) %in% groups]
  # summ_df <- summ_df[, c(tvar, groups)]
  the_report$summary = summ_df
  if(grepl("summary", return_report)) return(the_report)
  if(!isTRUE(return_report))
    write.csv(summ_df, file = paste0(path, "summary.csv"), row.names = FALSE)

  if(!is.null(results_group_specific$stats)){
    if(verbose) cat("--- Heatmap\n")
    annor = summ_df[summ_df[, "group"] != "NDE", "group", drop = FALSE]
    annor[, "group"] <- droplevels(annor[, "group"])
    columns <- grep("_percent", colnames(summ_df), value = TRUE)
    if(length(columns) == 0){ str(summ_df); stop("No stats columns found") }
    tvar <- apply( # keeping features with a % difference > 25
      X = summ_df[rownames(annor), columns],
      MARGIN = 1, FUN = function(x) abs(diff(range(x)))
    )
    labels_row_i <- rownames(annor)
    labels_row_i[!rownames(annor) %in% names(tvar[tvar > 25])] <- ""
    the_report$heatmap = list()
    for(i in c("mn", "p")){
      moments_i = moments[i]
      if(verbose) cat("    *", moments_i, "\n")
      anno_cnames = grepl(paste0(moments_i, "$"), colnames(summ_df))
      mat2plot <- as.matrix(summ_df[rownames(annor), anno_cnames])
      if(moments_i[[1]] == "_mean") mat2plot <- log2(mat2plot + 1)
      if(1){ # row-wise scale and get the range
        mat2plot <- t(scale(t(mat2plot)))
        topz <- max(c(min(abs(c(range(mat2plot), 2))), 1))
      }else{ topz <- 75 } # this was for percentages...
      mat2plot[mat2plot > topz] <- topz; mat2plot[mat2plot < (-topz)] <- -topz;
      colnames(mat2plot) <- sub(moments_i, "", colnames(mat2plot))
      tvar <- intersect(levels(annor[, "group"]), colnames(mat2plot))
      mat2plot <- mat2plot[, tvar]
      annoc = data.frame(
        Group = colnames(mat2plot),
        row.names = colnames(mat2plot))
      colors_i = c(as.list(annor), list(Group = colnames(mat2plot)))
      tmp <- v2cols(unname(unlist(colors_i)), colours)
      colors_i <- lapply(colors_i, v2cols, sour = tmp)
      if(!isTRUE(return_report)){
        fname <- paste0(path, "heatmap", moments_i)
        if(optimize){
          png(paste0(fname, ".png"), width = 1500, height = 1700, res = 250)
        }else{ pdf(paste0(fname, ".pdf"), width = 7, height = 9) }
      }
      cols <- grDevices::colorRampPalette(rev(c('yellow', 'black', 'blue')))(256)
      the_report$heatmap[[i]] = pheatmap(
        mat = mat2plot,
        color = cols, border_color = NA, scale = "none",
        cluster_cols = FALSE, cluster_rows = FALSE,
        annotation_row = annor, annotation_col = annoc,
        annotation_colors = colors_i,
        fontsize_row = 10/(4*log10(nrow(mat2plot))),
        annotation_names_col = FALSE, annotation_names_row = FALSE,
        show_colnames = FALSE, show_rownames = TRUE,
        gaps_col = cumsum(table(annoc$Group)),
        gaps_row = cumsum(table(annor[, "group"])),
        labels_row = labels_row_i
      ); graphics.off()
    }
    if(grepl(paste0("heatmap"), return_report)) return(the_report)
  }
  if(verbose) cat("\n%%%%%%%%%%%%%%%%%%%% %%%%% %%%%%% %%%%%%%%%%%%%%%%%%%%\n")
  if(isTRUE(return_report)) return(the_report)
  return(invisible(x = NULL))
}
