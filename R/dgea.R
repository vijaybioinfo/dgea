#!/usr/bin/R

# ------------------------------------------------------------------------------
# title: Differential gene expression analysis.
# purpose: This script sets the data for a DGE analysis.
# author: Ciro Ramírez-Suástegui
# email: ksuastegui@gmail.com
# date: 2021-02-21
# ------------------------------------------------------------------------------

lib3.5path = "/home/ciro/R/newer_packs_library/3.5"
if(grepl("3.5", getRversion()) && file.exists(lib3.5path))
  .libPaths(new = lib3.5path)

optlist <- list(
  optparse::make_option(
    opt_str = c("-m", "--method"), type = "character", default = "mast",
    help = "Method: for analysis mast, edgeR, limma, scde, or deseq2."
  ),
  optparse::make_option(
    opt_str = c("-a", "--metadata"), type = "character",
    help = "Meta data table: CSV or Rdata (Seurat) object."
  ),
  optparse::make_option(
    opt_str = c("-e", "--expression_data"), type = "character",
    help = paste0("Count table: 10x folder PATH, CSV or Rdata (Seurat).\n\t\t",
                  "If not given, 'metadata' will be taken (Seurat object).\n\t\t",
                  "It may be processed according to later parameters ('ctrans').")
  ),
  optparse::make_option(
    opt_str = c("-u", "--group1"), type = "character",
    help = "Group 1 for comparison: down-regulated genes."
  ),
  optparse::make_option(
    opt_str = c("-d", "--group2"), type = "character",
    help = paste0("Group 2 for comparison.: You can add OTHERS|REST\n\t\t",
                  "to combine the rest of the samples except 'group1.'")
  ),
  optparse::make_option(
    opt_str = c("-n", "--newnames"), type = "character", default = "none",
    help = "Rename: comparison folder name; GROUP1vsGROUP2 (used for plots)."
  ),
  optparse::make_option(
    opt_str = c("-c", "--hname"), type = "character",
    help = "Column name of variable to test the hypothesis on."
  ),
  optparse::make_option(
    opt_str = c("-o", "--output_dir"), type = "character",
    help = "Output directory."
  ),
  optparse::make_option(
    opt_str = c("--sepchar"), type = "character", default = '~',
    help = paste0("Character used to when multiple categories are needed.\n\t\t",
          "In other options' examples '~' is used as the default.")
  ),
  optparse::make_option(
    opt_str = c("-f", "--filters"), type = "character", default = 'mycolumn~myclass',
    help = paste0("Filtering: A single character of mycolumn~myclass or\n\t\t",
          "(quoted) \"list(c('column1','class2'), c('column14','-class1')),\"\n\t\t",
          "or an expression \"column1 == 'group1'\"\n\t\t",
          "names preceeded by '-' indicate exclusion of the group.")
  ),
  optparse::make_option(
    opt_str = c("-b", "--context"), type = "character", default = 'none',
    help = paste0("Sub-name: used for a set of comparisons\n\t\t",
                  "e.g., clusters, disease, etc. Useful for summarization.")
  ),
  optparse::make_option(
    opt_str = c("-l", "--covariates"), type = "character", default = 'cngeneson',
    help = "Covariates: separated by '~', e. g., sex~age~treatment."
  ),
  optparse::make_option(
    opt_str = c("-s", "--down_sample"), type = "character", default = "FALSE",
    help = "Take the same number of samples/cells from 'group1' and 'group2.'"
  ),
  optparse::make_option(
    opt_str = c("-g", "--genesetf"), type = "character", default = 'none',
    help = paste0("Feature filtering:\n\t\t",
           "1. If based on percentage of expressing cells = 'filt_pct'.\n\t\t",
           "2. Based on a regex patter, e. g., '^rps|^rpl|^mp-'.\n\t\t",
           "   Add '~' at the beginning to indicate removal ('~^rps|^rpl').\n\t\t",
           "3. It can also be a file (row names will be taken).\n\t\t",
           "You can also append a further filter at the end of the file\n\t\t",
           "name like in 'FILTERS,' e. g., mycolumn~myclass; or cluster~10.\n\t\t",
           "Note. Indicate filter before DGEA with: 'pre:PATTERN'.")
  ),
  optparse::make_option(
    opt_str = c("-t", "--ctrans"), type = "character", default = "log2",
    help = paste0("Data transformation: for visualisation.\n\t\t",
                  "You can specify separately 'log2' or 'cpm'\n\t\t",
                  "It can be a file, and if further transformations\n\t\t",
                  "you may specify after 'filename~log2'.")
  ),
  optparse::make_option(
    opt_str = c("-p", "--padj_threshold"), type = "numeric", default = 0.05,
    help = "Adjusted P-value threshold."
  ),
  optparse::make_option(
    opt_str = c("-r", "--fc_threshold"), type = "numeric", default = 0.25,
    help = "Fold change: ratio '0.25' usually for single-cell and '1' for bulk."
  ),
  optparse::make_option(
    opt_str = c("--thresh_min"), type = "numeric", default = 0,
    help = paste0("Threshold of expression: for cutoffs. Used to filter\n\t\t",
                  "low expression genes (in more than 25 % of the samples).")
  ),
  optparse::make_option(
    opt_str = c("--annotation"), type = "character",
    default = "/mnt/BioAdHoc/Groups/vd-ay/RNASeq_Workflow/Reference/GRCh37_annotation.csv",
    help = "Feature annotation: file to feature annotation."
  ),
  optparse::make_option(
    opt_str = c("--add_df"), type = "character", default = 'none',
    help = paste0("Extra results columns: path to tables to add to results\n\t\t",
                  "Separate them with '~'; 'path1~path2.'\n\t\t",
                  "Row 1 must match the gene names.")
  ),
  optparse::make_option(
    opt_str = c("--code"), type = "character",
    help = "R code you want to run inside. It runs just after loading data."
  ),
  optparse::make_option(
    opt_str = c("-v", "--verbose"), default = TRUE,
    help = "Verbose: Show progress."
  )
)
# Getting arguments from command line and setting their values to their respective variable names.
opt <- optparse::parse_args(optparse::OptionParser(option_list = optlist))
## Extra parameters ## ---------------------------------------------------------
# NOTE: if you give csv files, check row.names, default:
rnloc = 1 # row names in tables (usially for CSV/TXT files)
pseuc = 1 # [pseudo]count to add when log-transforming values
myseed = 27 # bringing determinism hehe
vs = 'vs' # string to separate comparison names
min.pct = 5 # minimum percentage of expression per group; 5 or 1 from Seurat
min.diff.pct = -Inf # minimum group's percentage difference

# /share/apps/R/3.5/bin/Rscript /mnt/BioHome/ciro/scripts/dgea/R/dgea.R \
#   --method=deseq2 \
#   --metadata=/home/ciro/large/asthma_airways/raw/clarke_GSE111892/metadata_clarke_bulk.csv \
#   --expression_data=/home/ciro/large/asthma_airways/raw/clarke_GSE111892/clarke_GSE111892_counts.rdata \
#   --group1=Stim --group2=Unstim \
#   --hname="Stim" \
#   --output_dir=/home/ciro/large/asthma_pjs/results/deseq2/clarke_bulk_test \
#   --filters="list(c('Stim', 'Stim', 'Unstim'), c('Tissue_Cell_type', 'Lung-TRM'))" \
#   --context=activation_lung \
#   --ctrans=/home/ciro/large/asthma_airways/raw/clarke_GSE111892/clarke_GSE111892_tpm.rdata~log2 \
#   --thresh_min=0
# Rscript /home/ciro/scripts/dgea/R/dgea_jobs.R -y /home/ciro/scripts/dgea/data/config_test.yaml -s TRUE

# Preparing coloured text
if(suppressMessages(require(crayon))){
  cyan = crayon::cyan; redb = crayon::red$bold
}else{ cyan = red = c }

if(opt$verbose){
  cat(cyan(    "-----------------------------------------\n"))
  cat(redb("- Differential gene expression analysis -\n"))
  cat(cyan(    "-----------------------------------------\n"))
}

if(opt$verbose) cat(cyan('\n----------------------- Packages and functions\n'))
packages_funcs = c(
  "/home/ciro/scripts/handy_functions/devel/file_reading.R", # readfile
  "/home/ciro/scripts/handy_functions/devel/filters.R",
  # filters_pattern, filters_complex, filters_subset_df, sample_even
  "/home/ciro/scripts/handy_functions/devel/utilities.R",
  # features_find_symbols, mat_names, count_transformation
  "/home/ciro/scripts/clustering/R/utilities.R", # get_source_data
  "/home/ciro/scripts/handy_functions/R/stats_summary_table.R", # stats_summary_table
  "/home/ciro/scripts/dgea/R/methods.R", # fitting_groups, DEAnalysis
  "data.table"
)
loaded <- lapply(X = packages_funcs, FUN = function(x){
  cat("*", x, "\n")
  if(!file.exists(x)){
    suppressMessages(require(package = x, quietly = TRUE, character.only = TRUE))
  }else{ source(x) }
})
checkclass <- function(grpn, verbose = FALSE){ # to check groups when merging
  if(verbose) cat('For:', grpn,'\n--- ')
  classes <- unlist(base::strsplit(grpn, split = opt$sepchar))
  for(tvar in classes){
    if(!tvar %in% annot[, opt$hname]) stop('Category ', tvar, ' not found, merging is not possible.')
  }
  if(verbose) cat('all is well\n')
}
geneset_fun <- function(x, y, verbose = FALSE, sepchar = "~"){
  yy <- y
  genesetf <- unlist(strsplit(gsub(".*:", x), sepchar))
  # if you give a file make sure the gene names are in the first column
  if(file.exists(genesetf[[1]][1])){
    if(verbose) cat('File of subset of genes\n')
    yy <- readfile(genesetf[[1]][1], stringsAsFactors = FALSE, row.names = 1)
    print(head(yy[, head(colnames(yy)), drop = FALSE]))
    yy <- if(genesetf[2] %in% colnames(yy)){
      filters_subset_df(genesetf[-1], yy, v = T)
    }else{
      if(verbose) cat('Using all genes in file\n'); rownames(yy)
    }; yy <- show_found(yy, y, 'Features', v = TRUE)
  }else if(isTRUE(genesetf[[1]][1] != "none")){
    if(verbose) cat('Based on pattern\n')
    yy <- grep(pattern = gsub("^~", "", x),
      x = y, value = TRUE,
      invert = grepl("^~", x))
  }
  if(verbose)
    if(length(yy) != length(y)) cat("FILTERED\n") else cat("No filtered\n")
  return(yy)
}

if(opt$verbose) cat(cyan('\n-------------------- Digesting parameters\n'))
opt$filters <- filters_pattern(opt$filters)
if(length(opt$filters) == 0) opt$filters <- "none"
str(opt)

if(opt$verbose) cat(cyan('\n------------------- Directories structure\n'))
if(!grepl("scratch|beegfs", getwd())){
  cat("No scratch folder involved; careful about temp files...\n")
  dir.create(opt$output_dir, recursive = TRUE); setwd(opt$output_dir)
}
if(opt$context != "none"){
  dir.create(opt$context, showWarnings = FALSE); setwd(opt$context)
}
compid <- if(grepl(vs, opt$newnames)){
  opt$newnames
}else{ paste0(opt$group1, vs, opt$group2) }
dir.create(compid, showWarnings = FALSE); setwd(compid)
cat(redb("Working in:", getwd(), "\n"))
dir.create('tmp', showWarnings = FALSE)

if(opt$verbose) cat(cyan('\n-------------------------------- Log file\n'))
global_prefix <- paste0(opt$method)
register_log <- !interactive() && !grepl("scratch|beegfs", getwd())
if(register_log){
  log_file <- paste0(getwd(), "/", global_prefix, '.log')
  if(opt$verbose) cat("Check log file:", log_file, "\n")
  if(!file.exists(log_file)) file.create(log_file)
  out_file <- file(log_file, open = 'wt')
  sink(out_file) ; sink(out_file, type = 'message')
}
if(opt$verbose) cat('Date and time:\n') ; st.time <- timestamp();

if(opt$verbose) cat(cyan('\n---------------------------- Loading data\n'))
Sys.time() # Find 10X directory (barcodes, gene names, counts),
# CSV file, TXT, RDS, Seurat object (@counts)
expr_data <- get_source_data(
  xpath = opt$expression_data,
  metadata = opt$metadata[1],
  verbose = opt$verbose
)
Sys.time()
opt$expression_data <- expr_data$source
annot <- remove.factors(expr_data$mdata)
cts <- expr_data$edata; rm(expr_data)
cat('Meta data:\n'); str(annot)
cat('Expression data:\n'); str(cts)

if(!is.null(opt$code)){
  cat(redb("Running code\n"), opt$code, "\n")
  eval(parse(text = opt$code))
}

if(opt$verbose) cat(cyan('\n------------------------------- Filtering\n'))
# Columns derived from features might not be accurate from the
# main (cts variable) table, so you might need to provide these columns in the
# metadata. A solution would be having the option of using the ctrans matrix
filtereddata <- meta_filtering(
  mdata = annot,
  filters = opt$filters,
  cname = opt$hname,
  sepchar = opt$sepchar,
  cts = cts,
  v = TRUE
)
annot <- filtereddata[[1]]; filters <- filters_summary(filtereddata[[2]])

if(opt$verbose) cat(cyan('\n---------------- Setting comparing groups\n'))
# Checking category and groups are in metadata files
if(!opt$hname %in% colnames(annot)){
  if(opt$verbose) cat('Warning: Category not found in meta data.\n')
  if(opt$verbose) cat('You gave  --->', opt$hname, '\n')
  opt$hname <- c(gsub(' ', '', opt$hname), make.names(opt$hname))
  opt$hname <- opt$hname[opt$hname %in% colnames(annot)]
  if(length(opt$hname) == 0) stop()
  if(opt$verbose) cat("... error corrected: '", opt$hname, "'\n", sep = "")
}

if(opt$verbose){
  cat('Groups in ', opt$hname, ':', sep = "")
  print(table(annot[, opt$hname], useNA = 'always'))
  cat("****", opt$hname, "->", "condition\n")
}
annot$condition <- annot[, opt$hname]

restvar <- if(grepl('^others$|^rest$', casefold(opt$group2))){
  if(opt$verbose) cat(opt$group1, 'vs OTHERS/REST comparison\n')
  tvar <- names(table(annot[, opt$hname])) # will ignore NAs
  group1_vec = unlist(strsplit(opt$group1, opt$sepchar))
  opt$group2 <- paste0(tvar[!tvar %in% group1_vec], collapse = opt$sepchar)
  "OTHERS"
}

# For merged comparisons
if(grepl(opt$sepchar, opt$group1) || grepl(opt$sepchar, opt$group2)){
  if(opt$verbose) cat("\nMerging within", opt$hname, "specified\n")
  checkclass(opt$group1, verbose = opt$verbose)
  checkclass(opt$group2, verbose = opt$verbose)
  grp1replace <- sub(paste0(vs, ".*"), "", compid)
  grp2replace <- sub(paste0(".*", vs), "", compid) # will find complete names ^NAME$
  groups_str <- function(x) paste0(paste0("^", strsplit(x, opt$sepchar)[[1]], "$"), collapse = "|")
  annot$condition <- sub(groups_str(opt$group1), grp1replace, annot[, opt$hname])
  annot$condition <- sub(groups_str(opt$group2), grp2replace, annot[, opt$hname])
  opt$group1 <- grp1replace; opt$group2 <- grp2replace # update the tested group names
  if(!is.null(restvar)){
    annot$condition <- sub(paste0("^", grp2replace, "$"), restvar, annot$condition)
    opt$group2 <- restvar
  }
}

grp1exist <- which(annot == opt$group1, arr.ind = TRUE)
grp2exist <- which(annot == opt$group2, arr.ind = TRUE)
if(any(c(nrow(grp1exist), nrow(grp2exist)) == 0)){ # if either is not found
  if(nrow(grp1exist) == 0) stop(opt$group1 , " not found.")
  if(nrow(grp2exist) == 0) stop(opt$group2 , " not found.")
}else if(opt$verbose) cat("Categories in place!\n")

cts <- cts[, show_found(rownames(annot), colnames(cts), v = opt$verbose)]
features <- rownames(cts)[which(Matrix::rowSums(cts > 0) > 0)];
cat(length(features), "/", nrow(cts), "features with nonzero expression across all samples\n")
cts = cts[rownames(cts) %in% features, ]
cat('-', ncol(cts), 'samples\n-', nrow(cts), 'features\n'); gc()

################################ Group fitting #################################
# In case we want to capture the variance of a certain group. This is mainly for deseq2.
if(grepl("deseq", opt$method)){
  ddsname <- paste0(getwd(), "/../.", opt$method, "_dds_", opt$hname,
    "_FILT", filters, "_", nrow(annot), "samples.Rdata")
  cat("\n"); void <- fitting_groups(
    edata = cts,
    annotati = annot,
    latents = opt$covariates,
    method = opt$method,
    ddsfname = ddsname,
    path = "tmp/"
  )
}
################################ ##### ####### #################################

#### Choosing two 2 population for comparison ####------------------------------
if(opt$verbose) cat(cyan("\n---------------- Building 2-groups matrix\n"))
grp1vector <- filters_subset_df(c("condition", opt$group1), annot, v = TRUE)
grp2vector <- filters_subset_df(c("condition", opt$group2), annot, v = TRUE)
minc <- min(length(grp2vector), length(grp1vector))
# to avoid skewed differential expression when groups are uneven
if(!opt$down_sample[1] %in% c("FALSE", "no", "n", "F")){
  if(opt$verbose) cat("Sampling cells <------------\n")
  if(all(file.exists(c("tmp/g1.csv", "tmp/g2.csv")))){
    grp1vector <- as.character(read.csv("tmp/g1.csv", row.names = 1)[, 1])
    grp2vector <- as.character(read.csv("tmp/g2.csv", row.names = 1)[, 1])
  }else if(any(opt$down_sample %in% c("TRUE", "yes", "y", "T"))){
    for(down_sample_i in opt$down_sample[opt$down_sample %in% colnames(annot)]){
      if(opt$verbose) cat("Even sample", down_sample_i, "\n")
      grp2vector <- sample_even(annot[grp2vector, ],
        down_sample_i, verbose = opt$verbose)
      grp1vector <- sample_even(annot[grp1vector, ],
        down_sample_i, verbose = opt$verbose)
    }
  }; minc <- min(length(grp2vector), length(grp1vector))
  minc <- ifelse(minc > 8000, 7500, minc)
  if(any(opt$down_sample %in% c("TRUE", "yes", "y", "T"))){
    set.seed(myseed); grp2vector <- sample(grp2vector, minc)
    set.seed(myseed); grp1vector <- sample(grp1vector, minc)
  }
}; minc <- min(length(grp2vector), length(grp1vector))
if(!file.exists("tmp/g1.csv")) write.csv(grp1vector, file = "tmp/g1.csv")
if(!file.exists("tmp/g2.csv")) write.csv(grp2vector, file = "tmp/g2.csv")

# Now finally filtering the matrix
if(opt$verbose) cat("Final number of samples\n")
tvar <- show_found(c(grp1vector, grp2vector), colnames(cts), v = opt$verbose)
cts <- as.matrix(cts[, tvar])
ecomp <- cbind(cts[, grp1vector, drop = FALSE], cts[, grp2vector, drop = FALSE])
rm(cts)
if(!exists('datavis')) datavis <- ecomp
# we shouldn't duplicate the matrix because it takes a lot of memory
# but we need to generate a stat report and this
# will help us filter genes
annot <- annot[colnames(ecomp), ]
# Avoid eliminating the gene by which you are creating the comparing groups
# e.g. FOXP3p vs FOXP3n cells/samples; at least a sample will express it in both groups
ecomp[, grp1vector[1]] <- ecomp[, grp1vector[1]] + 1
ecomp[, grp2vector[1]] <- ecomp[, grp2vector[1]] + 1
# add the group as a prefix to the sample/cell names
colnames(ecomp) <- c(
  paste(opt$group1, grp1vector, sep = opt$sepchar),
  paste(opt$group2, grp2vector, sep = opt$sepchar)
); str(ecomp)

# Determine matrix for visualisation
ctrans <- unlist(strsplit(opt$ctrans, opt$sepchar))
if(file.exists(ctrans[1])){
  opt$expression_data <- ctrans[1]
  if(opt$verbose) cat("Getting matrix for visualisations\n")
  datavis <- readfile(ctrans[1], row.names = rnloc, check.names = FALSE, verbose = TRUE)
  if(class(datavis) != "matrix") datavis <- as.matrix(datavis)
  ctrans <- if(length(ctrans[-1]) > 0) rev(ctrans) else "" # keep file name for norm type
}

datavis <- datavis[, c(grp1vector, grp2vector)]

# First we want to normalise data
dtype <- "CTS"; # currently only CPM transformation
if(ctrans[1] != 'FALSE'){
  dtype <- casefold(sub("log2", "", ctrans[1]), upper = TRUE) # extract transformation type
  if(grepl("cpm", casefold(ctrans[1])))
    datavis <- count_transformation(datavis, transf = "cpm", verbose = opt$verbose)[[1]]
  if(dtype == ""){ # if no info, try to identify the type of data
    csums <- which(c(1e6, 1e4) %in% mean(colSums(datavis))) # per million or 10,000 (used by Seurat)
    if(length(csums) == 0){ csums = "sad"; dtype = "CTS" }
    tvar <- any(grepl("tpm", casefold(c(opt$expression_data, basename(ctrans))))) # check in the name too
    dtype <- ifelse(tvar, "TPM", switch(csums, "1" = "CPM", "2" = "CP10K", dtype))
    if(opt$verbose) cat('Count autodetection:', dtype, '\n')
  }
}

if(opt$verbose) cat(cyan('\n------------------ Calculating statistics\n'))
dstats <- stats_summary_table(
  mat = datavis,
  groups = setNames(annot$condition, rownames(annot)),
  moments = c('b', 'mn', 'p'),
  expr_cutoff = opt$thresh_min, # for positive cells/sample percentage calculation
  datatype = dtype,
  verbose = opt$verbose
)
dstats$pct_diff <- apply(
  X = dstats[, paste0(c(opt$group1, opt$group2), "_percentage", dtype)],
  MARGIN = 1,
  FUN = diff
)
features <- rownames(dstats)[which(dstats$BexprFrac > 0)];
if(opt$verbose){
  cat(
    length(features), "/", nrow(datavis), "features expressed in at least one samples",
    "with minimum threshold of", opt$thresh_min, "\n"
  )
  cat("More filters: minimum % in a group >", min.pct, "and min diff % >", min.diff.pct, "\n")
}
tvar <- dstats[features, ]
tvar$min <- matrixStats::rowMins(as.matrix(tvar[, paste0(c(opt$group1, opt$group2), "_percentage", dtype)]))
features2 <- which(x = tvar$min > min.pct)
if(1){
  if(opt$verbose) cat("1.", length(x = features2), "features pass", min.pct, "threshold\n")
  features2 <- which(x = tvar$min > min.pct & abs(tvar$pct_diff) > min.diff.pct)
  if(opt$verbose) cat("2.", length(x = features2), "features pass", min.diff.pct, "threshold\n")
  if(opt$genesetf[[1]] == 'filt_pct' && length(x = features2) > 0){
    features <- rownames(dstats)[features2]
  }else if(opt$verbose) cat("NOTE: filters 1 and 2 were not applied.\n")
  if(opt$verbose) cat("Testing", length(features), "/", nrow(ecomp), "features\n")
}
if(isTRUE(grepl("^pre", opt$genesetf))){
  features <- geneset_fun(
    sub("^pre:", "", opt$genesetf), features, opt$verbose, sepchar = opt$sepchar
  )
}
# Filter the genes from the matrices and stat report
ecomp <- ecomp[rownames(ecomp) %in% features, ]
datavis <- datavis[rownames(datavis) %in% features, ]
dstats <- dstats[rownames(dstats) %in% features, ]

##################################### DGEA #####################################
if(opt$verbose) cat('C1:', opt$group1, '- C2 (will have positive fold changes):', opt$group2, '\n')
cond <- factor(sub(paste0('^(', opt$group1, '|', opt$group2, ').*'), '\\1', colnames(ecomp)) )
names(cond) <- colnames(ecomp)
if(opt$verbose) cat('Relevel:', opt$group2, '\n'); cond <- relevel(cond, opt$group2)
if(opt$verbose){ cat('Comparison proportion:'); print(table(cond)) }
res <- DEAnalysis(
  mat = ecomp,
  cdt = cond,
  mysufix = compid,
  amethod = opt$method,
  latentvar = opt$covariates,
  path = "tmp"
)
##################################### #### #####################################

if(opt$verbose) cat(red('\n------------------------- Finishing table\n'))
setorder(res, -log2FoldChange)
if(opt$verbose) cat('Table of thresholds:\n')
threstable <- data.frame(mat_names(
  c(0, 0.25, 0.6, 1, 2), # fold changes
  c(0.2, 0.05, 0.01, 0.001, 0.0001, 0.00001) # p-values
), check.names = 0)
for(p in colnames(threstable)){
  for(fc in rownames(threstable)){
    threstable[fc, p] <- length(getDEGenes(res, pv = as.numeric(p), fc = as.numeric(fc), gene_name = "gene"))
  }
}
# threstable[is.na(threstable)] <- rnorm(n = prod(dim(threstable)), mean = 100, sd = 20)
ddf <- reshape2::melt(
  cbind(log2FoldChange = rownames(threstable), threstable),
  variable.name = 'padj', value.name = 'nDEG'
)
threstable <- rbind(
  threstable,
  mat_names("", colnames(threstable)), round(threstable / nrow(res) * 100, 2)
)
if(opt$verbose) print(threstable)
write.csv(threstable, file = paste0(global_prefix, '_nDEGs.csv'), quote = FALSE)

if(opt$verbose) cat('Add groups, hyphened genes, and stats\n')
genesitos <- as.character(res[, gene])
temp <- genesitos %in% res[which(res[, log2FoldChange] > 0), gene]
temp <- setNames(ifelse(temp, opt$group1, opt$group2), genesitos)
degs <- filters_thresholds(
  x = data.frame(res, stringsAsFactors = FALSE),
  pv = opt$padj_threshold, fc = opt$fc_threshold, gene_name = "gene"
); temp[!names(temp) %in% degs] <- ""
res <- cbind(
  res, excel_name = paste0("'", genesitos),
  group = temp,
  dstats[res[, gene], ]
)
# Feature must be expressed (thresh_min) in at least 25 % of the samples (min_samp)
min_samp <- ifelse(min(table(annot[, opt$hname])) > 3, round(ncol(datavis) / 4), 3)
tvar <- rowSums(datavis > opt$thresh_min, na.rm = TRUE) > min_samp
if(sum(tvar) > 1){
  if(opt$verbose){
    cat("****", sum(tvar), "features with expression >",
    opt$thresh_min, "in >", min_samp, "samples\n")
  }
  res$minExp <- as.character(res[, gene]) %in% rownames(datavis)[tvar]
  tvar <- paste0("minExp", opt$thresh_min, "in", min_samp, "samples")
  colnames(res) <- gsub("minExp", tvar, colnames(res))
}
if(file.exists(opt$annotation)){
  tvar <- features_find_symbols(
    features_init = genesitos,
    annotation = opt$annotation
  )
  res <- if(!is.null(tvar$annotation)){
    joindf(tvar$annotation, data.frame(res, stringsAsFactors = FALSE))
  }else{ res }
}
if(opt$verbose) str(res)

tmp <- res$gene %in% geneset_fun(opt$genesetf, res$gene, opt$verbose)
if(!all(tmp)) res$filtered <- tmp

# excessive info to add — just read this and have no idea what I meant
fname <- paste0(global_prefix, '_results.csv')
if(opt$verbose) cat('Writing results:', fname, '\n')
write.csv(res, file = fname, quote = FALSE, row.names = FALSE)

if(!is.null(res$filtered)) res <- res[res$filtered,]

if(opt$verbose) cat(cyan('\n------ [Symbolically] linking input files\n'))
tvar <- tools::file_ext(opt$expression_data)
if(tvar != "") tvar <- paste0(".", tvar)
fnames <- list(
  c(opt$expression_data, paste0('expression_data', tvar)),
  c(opt$metadata, paste0('metadata.', tools::file_ext(opt$metadata)))
); if(opt$verbose) cat("Adding source files:\n")
void <- sapply(fnames, function(x){
  command <- paste(c('ln -s', x), collapse = " ")
  cat(command); if(!file.exists(x[2])) system(command); cat("\n")
})

source("/home/ciro/scripts/dgea/R/report.R")
try(
  dgea_report(
    results = res,
    group1 = opt$group1,
    group2 = opt$group2,
    padjthr = opt$padj_threshold,
    fcthr = opt$fc_threshold,
    output = paste0(global_prefix, "_"),
    verbose = opt$verbose,
    ##############
    pseuc = pseuc,
    dtype = dtype,
    # heatmap parameters
    mdata = annot,
    hname = opt$hname,
    edata = datavis,
    ctrans = ctrans,
    couls = NULL
  )
)

if(opt$verbose){
  cat('\n\n*******************************************************************\n')
  cat('Starting time:\n'); cat(st.time, '\n')
  cat('Finishing time:\n'); timestamp()
  cat('*******************************************************************\n')
  cat('SESSION INFO:\n'); print(sessionInfo()); cat("\n")
  cat('Pipeline finished successfully\n')
}
if(register_log){
  sink(type = 'message')
  sink()
}
