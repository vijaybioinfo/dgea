#!/usr/bin/R
################
# DGEA methos ##
################

# This methods are the ones approved to perform on transcriptomic data
# you can add 'log2cpm' to the method name to perform transformations

##################################### DGEA #####################################
DEAnalysis <- function(
  mat,
  cdt = NULL, # named vector of factors
  mysufix = 'G1vsG2',
  amethod = c('mastlog', 'edgeR', 'edgeR_classic', 'limma', 'scde', 'deseq2', 'DEsingle'),
  latentvar = NULL,
  path = "a1_rdata",
  sepchar = '~'
){
  setwd(path) #; amethod <- match.arg(amethod)
  cat('---------------------------------\n')
  cat('Method:', amethod, '\n')

  if(is.null(cdt)) cdt <- factor(sub(paste0('^(.*)', sepchar, '.*'), '\\1', colnames(mat)) )
  # The input for MAST is a matrix of log2+1 transcripts per million
  mat <- as.matrix(mat)
  if(grepl('cpm', amethod)){ cat('Converting to CPMs\n'); mat <- sweep(mat, 2, colSums(mat), '/') * 1e6}
  if(grepl('log2', amethod)){ cat('Getting shifted log2\n'); mat <- log2(mat + 1) }
  cat("Dimensions:", dim(mat), "\n")
  tvar <- create_formula(latents = latentvar, annotati = annot, amethod = amethod)
  myformula <- tvar[1]; latentvar <- tvar[-1]
  cat('Formula:', myformula, '\n')
  # print(grep("CD4$|CD8B", rownames(mat), value = TRUE))

  if(grepl('mast', amethod)){
    suppressPackageStartupMessages(library(MAST)); packageVersion('MAST')
    # The mast package fits two-part, generalized linear models that are specially adapted for
    # bimodal and/or zero-inflated single cell gene expression data.
    # Check the median of the amount of genes in the cells
    if(ncol(mat) < 30000){
      medgs <- vector()
      for (g in levels(cdt)) {
        medgs[[g]] <- median( colSums(mat[, grepl(paste0('^', g), colnames(mat)), drop = FALSE] > 0) )
      }
      cat('Median amount of features per group:\n'); print(medgs)
    }

    #Calculating the ngeneOn
    cat('Number of features expressed per cell, a.k.a, CDR:\n')
    cdr <- colSums(mat > 0)
    cat('Median:', median(cdr), '\n')
    cat('Mean:', mean(cdr), '\n\n')

    # exprsArray: matrix or array, columns are cells, rows are genes
    # In order to use MAST we need to construct a SingleCellAssay from a matrix or array of expression
    # It basically requires 4 parameters: exprsArray, cData, fData and class
    # cData: cellData an object that can be coerced to a DataFrame. Must have as many rows as ncol(exprsArray)
    # fData: an object that can be coerced to a DataFrame. Must have as many rows as nrow(exprsArray)

    cat('Building object...\n')
    datafMAST <- FromMatrix(
      class = "SingleCellAssay", as.matrix(mat),
      cData = data.frame(
        "wellKey" = colnames(mat),
        "condition" = cdt,
        "ngeneson" = cdr,
        "cngeneson" = scale(cdr)
      ),
      fData = data.frame("primerid" = rownames(mat)))
    rm(mat)

    if(length(latentvar) > 0){
      tvar <- latentvar[latentvar %in% colnames(annot)]
      tvar <- tvar[sapply(tvar, function(x) length(unique(annot[, x])) > 1 )]
      cat('Adding latent variables:', show_commas(tvar), '\n')
      tmp <- sub("^.*~", "", rownames(colData(datafMAST)))
      additional_data <- annot[tmp, tvar, drop = FALSE]
      cat('Scale latent variables\n'); print(head(additional_data))
      additional_data <- apply(X = additional_data, MARGIN = 2, FUN = function(x){
        if(is.numeric(x)) scale(x) else scale(as.numeric(factor(x)))
      }); rownames(additional_data) <- tmp; print(head(additional_data))
      colData(datafMAST) = cbind(colData(datafMAST), additional_data)
      print(sapply(tvar, function(x) head(table(colData(datafMAST)[, x]), 10) ))
    }

    #### Adaptive thresholding ####
    # an adaptive scheme to threshold values below a cut-off that
    # depends on the intensity of the signal cluster
    # from the gene (determined from the median expression value)
    # threshold vs genes binned by median expression value
    cat('Adaptative threshold\n')
    thres <- try( suppressMessages(
      # thresholdSCRNACountMatrix(assay(datafMAST), nbins = 10, min_per_bin = 5)
      thresholdSCRNACountMatrix(assay(datafMAST), nbins = 20, min_per_bin = 30)
    ) ); save(thres, file = paste0('01_thresholdSCRNACountMatrix_', mysufix, '.rdata'))
    if(class(thres)[1] != 'try-error'){
      assays(datafMAST, withDimnames = FALSE) <- list(thresh=thres$counts_threshold, tpm=assay(datafMAST))
    }; threshold_file <- "01_threshold_features.pdf"
    if(!file.exists(threshold_file)){
      thresl <- length(thres$cutpoint)
      if(thresl %in% 1:2) ssize <- 1 else ssize <- round(sqrt(thresl)+.5)
      if(thresl==1){ tvar <- 1 }else{ tvar <- round(thresl/ssize) }
      if(tvar*ssize < thresl) tvar <- tvar + 1
      if(!sqrt(thresl)%%1) tvar <- ssize <- sqrt(thresl)
      if(class(thres)[1] != 'try-error'){
        pdf(threshold_file, width = 5 * tvar, height = 3 * ssize)
        par(mfrow=c(ssize,tvar)); plot(thres); graphics.off()
      }
    }; testpct = 0.2
    expressed_genes <- freq(datafMAST) > testpct # not applied
    # Filter out features < 10 % is suggested by a maintainer
    # https://github.com/RGLab/MAST/issues/135
    cat('If pct > ', testpct, ': ', nrow(datafMAST), ' -> ', sum(expressed_genes),
        ' features to test\n', sep = "")
    cat("Testing", nrow(datafMAST), "features on", ncol(datafMAST), "samples\n")

    ncores <- parallel::detectCores()
    ncores <- ifelse(ncores > 3, 3, ncores)
    cat('Using', ncores, 'cores\n\n')
    options(mc.cores = ncores)

    fcHurdle_file <- paste0('04_fcHurdle_', mysufix, '.rdata')
    if(!file.exists(fcHurdle_file)){
      cat('Regression\n')
      regression_file <- paste0('02_regression_', mysufix, '.rdata')
      if(!file.exists(regression_file)){
        zlm_res <- zlm(as.formula(myformula), datafMAST)
        save(zlm_res, file = regression_file)
      }else{ load(regression_file) }

      condtrast <- paste0('condition', levels(cdt)[2])
      cat('\nSummary - diLRT:', condtrast,'\n')
      # only test the condition coefficient.
      summary_file <- paste0('03_summary_', mysufix, '.rdata')
      if(!file.exists(summary_file)){
        summaryCond <- summary(zlm_res, doLRT = condtrast)
        save(summaryCond, file = summary_file)
      }else{ load(summary_file) }
      cat('Top 4 features by contrast using the logFC\n')
      print(summaryCond, n=4)
      cat('by discrete Z-score\n')
      print(summaryCond, n=4, by='D')
      cat('by continuous Z-score\n')
      print(summaryCond, n=4, by='C'); cat('\n')

      #It looks that our results are in summaryCond
      # see how logFC is calculated at https://rdrr.io/bioc/MAST/man/logFC.html and
      # https://github.com/RGLab/MAST/issues/100 for the base
      cat('Getting results\n')
      summaryDt <- summaryCond$datatable
      tvar <- contrast == condtrast & component == 'H', .(primerid, `Pr(>Chisq)`)
      tmp <- contrast == condtrast & component == 'logFC', .(primerid, coef, ci.hi, ci.lo)
      fcHurdle <- merge(
        z = summaryDt[tvar], # hurdle P values
        y = summaryDt[tmp], # logFC
        by = 'primerid')
      cat('Saving...\n')
      save(fcHurdle, file = fcHurdle_file)
    }else{
      cat('Results provided!\n')
      load(fcHurdle_file)
    }
    fcHurdle[, fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    fcHurdle[, bfr:=p.adjust(`Pr(>Chisq)`, 'bonferroni')]

    # Can't estimate the LFC if one class has no expression for the feature
    # https://support.bioconductor.org/p/99244
    results <- merge(
      x = fcHurdle[fdr <= 1 & abs(coef) > 0], # removes features used to create groups
      y = as.data.table(mcols(datafMAST)),
      by = 'primerid')
    names(results) <- c("gene", "pvalue", "log2FoldChange", "ub", "lb", "padj", "pcor")
    results <- results[, c("gene", "lb", "log2FoldChange", "ub", "pvalue", "padj", "pcor")]
  }
  ##################################### MAST #####################################


  ##################################### edgeR ####################################
  if(sum(grepl('edgeR', amethod))){
    suppressPackageStartupMessages(library(edgeR)); packageVersion('edgeR')
    dge <- DGEList(mat, group = cdt)
    keep <- rowSums(cpm(dge) > 1) >= 2 # filtering
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    tvar <- paste0('calcNormEdgeR_', mysufix, '.rdata')
    if(!file.exists(tvar)){
      cat('calcNormFactors\n')
      dge <- calcNormFactors(dge)
      save(dge, file = tvar)
    }else{
      load(tvar)
    }
    cdr <- scale(colMeans(mat > 0))
    if(sum(amethod == 'edgeR_classic')){
      design <- model.matrix(~ cdt)
    }else{
      design <- model.matrix(~ cdr + cdt)
    }
    tvar <- paste0('estimateDisp_', mysufix, '.rdata')
    if(!file.exists(tvar)){
      cat('estimateDisp\n')
      dge <- estimateDisp(dge, design = design)
      save(dge, file = tvar)
    }else{
      load(tvar)
    }
    tvar <- paste0('glmQLFit_', mysufix, '.rdata')
    if(!file.exists(tvar)){
      cat('glmQLFit\n')
      fit <- glmQLFit(dge, design = design)
      save(fit, file = tvar)
    }else{
      load(tvar)
    }
    tvar <- paste0('glmQLFTest_', mysufix, '.rdata')
    if(!file.exists(tvar)){
      cat('glmQLFTest\n')
      qlf <- glmQLFTest(fit)
      save(qlf, file = tvar)
    }else{
      load(tvar)
    }
    tt <- topTags(qlf, n = Inf)
    results <- tt@.Data[[1]]

    names(results) <- c("log2FoldChange", "logCPM", "F", "pvalue", "padj")
    results$gene <- rownames(results)
    results <- data.table(results[, c("gene", "logCPM", "log2FoldChange", "pvalue", "padj")])
  }
  ##################################### edgeR ####################################

  ##################################### limma ####################################
  if(grepl('limma', amethod)){
    suppressPackageStartupMessages(library(limma)); packageVersion('limma')
    dge <- DGEList(mat, group = cdt)
    tvar <- paste0('calcNormLimma_', mysufix, '.rdata')
    if(!file.exists(tvar)){
      cat('calcNormFactors\n')
      dge <- calcNormFactors(dge)
      save(dge, file = tvar)
    }else{
      load(tvar)
    }
    cdr <- scale(colMeans(mat > 0))
    design <- model.matrix(~ cdt)
    tvar <- paste0('lmFit_', mysufix, '.rdata')
    if(!file.exists(tvar)){
      cat('lmFit\n')
      y <- new("EList")
      y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
      fit <- lmFit(y, design = design)
      save(fit, file = tvar)
    }else{
      load(tvar)
    }
    tvar <- paste0('eBayes_', mysufix, '.rdata')
    if(!file.exists(tvar)){
      cat('eBayes\n')
      fit <- eBayes(fit, trend = TRUE, robust = TRUE)
      save(fit, file = tvar)
    }else{
      load(tvar)
    }

    pdf('mds.pdf', width = 8, height = 8)
    void <- try(limma::plotMDS(dge, col = as.numeric(cdt), pch = 19))
    void <- try(plotMD(fit))
    dev.off()

    results <- topTable(fit, n = Inf, adjust.method = "BH")
    head(results)
    names(results) <- c("log2FoldChange", "AveExpr", "t", "pvalue", "padj", "B")
    results$gene <- rownames(results)
    results <- data.table(results[, c("gene", "AveExpr", "B", "log2FoldChange", "pvalue", "padj")])
  }
  if(grepl('scde', amethod)){
    suppressPackageStartupMessages(library(scde)); packageVersion('scde')
    cd <- clean.counts(mat, min.lib.size = 100, min.reads = 1, min.detected = 3)
    cat('Taking', nrow(cd), 'out of', nrow(mat), 'features\n')
    cat('And', ncol(cd), 'out of', ncol(mat), 'samples\n')
    print(t(headmat(cd, n = 10)))

    #### Fitting error models ####
    # All subsequent calculations will rely on the model
    # The fitting process relies on a subset of robust genes that are detected in
    # multiple cross-cell comparisons.
    # Note: this step takes a considerable amount of time unless multiple cores are used.

    # determine the number of cores
    ncores <- parallel::detectCores()
    ncores <- ifelse(ncores > 3, 3, ncores)
    cat('Using', ncores, 'cores\n\n')
    # determine min nonfailed
    # https://hms-dbmi.github.io/scde/help.html?place=msg%2Fsinglecellstats%2FMPbFH77VIXE%2FM9flYtgyAgAJ
    if(ncol(cd) < 100){
      mnf <- 3
      nk <- ncol(cd) / 4
    }else if(ncol(cd) < 1000){
      mnf <- 5
      nk <- ncol(cd) / 10
    }else{
      mnf <- 10
      nk <- ncol(cd) / 100
    }

    memt <- as.numeric(system("awk '/MemTotal/ {print $2}' /proc/meminfo",intern=TRUE))/1e3
    memf <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo",intern=TRUE))/1e3
    mema <- as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo",intern=TRUE))/1e3
    cat('\nIn system:\n-',memt,'MB of total memory (from /proc/meminfo)\n- ')
    cat(memf,'MB of free memory\n- ')
    cat(mema,'MB of available memory\n')
    cat('R objects used memory: ');print(pryr::mem_used());cat('\n')

    # calculate models
    cat('Calculating error model:\n')
    cd <- apply(cd, 2, function(x){ storage.mode(x) <- 'integer'; x })
    # less than 1000 cells, probably should be knn.error.models directly
    timestamp()
    ifm <- paste0('ifm_', mysufix, '.rdata')
    if(!file.exists(ifm)){
      cat('Using SCDE standard function...\n')
      o.ifm <- try(scde.error.models(counts = cd, groups = cdt, n.cores = ncores, min.nonfailed = mnf,
        save.model.plots = F, min.count.threshold = 4, verbose = 0, save.crossfit.plots = F)
      )
      # Fit single-cell error/regression models
      cat('! done !\n')
      if(class(o.ifm) == 'try-error'){
        cat('\n\nTrying KNN fitting:\n')
        # Build error models for heterogeneous cell populations, based on K- nearest neighbor cells.
        o.ifm <- knn.error.models(counts = cd, groups = cdt, n.cores = ncores, min.nonfailed = mnf,
                      save.model.plots = F, min.count.threshold = 4, verbose = 0, k = nk)
      }; timestamp()
      save(o.ifm, file = ifm)
    }else{ load(ifm) }

    # Output
    # corr.a and corr.b: slope and intercept of the correlated component fit
    # conc.*: the concomitant fit
    # corr.theta: NB over-dispersion
    # fail.r: the background Poisson rate (fixed)

    # Poor cells may result in abnormal fits, most commonly showing negative corr.a, and should be removed.
    valid.cells <- o.ifm$corr.a > 0
    cat('\nNormal fits:', sum(valid.cells), '\n')
    o.ifm <- o.ifm[valid.cells, ]

    # Define an expression magnitude prior for the genes
    # grid of values on which the numerical calculations will be carried out
    pdf('priorDist.pdf', width=10, height = 10)
    o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = T)
    dev.off()

    # Add variables - variables here are defined in parent environment
    if(length(latentvar) > 0){
      cat('Batch correction:', show_commas(latentvar[1]), '\n')
      batchy <- as.factor(c(batch1, latentvar[1]))
    }else{ batchy <- NULL }

    dir.create('scde_GED');
    for(gene in c('IL5', 'IL4', 'IL13', 'IL31', 'CD4', 'GZMB', 'GZMA', 'B2M', 'BCL6', 'CXCR5', 'FOXP3')){
      if(gene %in% rownames(cd)){
        pdf(paste0('scde_GED/test_', gene, '.pdf'))
        scde.test.gene.expression.difference(gene, models = o.ifm, counts = cd,
          prior = o.prior, batch = batchy)
        dev.off()
      }
    }

    #### Testing for differential expression ####
    # first define a factor that specifies which two groups of cells are to be compared
    # factor elements correspond to the rows of the model matrix
    # NA for ells that won't be included in either group
    # define two groups of cells
    cdt <- cdt[row.names(o.ifm)]
    # run differential expression tests on all genes
    cat('Differential expression tests:\n')
    ifm <- paste0('ediff_', mysufix, '.rdata')
    if(!file.exists(ifm)){
      timestamp()
      ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups = cdt, batch = batchy,
                n.randomizations  =  100, n.cores  =  12, return.posteriors = F, verbose  =  1)
      timestamp() ## Output
      save(ediff, file = ifm)
    }else{ load(ifm) }
    # mle: maximum likelihood estimate lb, ub: Bounds of the 95 ce conservative
    # estimate of expression-fold change (equals to the min(abs(c(lb, ub)))
    # 0 if the CI crosses the 0 Z uncorrected Z-score of expression difference

    #### Getting p-values and performing FDR ####
    cat('Getting pvalues and adjusted pvalues from Z-scores\n')
    pvalues <- 2 * pnorm(abs(ediff$Z), lower.tail = F) # 2-tailed p-value
    padj <- 2 * pnorm(abs(ediff$cZ), lower.tail = F) # Adjusted to control for FDR

    results <- cbind(ediff, pvalues, padj)
    colnames(results) <- c('lb', 'log2FoldChange', 'ub', 'ce', 'z', 'cz', 'pvalue', 'padj')
    results$log2FoldChange <- -1 * results$log2FoldChange
    results$gene <- rownames(results); results <- data.table(results)
  }
  if(grepl('deseq2', amethod)){
    cat("Negative Binomial (a.k.a. Gamma-Poisson) distribution\n")
    suppressWarnings(suppressPackageStartupMessages(library(DESeq2)))
    suppressPackageStartupMessages(library(BiocParallel))
    ncores <- system("nproc", intern = TRUE)
    cat("Registering", ncores, "cores\n")
    BiocParallel::register(BiocParallel::MulticoreParam(ncores))
    ddsfname <- readLines("ddslocation")
    cat("Getting deseq object", basename(ddsfname), "\n")
    load(ddsfname)

    cat("-- Contrast and shrinkage --\n")
    cat("Results names:", paste0(resultsNames(dds), collapse = "\n"), "\n")
    # NOTE:  a gene can have a small p-value although the change in expression
    # is not great, as long as the standard error associated with the estimated LFC is small.
    # CONTRAST: the name of the variable, the name of the factor level for the
    # numerator of the log2 ratio, and the name of the factor level for the denominator.
    ctrast <- c("condition", rev(levels(cdt))); cat("Extracting:", show_commas(ctrast), "\n")
    cat("c;"); tvar <- paste0('contrast_', mysufix, '.rdata')
    if(!file.exists(tvar)){
      timestamp(); res <- results(dds, contrast = ctrast, alpha = 0.05)
      save(res, file = tvar)
    }else{ load(tvar) }
    cat("s;"); tvar <- paste0('lfcShrink_', mysufix, '.rdata')
    if(!file.exists(tvar)){
      # ctrast[2:3] %in% levels(colData(dds)[[ctrast[1]]])
      timestamp(); resLFC <- lfcShrink(dds, contrast = ctrast, type = "normal")
      save(resLFC, file = tvar)
    }else{ load(tvar) }
    # final_res <- paste0('final_res_', mysufix, '.rdata')
    pdf("../deseq2_ma_plot.pdf", 9, 9)
    DESeq2::plotMA(res, alpha = 0.05, main = gsub("vs", " vs ", mysufix))
    DESeq2::plotMA(resLFC, alpha = 0.05, main = paste0("Shrunken: ", gsub("vs", " vs ", mysufix)))
    graphics.off(); #rm(res, resLFC)

    # cnames <- c("shrunkenL2FC" = "log2FoldChange", "shrunken_lfcSE" = "lfcSE")
    colnames(res)[-1] <- paste0("unshr_", colnames(res)[-1])
    cnames <- colnames(resLFC)[!colnames(resLFC) %in% colnames(res)]
    names(cnames) <- cnames
    cat("\nAdding columns:", show_commas(cnames, Inf), "\n")
    for(addthis in 1:length(cnames)){
      cat(cnames[addthis], "=>", names(cnames)[addthis], "\n")
      res$adding123 <- resLFC[, cnames[addthis]]
      colnames(res) <- sub("adding123", names(cnames)[addthis], colnames(res))
    }
    summary(res)
    results <- data.frame(res, stringsAsFactors = FALSE)
    results$padj[is.na(results$padj)] <- 2

    results$gene <- rownames(results)
    results <- data.table(results[, c("gene", names(results)[!names(results) %in% "gene"])])
  }
  if(grepl('DEsingle', amethod)){
    # Set the parameters and register the back-end to be used
    ncores <- parallel::detectCores()
    ncores <- ifelse(ncores > 3, 3, ncores)
    param <- SnowParam(workers = ncores, type = "SOCK", progressbar = TRUE)
    BiocParallel::register(param)

    # Detecting the DE genes
    resultf <- paste0('detype_', mysufix, '.rdata')
    if(!file.exists(tvar)){
      tvar <- paste0('contrast_', mysufix, '.rdata')
      if(!file.exists(tvar)){
        results_des <- DEsingle(counts = mat, group = cdt, parallel = TRUE, BPPARAM = param)
        save(results_des, file = tvar)
      }else{ load(tvar) }
      # Dividing the DE genes into 3 categories at threshold of FDR < 0.05
      results.classified <- DEtype(results = results_des, threshold = 0.05)
      save(results_des, file = resultf)
    }else{ load(resultf) }
    results <- results.classified
    results$log2FoldChange <- log2(results$foldChange)
    colnames(results) <- gsub("pvalue.adj.FDR", "padj", colnames(results))
  }

  save(results,  file = paste0('results_', mysufix, '_', amethod, '.rdata'))
  setwd('..')
  cat('---------------------------------\n')
  return(results)
}
##################################### DGEA #####################################

# Add latent variable to formula
create_formula <- function(
  latents, # colum1~variable3~column2
  annotati,
  amethod = NULL
){
  cat("Creating formula\n")
  latents <- as.character(unlist(strsplit(latents, '~')))
  latents_tmp <- latents[latents %in% colnames(annotati)]
  if(length(latents_tmp)){
    cat('Adding latent variables:', show_commas(latents_tmp), '\n')
    latents_tmp <- latents_tmp[sapply(latents_tmp, function(x) length(unique(annotati[, x])) > 1 )]
    cat('Found:', show_commas(latents_tmp), '\n')
  } # design formula
  # 2021-02-03: at some point I removed the cngeneson from the formula!! Re-check your results
  latents_tmp <- latents_tmp#unique(c(latents[grepl("^cngeneson$", latents)], latents_tmp))
  if((!'cngeneson' %in% latents_tmp) && isTRUE(grepl("mast", amethod))){
    latents_tmp <- c("cngeneson", latents_tmp) # add if using mast and it's not included in formula
  }
  latents <- paste(latents_tmp, collapse = " + ")
  variables <- c(latents, "+", "condition")
  variables <- if(isTRUE(grepl("mast", amethod))) rev(variables) else variables
  thisformula <- if(length(latents) > 0 && latents != ""){
    paste(c("~", variables), collapse = " ")
  }else{ "~ condition" }
  # myformula <- paste("~ condition +", latents) # used to use this order, but...
  # Devon Ryan said this one should be used: https://www.biostars.org/p/278684/
  # NOTE: In order to benefit from the default settings of the package,
  # you should put the variable of interest at the end of the formula and make
  # sure the control level is the first level.
  # http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
  # https://sites.tufts.edu/biotools/files/2019/04/bioinformatics_for_rnaseq_day2.pdf
  # I think for MAST, the order doesn't matter because you specify the condition in the doLRT option
  # I may have tested the order and found no difference!
  return(c(thisformula, latents_tmp[!latents_tmp %in% "cngeneson"]))
}

fitting_groups <- function(
  edata,
  annotati,
  latents,
  ddsfname = "dds.Rdata",
  path = "a1_rdata/"
){
  cat("DESeq group fitting for", nrow(annotati), "libraries\n")
  cat("File:", ddsfname, "\n")
  write.table(ddsfname, file = paste0(path, 'ddslocation'), row.names = FALSE, col.names = FALSE, quote = FALSE)
  myformula <- create_formula(latents = latents, annotati = annotati)[1]
  cat("Formula:", myformula, "\n")
  if(file.exists(ddsfname)){ cat("Existing DESeq object.\n"); return(NULL) }

  suppressWarnings(suppressPackageStartupMessages(library(DESeq2)))
  suppressPackageStartupMessages(library(BiocParallel))
  ncores <- system("nproc", intern = TRUE)
  cat("Registering", ncores, "cores\n")
  BiocParallel::register(BiocParallel::MulticoreParam(ncores))
  waiting <- seq(Sys.time(), by = paste(sample(1:20, 1), "sec"), length.out = 2)[2]
  print(waiting - Sys.time())
  while(Sys.time() < waiting){}
  i <- 1; flagfile <- sub(".Rdata", "_running", ddsfname) # ideally wait a random amount of time so just one goes first
  while(file.exists(flagfile)){
    if(i == 1) cat('Waiting for dds to finish'); i <- i + 1
  } # maybe add geneset to the ID too
  if(!file.exists(ddsfname)){
    write.table("Espera...", file = flagfile)
    sttime <- proc.time(); timestamp()
    dds <- DESeqDataSetFromMatrix(countData = edata, colData = annotati, design = as.formula(myformula))
    # The order of the variables of the design do not matter so long as the user
    # specifies the comparison to build a results table for, using the name or
    # contrast arguments of results.
    # dds$condition <- relevel(dds$condition, ref = levels(cdt)[2])
    # dds$condition <- droplevels(dds$condition)
    dds <- DESeq(dds); timestamp()
    # # DESeq runs the following functions in order:
    # dds <- estimateSizeFactors(dds)
    # dds <- estimateDispersions(dds)
    # dds <- nbinomWaldTest(dds)
    print(proc.time() - sttime); cat("Elapsed minutes:", c(proc.time() - sttime)[3] / 60, "\n")
    save(dds, file = ddsfname)
    file.remove(flagfile)
  }
  cat("Finished.\n")
}
