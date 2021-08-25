#!/bin/R

library(optparse)

optlist <- list(
  make_option(
    opt_str = c("-y", "--yaml"), type = "character", default = "config_project.yaml",
    help = "Configuration file."
  ),
  make_option(
    opt_str = c("-s", "--submit"), type = "logical", default = FALSE,
    help = "Submit jobs."
  ),
  make_option(
    opt_str = c("-v", "--verbose"), type = "logical", default = TRUE,
    help = "Verbose."
  )
)
optparse <- OptionParser(option_list = optlist)
opt <- parse_args(optparse)
if(interactive()) opt$yaml = "/home/ciro/scripts/dgea/config.yaml"

## Functions ## ----------------------------------------------------------------
dirnamen <- function(x, n = 1){
  for(i in 1:n) x <- dirname(x)
  return(x)
}
running_jobs <- function(){
  system("qstat -fu ${USER} | grep -E 'Job_Name|Job Id|job_state' | sed 's/Id: /Id_/g; s/ = /: /g; s/.herman.*/:/g' > ~/.tmp")
  jobs_yaml = yaml::read_yaml("~/.tmp")
  jobs_yaml <- jobs_yaml[sapply(jobs_yaml, function(x) x[['job_state']] ) != "C"]
  jobs_df <- data.frame(
    id = gsub(".*Id_", "", names(jobs_yaml)),
    Name = sapply(jobs_yaml, function(x) x[["Job_Name"]] ),
    stringsAsFactors = FALSE
  ); rownames(jobs_df) <- NULL
  jobs_df
}

## Reading files ## ------------------------------------------------------------
config_file = yaml::read_yaml(opt$yaml)

dir.create(config_file$output_dir, showWarnings = FALSE)
setwd(config_file$output_dir)
dir.create(config_file$project, showWarnings = FALSE)
setwd(config_file$project)
dir.create("scripts", showWarnings = FALSE)
if(opt$verbose){
  cat(crayon::cyan("-------------------------------------\n"))
  cat(crayon::red$bold("------------ DGE Analysis\n"))
  cat("Configuration:")
  str(config_file[!names(config_file) %in% "job"])
  cat("Working at:", getwd(), "\n")
  system("ls -loh")
  cat("\n")
}

if(opt$verbose) cat("Getting job template\n")
template_pbs_con <- file(description = config_file$job$template, open = "r")
template_pbs <- readLines(con = template_pbs_con)
close(template_pbs_con)

# Setting lements to iterate through
if(opt$verbose) cat("Fetching comparisons\n")
if(is.null(config_file$comparisons$file)) config_file$comparisons$file = "none"
comparisons = config_file$comparisons[!names(config_file$comparisons) %in% "file"]
if(file.exists(config_file$comparisons$file)){
  if(opt$verbose) cat("- from file\n")
  comps_tab = read.csv(config_file$comparisons$file, stringsAsFactors = FALSE)
  more_comparions <- lapply(
    X = 1:nrow(comps_tab),
    FUN = function(x){
      c(list(contrast = unname(unlist(comps_tab[x, 1:2]))),
        as.list(comps_tab[x, -c(1:2)]))
  })
  comparisons <- c(comparisons, more_comparions)
}
# expanding those with more than two group to compare
tvar <- unname(which(x = sapply(comparisons, function(x) length(x$contrast) ) > 2))
for(i in tvar){
  more_comparions = apply(
    X = combn(x = comparisons[[i]]$contrast, m = 2),
    MARGIN = 2, FUN = function(x){
      y <- comparisons[[i]]
      y$contrast <- x
      return(y)
  }); comparisons <- c(comparisons, more_comparions)
}; if(length(tvar) > 0) comparisons <- comparisons[-tvar]
# Getting contexts
# tvar <- names(comparisons) == ""
# if(sum(tvar)) names(comparisons)[tvar] <- paste0("comparisons", 1:sum(tvar))
# for(i in 1:length(comparisons)){
#   tvar <- comparisons[[i]]$context
#   comparisons[[i]]$context <- ifelse(is.null(tvar), names(comparisons[i]), tvar)
# }

if(opt$verbose) cat("Processing", length(comparisons), "comparisons\n")
if(opt$verbose) cat("--------------------------------------\n")
for(m in config_file$method){
  for(config_comp in comparisons){
    if(opt$verbose) cat(config_comp$contrast, "in", config_comp$test_column)
    my_comparison <- if(isTRUE(grepl("vs", config_comp$name))){
      config_comp$name
    }else{ paste0(config_comp$contrast, collapse = "vs") }
    res_file0 <- paste0(config_comp$context, "/", my_comparison, "/", m, "_results.csv")
    res_file <- paste0(getwd(), "/", res_file0)
    re_run <- isTRUE(any(grepl("^f$|force", c(config_file$job$submit, opt$submit))))
    if(file.exists(res_file) && !re_run){
      if(opt$verbose) cat(crayon::green$bold(" - done\n")); next
    }
    my_filter = if(!is.null(names(config_comp$filter))){
      unlist(lapply(
        X = names(config_comp$filter),
        FUN = function(x){
          if(all(!grepl("^expr", config_comp$filter[[x]][1]))){
            paste0("c('", x, "', '", paste0(config_comp$filter[[x]], collapse = "', '"), "')")
          }else{ config_comp$filter[[x]] }
        }
      ))
    }else{ if(isTRUE(config_comp$filter != "none")) config_comp$filter }
    tvar <- all(!grepl("^expr", my_filter)) && !any(grepl("list|~", my_filter))
    if(tvar) my_filter <- paste0("list(", paste0(my_filter, collapse = ", "), ")")
    # tvar <- gsub("(^..).*(..$)", "\\1\\2", digest::digest(my_filter, "md5", serialize = FALSE))
    routine_pbs_file = paste0(c(
      gsub("(^..).*", "\\1", m), config_comp$context,
      config_comp$test_column, my_comparison#, tvar
    ), collapse = "_")

    # Parameters
    params <- paste0(
      config_file$exec,
      " ", config_file$script,
      " --method=", m,
      " --metadata=", config_file$metadata,
      " --expression_data=", config_file$expression_data,
      " --group1=", config_comp$contrast[1],
      " --group2=", config_comp$contrast[2],
      " --newnames=", my_comparison,
      " --hname=", config_comp$test_column,
      " --output_dir=", getwd(),
      " --filters=\"", my_filter, "\"",
      " --context=", ifelse(is.null(config_comp$context), "none", config_comp$context),
      " --covariates=", config_file$covariates,
      " --down_sample=", config_file$down_sample,
      " --ctrans=", config_file$ctrans,
      " --padj_threshold=", config_file$padj_threshold,
      " --fc_threshold=", config_file$fc_threshold,
      " --annotation=", config_file$annotation,
      " --thresh_min=", 0,
      " --verbose=", TRUE,
      " "
    )

    # Output directory
    output_dir <- getwd()

    pbs <- gsub("cellranger", "dgea", template_pbs)
    pbs <- gsub("\\{username\\}", Sys.info()[["user"]], pbs)
    pbs <- gsub("\\{sampleid\\}", my_comparison, pbs)
    pbs <- gsub("\\{routine_pbs\\}", routine_pbs_file, pbs)
    pbs <- gsub("\\{outpath\\}", output_dir, pbs)
    pbs <- gsub("\\{routine_params\\}", params, pbs)
    if(file.exists(res_file)){
      pbs <- gsub("# \\{after_copy\\}", paste("mkdir --parents", dirname(res_file0)), pbs)
      tvar <- paste0("cp -R ", dirname(res_file), "/. ./", dirname(res_file0), "/")
      pbs <- gsub("# \\{pre_routine\\}", tvar, pbs)
    }
    tvar <- if(is.null(config_comp$job)) config_file$job$main else config_comp$job
    for(i in names(config_file$job$main)){
      job_parm <- if(i %in% names(tvar)) tvar[[i]] else config_file$job$main[[i]]
      job_parm <- if(m %in% names(job_parm)) job_parm[[m]] else job_parm[[1]]
      pbs <- gsub(paste0("\\{", i, "\\}"), job_parm, pbs)
    }
    pbs <- gsub("\\.\\.\\/scripts", "scripts", pbs)

    running <- try(running_jobs(), silent = TRUE)
    if(class(running) == "try-error") running <- list(Name = "none", id = "X124")
    do_we_submit <- any(c(opt$submit, re_run))
    if(any(grepl(routine_pbs_file, running$Name)) && isFALSE(do_we_submit)){
      if(opt$verbose) cat(" - running\n"); next
    }

    pbs_file <- paste0(getwd(), "/scripts/", routine_pbs_file, ".sh")
    writeLines(text = pbs, con = pbs_file)
    if(isTRUE(any(c(config_file$job$submit, do_we_submit)))){
      depend <- if(isTRUE(config_file$job$depend %in% running$id)) paste0("-W depend=afterok:", config_file$job$depend)
      pbs_command <- paste("qsub", depend, pbs_file)
      if(opt$verbose) cat("\n", pbs_command, "\n"); system(pbs_command)
      void <- suppressWarnings(file.remove(gsub(".sh$", paste0("_", my_comparison, ".out.txt"), pbs_file)))
    }
    if(opt$verbose) cat("\n")
  }
}
if(opt$verbose) cat("--------------------------------------\n")
if(opt$verbose) cat("PBS files at:", paste0(getwd(), "/scripts/"), "\n")
