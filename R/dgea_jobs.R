#!/usr/bin/R

optlist <- list(
  optparse::make_option(
    opt_str = c("-y", "--yaml"), type = "character", default = "config_project.yaml",
    help = "Configuration file."
  ),
  optparse::make_option(
    opt_str = c("-s", "--submit"), type = "logical", default = FALSE,
    help = "Submit jobs."
  ),
  optparse::make_option(
    opt_str = c("-p", "--params"), type = "character",
    help = "Override dgea.R parameters. Please always quote character types."
  ),
  optparse::make_option(
    opt_str = c("-v", "--verbose"), type = "logical", default = TRUE,
    help = "Verbose."
  )
)
optparse <- optparse::OptionParser(option_list = optlist)
opt <- optparse::parse_args(optparse)
if(interactive()){ # Example/manually
  opt$yaml = "https://raw.githubusercontent.com/vijaybioinfo/dgea/main/config.yaml"
}
config = yaml::read_yaml(opt$yaml)

if(requireNamespace("crayon", quietly = TRUE)){
  cyan = crayon::cyan; redb = crayon::red$bold; gren = crayon::green; yelo = crayon::yellow
}else{ cyan = redb = gren = yelo = c }

if(opt$verbose){
  cat(cyan("########################################################\n"))
  cat(redb("########### Differential Expression Analysis ###########\n"))
  cat(cyan("########################################################\n"))
}

if(opt$verbose) cat(cyan('\n------- Directories structure -----------------\n'))
output_dir <- paste0(sub("\\/$", "", config$output_dir), "/", config$project, "/")
if(!grepl("scratch|beegfs", getwd())){
  cat("WARNING: No scratch folder involved; careful about temp files.\n")
  dir.create(output_dir, recursive = TRUE); setwd(output_dir)
}
if(opt$verbose) cat(redb("Working in:", getwd(), "\n"))
dir.create("scripts", showWarnings = FALSE)
if(opt$verbose){
  cat("Configuration:")
  str(config[!names(config) %in% "job"], max.level = 2)
  cat("Working at:", getwd(), "\n")
  system("ls -loh")
  cat("\n")
}

if(opt$verbose) cat(cyan('\n------- Packages and functions ----------------\n'))
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
parameters2set_fun <- function(x, config_name, verbose = FALSE){
  parameters2set <- NULL
  if(isTRUE(grepl("\\-\\-", x))){
    if(verbose) cat("Passed as extras\n")
    tmp <- setdiff(unlist(strsplit(x, "\\-\\-")), "")
    parameters2set <- setNames(gsub(".*=", "", tmp), gsub("=.*", "", tmp))
  }
  if(!is.null(names(parameters2set))){
    for(i in 1:length(parameters2set)){
      param = paste0("config_name", "$", names(parameters2set)[i])
      command = paste0(param, " = ", parameters2set[[i]])
      if(is.null(eval(parse(text = param)))){
        if(verbose) cat(" ", command, "\n")
        eval(parse(text = command))#, envir = .GlobalEnv)
      }
    }
  }; return(config_name)
}

if(opt$verbose) cat(cyan('\n------- Digesting parameters ------------------\n'))
# opt$params = "--genesetf='pre:-PR[LS]' --ctrans='log2' --thresh_min=0.2"
config <- parameters2set_fun(opt$params, config, verbose = opt$verbose)

if(opt$verbose) cat("Getting job template\n")
template_pbs_con <- file(description = config$job$template, open = "r")
template_pbs <- readLines(con = template_pbs_con)
close(template_pbs_con)

# Setting lements to iterate through
if(opt$verbose) cat("Fetching comparisons\n")
if(is.null(config$comparisons$file)) config$comparisons$file = "none"
comparisons = config$comparisons[!names(config$comparisons) %in% "file"]
if(file.exists(config$comparisons$file)){
  if(opt$verbose) cat("- from file\n")
  comps_tab = read.csv(config$comparisons$file, stringsAsFactors = FALSE)
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
for(m in config$method){
  for(config_comp in comparisons){
    if(opt$verbose)
      cat(config_comp$contrast, "\n", config_comp$test_column, config_comp$context)
    my_comparison <- if(isTRUE(grepl("vs", config_comp$name))){
      config_comp$name
    }else{ paste0(config_comp$contrast, collapse = "vs") }
    res_file0 <- paste0(config_comp$context, "/", my_comparison, "/", m, "_results.csv")
    res_file <- paste0(getwd(), "/", res_file0)
    re_run <- isTRUE(any(grepl("^f$|force", c(config$job$submit, opt$submit))))
    if(file.exists(res_file) && !re_run){
      if(opt$verbose) cat(gren(" - done\n")); next
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
      config$exec,
      " ", config$script,
      " --method=", m,
      " --metadata=", config$metadata,
      " --expression_data=", config$expression_data,
      " --group1=", config_comp$contrast[1],
      " --group2=", config_comp$contrast[2],
      " --newnames=", my_comparison,
      " --hname=", config_comp$test_column,
      " --output_dir=", getwd(),
      " --filters=\"", my_filter, "\"",
      " --context=", ifelse(is.null(config_comp$context), "none", config_comp$context),
      " --covariates=", config$covariates,
      " --down_sample=", config$down_sample,
      " --genesetf=", ifelse(is.null(config$genesetf), "none", config$genesetf),
      " --ctrans=", config$ctrans,
      " --padj_threshold=", config$padj_threshold,
      " --fc_threshold=", config$fc_threshold,
      " --annotation=", config$annotation,
      " --thresh_min=", 0,
      " --add_df=", ifelse(is.null(config$add_df), "none", config$add_df),
      " --verbose=", TRUE,
      " "
    )
    if(!is.null(config_comp$code)){
      params <- paste0(params, '--code="', config_comp$code, '"')
    }

    # Output directory
    output_dir <- getwd()

    pbs <- gsub("cellranger", "dgea", template_pbs)
    pbs <- gsub("\\{username\\}", Sys.info()[["user"]], pbs)
    pbs <- gsub("\\{routine_pbs\\}", routine_pbs_file, pbs)
    pbs <- gsub("\\{sampleid\\}", my_comparison, pbs)
    pbs <- gsub("\\{outpath\\}", output_dir, pbs)
    pbs <- gsub("\\{routine_params\\}", params, pbs)
    if(file.exists(res_file)){
      # Copying partial results to continue from them
      pbs <- gsub("# \\{after_copy\\}", paste("mkdir --parents", dirname(res_file0)), pbs)
      tvar <- paste0("cp -R ", dirname(res_file), "/. ./", dirname(res_file0), "/")
      pbs <- gsub("# \\{pre_routine\\}", tvar, pbs)
    }
    tvar <- if(is.null(config_comp$job)) config$job$main else config_comp$job
    for(i in names(config$job$main)){
      job_parm <- if(i %in% names(tvar)) tvar[[i]] else config$job$main[[i]]
      job_parm <- if(m %in% names(job_parm)) job_parm[[m]] else job_parm[[1]]
      pbs <- gsub(paste0("\\{", i, "\\}"), job_parm, pbs)
    }
    pbs <- gsub("\\.\\.\\/scripts", "scripts", pbs)

    # Checking if it's currently on the run
    running <- running_jobs()
    do_we_submit <- any(c(opt$submit, re_run))
    job_name <- routine_pbs_file#paste0(routine_pbs_file, "_", my_comparison, "$")
    if(any(grepl(job_name, running$Name)) && isFALSE(do_we_submit)){
      if(opt$verbose) cat(" - running\n"); next
    }

    pbs_file <- paste0(getwd(), "/scripts/", routine_pbs_file, ".sh")
    if(opt$verbose) cat(" - creating job file");
    writeLines(text = pbs, con = pbs_file)

    # Ready? Check if you really really want to... and for dependencies
    if(isTRUE(any(c(config$job$submit, do_we_submit)))){
      depend <- if(isTRUE(config$job$depend %in% running$id))
        paste0("-W depend=afterok:", config$job$depend)
      pbs_command <- paste("qsub", depend, pbs_file)
      if(opt$verbose) cat("\n", pbs_command, "\n"); system(pbs_command)
      # system2(pbs_command, stdout = opt$verbose, stderr = opt$verbose)
      void <- suppressWarnings(file.remove(gsub(".sh$", paste0("_", my_comparison, ".out.txt"), pbs_file)))
      void <- suppressWarnings(file.remove(gsub("sh$", "out.txt", pbs_file)))
    }
    if(opt$verbose) cat("\n")
  }
}
if(opt$verbose) cat("--------------------------------------\n")
if(opt$verbose) cat("PBS files at:", paste0(getwd(), "/scripts/"), "\n")
