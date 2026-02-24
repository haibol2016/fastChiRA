#!/usr/bin/env Rscript
# R script to submit IntaRNA hybridization jobs using batchtools (LSF cluster)
# Usage: Rscript submit_intarna_batchtools.R <config_json> <chunks_json>
#
# Same implementation strategy as submit_chunks_batchtools.R (chira_map.py):
# - One LSF job per chunk; each job runs process_intarna_chunk_batchtools.py (Python)
# - Python worker reads query.fa and target.fa from chunk_dir, runs IntaRNA once with
#   --outPairwise (one result per same-index pair), writes result.csv

suppressPackageStartupMessages({
  library(batchtools)
  library(jsonlite)
})

count_status <- function(x) {
  if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) sum(!is.na(x)) else sum(as.logical(x), na.rm = TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript submit_intarna_batchtools.R <config_json> <chunks_json>")
}

config_file <- args[1]
chunks_file <- args[2]

if (!file.exists(config_file)) {
  stop(sprintf("ERROR: Config file does not exist: %s", config_file))
}
if (!file.exists(chunks_file)) {
  stop(sprintf("ERROR: Chunks file does not exist: %s", chunks_file))
}

tryCatch({
  config <- fromJSON(config_file)
}, error = function(e) {
  cat(sprintf("ERROR: Failed to parse config JSON: %s\n", config_file))
  stop(e)
})

tryCatch({
  chunks <- fromJSON(chunks_file)
}, error = function(e) {
  cat(sprintf("ERROR: Failed to parse chunks JSON: %s\n", chunks_file))
  stop(e)
})

# Ensure chunks is a data frame with chunk_idx
if (is.list(chunks) && !is.data.frame(chunks)) {
  chunks <- as.data.frame(do.call(rbind, lapply(chunks, as.list)))
}
if (!"chunk_idx" %in% names(chunks)) {
  stop("ERROR: chunks JSON must contain 'chunk_idx' for each chunk")
}

reg_dir <- config$reg_dir
queue <- config$queue
cores_per_job <- config$cores_per_job
memory_per_job <- config$memory_per_job
walltime <- config$walltime
conda_env <- config$conda_env
intarna_params <- config$intarna_params
job_name_prefix <- config$job_name_prefix
max_parallel <- ifelse(is.null(config$max_parallel), nrow(chunks), config$max_parallel)
batchtools_work_dir <- normalizePath(config$batchtools_work_dir, mustWork = FALSE)
python_script <- normalizePath(config$python_script, mustWork = FALSE)

template_file <- if (is.null(config$template_file) || config$template_file == "") {
  if (file.exists("lsf_custom.tmpl")) {
    normalizePath("lsf_custom.tmpl")
  } else {
    "lsf-simple"
  }
} else {
  config$template_file
}

if (dir.exists(reg_dir)) {
  unlink(reg_dir, recursive = TRUE, force = TRUE)
}

reg <- makeRegistry(
  file.dir = reg_dir,
  conf.file = NA,
  work.dir = getwd(),
  seed = 1
)

if (template_file != "lsf-simple" && !file.exists(template_file)) {
  stop(sprintf("Template file '%s' does not exist.", template_file))
}
cat(sprintf("Configuring LSF cluster functions with template: %s\n", template_file))
reg$cluster.functions <- makeClusterFunctionsLSF(
  template = template_file,
  scheduler.latency = 1,
  fs.latency = 65
)
if (is.null(reg$cluster.functions)) {
  stop("ERROR: Cluster functions are NULL after configuration!")
}

# Job function: run Python worker for one chunk (same pattern as process_chunk_job in submit_chunks_batchtools.R)
process_intarna_chunk_job <- function(chunk_idx, batchtools_work_dir, python_script,
                                      intarna_params, conda_env, ncpus) {
  chunk_dir <- file.path(batchtools_work_dir, as.character(chunk_idx))
  cmd <- paste(
    "python3", shQuote(python_script),
    shQuote(chunk_dir),
    shQuote(intarna_params),
    ncpus,
    if (!is.null(conda_env) && conda_env != "") shQuote(conda_env) else ""
  )
  if (!is.null(conda_env) && conda_env != "") {
    conda_init <- 'eval "$(conda shell.bash hook)" 2>/dev/null || source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true'
    conda_activate <- paste("conda activate", shQuote(conda_env))
    cmd <- paste(conda_init, "&&", conda_activate, "&&", cmd)
  }
  system(cmd)
}

walltime_seconds <- if (grepl(":", walltime)) {
  parts <- strsplit(walltime, ":")[[1]]
  as.numeric(parts[1]) * 3600 + as.numeric(parts[2]) * 60
} else {
  as.numeric(walltime) * 60
}

resources <- list(
  ncpus = cores_per_job,
  mpp = memory_per_job,
  queue = queue,
  walltime = walltime_seconds,
  job.name = paste0(job_name_prefix, "_intarna")
)

ids <- batchMap(
  fun = process_intarna_chunk_job,
  chunk_idx = chunks$chunk_idx,
  more.args = list(
    batchtools_work_dir = batchtools_work_dir,
    python_script = python_script,
    intarna_params = intarna_params,
    conda_env = conda_env,
    ncpus = cores_per_job
  ),
  reg = reg
)

if (max_parallel < length(ids)) {
  num_batches <- ceiling(length(ids) / max_parallel)
  remaining_ids <- ids
  all_submitted <- c()

  for (batch_idx in 1:num_batches) {
    batch_size <- min(max_parallel, length(remaining_ids))
    if (batch_size == 0L) break
    batch_ids <- remaining_ids[1:batch_size]
    remaining_ids <- remaining_ids[-(1:batch_size)]

    cat(sprintf("Submitting batch %d/%d: %d jobs...\n", batch_idx, num_batches, length(batch_ids)))
    tryCatch({
      submitJobs(batch_ids, resources = resources, reg = reg)
      all_submitted <- c(all_submitted, batch_ids)
    }, error = function(e) {
      cat(sprintf("ERROR submitting batch %d: %s\n", batch_idx, as.character(e)))
      stop(e)
    })

    cat(sprintf("Waiting for batch %d to complete...\n", batch_idx))
    while (!waitForJobs(batch_ids, timeout = walltime_seconds, progressbar = TRUE,
           sleep = 120, stop.on.error = TRUE, reg = reg)) {
      Sys.sleep(120)
    } 
    cat(sprintf("Batch %d completed.\n", batch_idx))
  }
  submitted_ids <- all_submitted
} else {
  cat(sprintf("Submitting all %d jobs at once...\n", length(ids)))
  tryCatch({
    submitJobs(ids, resources = resources, reg = reg)
    submitted_ids <- ids

    cat("Waiting for all jobs to complete...\n")
    while (!waitForJobs(ids, timeout = walltime_seconds, progressbar = TRUE,
            sleep = 120, stop.on.error = TRUE, reg = reg)) {
      Sys.sleep(120)
    }
    cat("All jobs completed.\n")
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", as.character(e)))
    stop(e)
  })

}

job_table <- getJobTable(ids = submitted_ids, reg = reg)
n_error <- count_status(job_table$error)
n_expired <- count_status(job_table$expired)
if (n_error > 0L || n_expired > 0L) {
  error_or_expired <- !is.na(job_table$error) | !is.na(job_table$expired)
  failed <- submitted_ids[error_or_expired]
  cat(sprintf("\nERROR: %d job(s) failed or expired. Failed job IDs: %s\n", n_error + n_expired, 
              paste(failed, collapse = ", ")), file = stderr())
  quit(status = 1L, save = "no")
}

lsf_ids <- job_table$batch.id[!is.na(job_table$batch.id)]
cat(sprintf("\nRegistry: %s\n", reg_dir))
cat(sprintf("LSF job IDs: %s\n", ifelse(length(lsf_ids) > 0, paste(lsf_ids, collapse = ", "), "NONE")))

writeLines(as.character(ids), file.path(reg_dir, "job_ids.txt"))
writeLines(reg_dir, file.path(reg_dir, "registry_path.txt"))
