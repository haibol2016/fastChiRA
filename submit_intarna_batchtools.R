#!/usr/bin/env Rscript
# R script to submit IntaRNA hybridization jobs using batchtools (LSF cluster)
# Usage: Rscript submit_intarna_batchtools.R <config_json> <jobs_json>
#
# Follows the same pattern as submit_chunks_batchtools.R (chira_map.py) and
# InPAS-style batchtools (e.g. 12.search_CPs.R: makeRegistry, batchMap, LSF template, polling).
# One batchtools job per chunk. Each job runs IntaRNA once per locus pair (parallel within job via future_lapply).

# Polling interval (seconds) while waiting for jobs; config can override via poll_sleep_seconds
.DEFAULT_POLL_SLEEP_SECONDS <- 120L

suppressPackageStartupMessages({
  library(batchtools)
  library(jsonlite)
})

# batchtools getJobTable() may return done/running/error/expired as POSIXct (timestamps), not logical.
# sum() on POSIXct fails; count non-NA for timestamps, sum(..., na.rm=TRUE) for logical.
count_status <- function(x) {
  if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) sum(!is.na(x)) else sum(as.logical(x), na.rm = TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript submit_intarna_batchtools.R <config_json> <jobs_json>")
}

config_file <- args[1]
jobs_file <- args[2]

if (!file.exists(config_file)) {
  stop(sprintf("ERROR: Config file does not exist: %s", config_file))
}
if (!file.exists(jobs_file)) {
  stop(sprintf("ERROR: Jobs file does not exist: %s", jobs_file))
}

tryCatch({
  config <- fromJSON(config_file)
}, error = function(e) {
  cat(sprintf("ERROR: Failed to parse config JSON: %s\n", config_file))
  cat(sprintf("Error: %s\n", as.character(e)))
  stop(e)
})

tryCatch({
  jobs <- fromJSON(jobs_file)
}, error = function(e) {
  cat(sprintf("ERROR: Failed to parse jobs JSON: %s\n", jobs_file))
  cat(sprintf("Error: %s\n", as.character(e)))
  stop(e)
})

reg_dir <- config$reg_dir
queue <- config$queue
cores_per_job <- config$cores_per_job
memory_per_job <- config$memory_per_job
walltime <- config$walltime
conda_env <- config$conda_env
intarna_params <- config$intarna_params
job_name_prefix <- config$job_name_prefix
max_parallel <- ifelse(is.null(config$max_parallel), nrow(jobs), config$max_parallel)
POLL_SLEEP_SECONDS <- if (is.null(config$poll_sleep_seconds)) .DEFAULT_POLL_SLEEP_SECONDS else max(10L, as.integer(config$poll_sleep_seconds)[1L])
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

# Workers load these R packages when running the job; install in the same conda env as IntaRNA (e.g. conda env used via --batchtools_conda_env).
reg <- makeRegistry(
  file.dir = reg_dir,
  conf.file = NA,
  work.dir = getwd(),
  seed = 1,
  packages = c("future", "future.apply")
)

# InPAS-style: require template to exist when not using built-in lsf-simple (see 12.search_CPs.R cluster_type == "lsf")
if (template_file != "lsf-simple" && !file.exists(template_file)) {
  stop(sprintf("Template file '%s' does not exist. Use --batchtools_template to specify a valid template.", template_file))
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

# Run IntaRNA once per locus pair (only real chimeric pairs; no all-vs-all).
# Creates result.csv with at least a header at job start so the file always exists (even if the job fails).
# Writes per-pair temp FASTA files and calls IntaRNA with -q / -t paths.
# Runs multiple pairs in parallel within the job when ncpus > 1 (uses future.apply::future_lapply).
# Helper and constants are defined inside so they are serialized with the job and available on batchtools workers.
run_intarna_job <- function(n, query_fa, target_fa, output_csv, params, conda_env = NULL, ncpus = 1L) {
  # Must be inside this function so batchtools workers have it when the job runs on the cluster
  parse_fasta_blocks <- function(path) {
    if (!file.exists(path)) return(list())
    lines <- readLines(path, warn = FALSE)
    starts <- which(grepl("^>", lines))
    blocks <- list()
    for (i in seq_along(starts)) {
      start <- starts[i]
      end <- if (i < length(starts)) starts[i + 1L] - 1L else length(lines)
      id <- sub("^>\\s*", "", lines[start])
      seq <- paste(lines[(start + 1L):end], collapse = "")
      blocks[[i]] <- list(id = id, seq = seq)
    }
    blocks
  }
  intarna_csv_header <- "id1;id2;start1;end1;start2;end2;hybridDPfull;E"

  chunk_dir <- dirname(query_fa)
  pairs_tsv <- file.path(chunk_dir, "pairs.tsv")

  # Create result.csv immediately with header so the file exists even if we stop() later (e.g. pairs_tsv not found on compute node).
  writeLines(intarna_csv_header, output_csv)

  if (!file.exists(pairs_tsv)) {
    stop(sprintf("pairs.tsv not found: %s (check that batchtools work dir is on a shared filesystem visible from compute nodes)", pairs_tsv))
  }
  pairs <- read.table(pairs_tsv, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  if (nrow(pairs) == 0L) {
    return(invisible(NULL))
  }

  base_args <- strsplit(params, " ")[[1]]
  base_args <- base_args[base_args != ""]
  run_one_intarna <- function(q_fa, t_fa, out_file) {
    args <- c(base_args, "-q", q_fa, "-t", t_fa, "--out", out_file)
    if (!is.null(conda_env) && conda_env != "") {
      conda_init <- 'eval "$(conda shell.bash hook)" 2>/dev/null || source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true'
      conda_activate <- paste("conda activate", shQuote(conda_env))
      intarna_cmd <- paste("IntaRNA", paste(sapply(args, shQuote), collapse = " "))
      cmd <- paste(conda_init, "&&", conda_activate, "&&", intarna_cmd)
      ret <- system(cmd, intern = TRUE)
    } else {
      ret <- system2("IntaRNA", args = args, stdout = TRUE, stderr = TRUE)
    }
    ret
  }

  query_blocks <- parse_fasta_blocks(query_fa)
  target_blocks <- parse_fasta_blocks(target_fa)
  if (length(query_blocks) != nrow(pairs) || length(target_blocks) != nrow(pairs)) {
    stop(sprintf("Chunk %s: pair count mismatch (pairs=%d, query blocks=%d, target blocks=%d)",
                 n, nrow(pairs), length(query_blocks), length(target_blocks)))
  }
  # Run IntaRNA in parallel (ncpus workers); use future.chunk.size so each worker gets multiple subjobs
  ncpus_int <- as.integer(ncpus)[1L]
  if (is.na(ncpus_int) || ncpus_int < 1L) ncpus_int <- 1L
  nworkers <- max(1L, min(ncpus_int, nrow(pairs)))
  pair_indices <- seq_len(nrow(pairs))
  # Average number of pairs per future (subjob); multiple subjobs per worker to reduce overhead
  n_pairs <- length(pair_indices)
  future_chunk_size <- max(1L, ceiling(n_pairs / (nworkers * 4L)))
  if (nworkers > 1L) {
    if (.Platform$OS.type == "unix") {
      future::plan(future::multicore, workers = nworkers)
    } else {
      future::plan(future::multisession, workers = nworkers)
    }
    on.exit(future::plan(future::sequential), add = TRUE)
  }
  results <- future.apply::future_lapply(pair_indices, function(i) {
    tmp_out <- tempfile(pattern = "out_", fileext = ".csv")
    on.exit(unlink(tmp_out, force = TRUE), add = TRUE)
    q_tmp <- tempfile(pattern = "q_", fileext = ".fa")
    t_tmp <- tempfile(pattern = "t_", fileext = ".fa")
    on.exit(unlink(c(q_tmp, t_tmp), force = TRUE), add = TRUE)
    writeLines(c(paste0(">", query_blocks[[i]]$id), query_blocks[[i]]$seq), q_tmp)
    writeLines(c(paste0(">", target_blocks[[i]]$id), target_blocks[[i]]$seq), t_tmp)
    ret <- run_one_intarna(q_tmp, t_tmp, tmp_out)
    if (!is.null(attr(ret, "status")) && attr(ret, "status") != 0) {
      stop(sprintf("IntaRNA failed for pair %d in chunk %s: %s", i, n, paste(ret, collapse = "\n")))
    }
    if (file.exists(tmp_out)) {
      lines_out <- readLines(tmp_out, warn = FALSE)
      if (length(lines_out) > 1L) lines_out[-1L] else character(0L)
    } else {
      character(0L)
    }
  }, future.seed = TRUE, future.globals = TRUE, future.chunk.size = future_chunk_size)
  # Write all data lines in order (header already written at start)
  out_con <- file(output_csv, "a")
  on.exit(close(out_con), add = TRUE)
  for (r in results) if (length(r) > 0L) writeLines(r, out_con)
  invisible(NULL)
}

cat(sprintf("Submitting %d IntaRNA jobs to LSF cluster...\n", nrow(jobs)))

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
  fun = run_intarna_job,
  n = jobs$n,
  query_fa = jobs$query_fa,
  target_fa = jobs$target_fa,
  output_csv = jobs$output_csv,
  more.args = list(params = intarna_params, conda_env = conda_env, ncpus = cores_per_job),
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

    # Wait for this batch by polling (avoid waitForJobs() which can trigger sum() on POSIXct in some batchtools versions)
    cat(sprintf("Waiting for batch %d to complete before submitting next batch...\n", batch_idx))
    repeat {
      job_table <- getJobTable(ids = batch_ids, reg = reg)
      batch_done <- count_status(job_table$done)
      batch_running <- count_status(job_table$running)
      batch_error <- count_status(job_table$error)
      batch_expired <- count_status(job_table$expired)
      completed <- batch_done + batch_error + batch_expired
      if (completed >= length(batch_ids)) break
      cat(sprintf("  Batch %d: %d/%d completed (%d done, %d running, %d error). Waiting...\n",
                  batch_idx, completed, length(batch_ids), batch_done, batch_running, batch_error))
      Sys.sleep(POLL_SLEEP_SECONDS)
    }
    job_table <- getJobTable(ids = batch_ids, reg = reg)
    batch_done <- count_status(job_table$done)
    batch_error <- count_status(job_table$error)
    batch_expired <- count_status(job_table$expired)
    cat(sprintf("Batch %d completed (%d done, %d error, %d expired). Proceeding to next batch.\n",
                batch_idx, batch_done, batch_error, batch_expired))
  }
  submitted_ids <- all_submitted
} else {
  # Submit all jobs at once (no limit)
  cat(sprintf("Submitting all %d jobs at once (no max_parallel limit)...\n", length(ids)))
  tryCatch({
    submitJobs(ids, resources = resources, reg = reg)
    submitted_ids <- ids
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", as.character(e)))
    stop(e)
  })

  cat("Waiting for all jobs to complete...\n")
  repeat {
    job_table <- getJobTable(ids = ids, reg = reg)
    batch_done <- count_status(job_table$done)
    batch_running <- count_status(job_table$running)
    batch_error <- count_status(job_table$error)
    batch_expired <- count_status(job_table$expired)
    completed <- batch_done + batch_error + batch_expired
    if (completed >= length(ids)) break
    cat(sprintf("  All jobs: %d/%d completed (%d done, %d running, %d error). Waiting...\n",
                completed, length(ids), batch_done, batch_running, batch_error))
    Sys.sleep(POLL_SLEEP_SECONDS)
  }
  job_table <- getJobTable(ids = ids, reg = reg)
  cat(sprintf("All jobs completed (%d done, %d error, %d expired).\n",
              count_status(job_table$done), count_status(job_table$error), count_status(job_table$expired)))
}

job_table <- getJobTable(ids = submitted_ids, reg = reg)
n_error <- count_status(job_table$error)
n_expired <- count_status(job_table$expired)
if (n_error > 0L || n_expired > 0L) {
  # Handle both logical and POSIXct columns (batchtools may return timestamps for error/expired)
  error_or_expired <- !is.na(job_table$error) | !is.na(job_table$expired)
  failed <- submitted_ids[error_or_expired]
  cat(sprintf("\nERROR: %d job(s) failed or expired. Failed job IDs: %s\n", n_error + n_expired, paste(failed, collapse = ", ")), file = stderr())
  cat("Check the batchtools registry logs (e.g. in batchtools_work/registry/) for details.\n", file = stderr())
  quit(status = 1L, save = "no")
}
lsf_ids <- job_table$batch.id[!is.na(job_table$batch.id)]
cat(sprintf("\nRegistry: %s\n", reg_dir))
cat(sprintf("LSF job IDs: %s\n", ifelse(length(lsf_ids) > 0, paste(lsf_ids, collapse = ", "), "NONE")))

writeLines(as.character(ids), file.path(reg_dir, "job_ids.txt"))
writeLines(reg_dir, file.path(reg_dir, "registry_path.txt"))
