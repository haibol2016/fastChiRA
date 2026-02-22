#!/usr/bin/env Rscript
# R script to submit IntaRNA hybridization jobs using batchtools (LSF cluster)
# Usage: Rscript submit_intarna_batchtools.R <config_json> <jobs_json>
#
# Follows the same pattern as submit_chunks_batchtools.R (chira_map.py).
# Each job runs IntaRNA once per chunk: multi-FASTA query.fa + target.fa -> result.csv (default). With --intarna_per_pair_only, runs once per pair instead.

suppressPackageStartupMessages({
  library(batchtools)
  library(jsonlite)
})

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
try_multi_fasta_first <- if (!is.null(config$intarna_try_multi_fasta_first)) config$intarna_try_multi_fasta_first else TRUE
max_parallel <- ifelse(is.null(config$max_parallel), nrow(jobs), config$max_parallel)
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

if (template_file == "lsf-simple" || file.exists(template_file)) {
  cat(sprintf("Configuring LSF cluster functions with template: %s\n", template_file))
  reg$cluster.functions <- makeClusterFunctionsLSF(
    template = template_file,
    scheduler.latency = 1,
    fs.latency = 65
  )
  if (is.null(reg$cluster.functions)) {
    stop("ERROR: Cluster functions are NULL after configuration!")
  }
} else {
  stop(sprintf("Template file '%s' does not exist. Use --batchtools_template to specify a valid template.", template_file))
}

# Parse multi-FASTA into list of {id, seq} in order (same order as pairs.tsv)
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

# Run IntaRNA. Default: one multi-FASTA run per chunk (query.fa + target.fa -> result.csv).
# If try_multi_fasta_first is FALSE (--intarna_per_pair_only), run once per pair and merge CSV.
run_intarna_job <- function(n, query_fa, target_fa, output_csv, params, conda_env = NULL, try_multi_fasta_first = TRUE) {
  chunk_dir <- dirname(query_fa)
  pairs_tsv <- file.path(chunk_dir, "pairs.tsv")
  if (!file.exists(pairs_tsv)) {
    stop(sprintf("pairs.tsv not found: %s", pairs_tsv))
  }
  pairs <- read.table(pairs_tsv, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  if (nrow(pairs) == 0L) {
    file.create(output_csv)
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

  # Default: one multi-FASTA run per chunk (same as before; result.csv filled by IntaRNA)
  if (try_multi_fasta_first) {
    ret <- run_one_intarna(query_fa, target_fa, output_csv)
    if (!is.null(attr(ret, "status")) && attr(ret, "status") != 0) {
      stop(sprintf("IntaRNA failed for chunk %s: %s", n, paste(ret, collapse = "\n")))
    }
    return(invisible(NULL))
  }

  # Optional: one IntaRNA run per pair (use --intarna_per_pair_only if multi-FASTA fails in your setup)
  query_blocks <- parse_fasta_blocks(query_fa)
  target_blocks <- parse_fasta_blocks(target_fa)
  if (length(query_blocks) != nrow(pairs) || length(target_blocks) != nrow(pairs)) {
    stop(sprintf("Chunk %s: pair count mismatch (pairs=%d, query blocks=%d, target blocks=%d)",
                 n, nrow(pairs), length(query_blocks), length(target_blocks)))
  }
  out_con <- file(output_csv, "w")
  on.exit(close(out_con), add = TRUE)
  header_written <- FALSE
  for (i in seq_len(nrow(pairs))) {
    q_tmp <- tempfile(pattern = "q_", fileext = ".fa")
    t_tmp <- tempfile(pattern = "t_", fileext = ".fa")
    tmp_out <- tempfile(pattern = "out_", fileext = ".csv")
    writeLines(c(paste0(">", query_blocks[[i]]$id), query_blocks[[i]]$seq), q_tmp)
    writeLines(c(paste0(">", target_blocks[[i]]$id), target_blocks[[i]]$seq), t_tmp)
    ret <- run_one_intarna(q_tmp, t_tmp, tmp_out)
    unlink(c(q_tmp, t_tmp))
    if (!is.null(attr(ret, "status")) && attr(ret, "status") != 0) {
      unlink(tmp_out)
      stop(sprintf("IntaRNA failed for pair %d in chunk %s: %s", i, n, paste(ret, collapse = "\n")))
    }
    if (file.exists(tmp_out)) {
      lines_out <- readLines(tmp_out, warn = FALSE)
      if (!header_written) {
        writeLines(lines_out, out_con)
        header_written <- TRUE
      } else {
        if (length(lines_out) > 1L) writeLines(lines_out[-1L], out_con)
      }
      unlink(tmp_out)
    }
  }
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
  more.args = list(params = intarna_params, conda_env = conda_env, try_multi_fasta_first = try_multi_fasta_first),
  reg = reg
)

if (max_parallel < length(ids)) {
  num_batches <- ceiling(length(ids) / max_parallel)
  remaining_ids <- ids
  all_submitted <- c()

  for (batch_idx in 1:num_batches) {
    batch_size <- min(max_parallel, length(remaining_ids))
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

    # Wait for this batch to complete before submitting next batch (if any remain)
    # Following InPAS pattern: use waitForJobs() for proper batchtools integration
    cat(sprintf("Waiting for batch %d to complete before submitting next batch...\n", batch_idx))
    # waitForJobs() returns TRUE when all jobs are done, FALSE if timeout/error
    while (!waitForJobs(ids = batch_ids, sleep = 30, timeout = Inf,
                        stop.on.error = FALSE, reg = reg)) {
      # Get status for progress reporting
      job_table <- getJobTable(ids = batch_ids, reg = reg)
      batch_done <- sum(job_table$done, na.rm = TRUE)
      batch_running <- sum(job_table$running, na.rm = TRUE)
      batch_error <- sum(job_table$error, na.rm = TRUE)
      batch_expired <- sum(job_table$expired, na.rm = TRUE)
      completed <- batch_done + batch_error + batch_expired

      cat(sprintf("  Batch %d: %d/%d completed (%d done, %d running, %d error). Waiting...\n",
                  batch_idx, completed, length(batch_ids), batch_done, batch_running, batch_error))
      Sys.sleep(30)  # Additional sleep for progress reporting
    }
    # Final status report
    job_table <- getJobTable(ids = batch_ids, reg = reg)
    batch_done <- sum(job_table$done, na.rm = TRUE)
    batch_error <- sum(job_table$error, na.rm = TRUE)
    batch_expired <- sum(job_table$expired, na.rm = TRUE)
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
  while (!waitForJobs(ids = ids, sleep = 30, timeout = Inf, stop.on.error = FALSE, reg = reg)) {
    job_table <- getJobTable(ids = ids, reg = reg)
    batch_done <- sum(job_table$done, na.rm = TRUE)
    batch_running <- sum(job_table$running, na.rm = TRUE)
    batch_error <- sum(job_table$error, na.rm = TRUE)
    batch_expired <- sum(job_table$expired, na.rm = TRUE)
    completed <- batch_done + batch_error + batch_expired
    cat(sprintf("  All jobs: %d/%d completed (%d done, %d running, %d error). Waiting...\n",
                completed, length(ids), batch_done, batch_running, batch_error))
    Sys.sleep(30)
  }
  job_table <- getJobTable(ids = ids, reg = reg)
  cat(sprintf("All jobs completed (%d done, %d error, %d expired).\n",
              sum(job_table$done, na.rm = TRUE), sum(job_table$error, na.rm = TRUE), sum(job_table$expired, na.rm = TRUE)))
}

job_table <- getJobTable(ids = submitted_ids, reg = reg)
lsf_ids <- job_table$batch.id[!is.na(job_table$batch.id)]
cat(sprintf("\nRegistry: %s\n", reg_dir))
cat(sprintf("LSF job IDs: %s\n", ifelse(length(lsf_ids) > 0, paste(lsf_ids, collapse = ", "), "NONE")))

writeLines(as.character(ids), file.path(reg_dir, "job_ids.txt"))
writeLines(reg_dir, file.path(reg_dir, "registry_path.txt"))
