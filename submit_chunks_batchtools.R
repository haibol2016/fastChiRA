#!/usr/bin/env Rscript
# R script to submit chunk processing jobs using batchtools for LSF cluster
# Usage: Rscript submit_chunks_batchtools.R <config_json> <chunks_json>
#
# This script uses batchtools to submit LSF jobs for each chunk.
# Each job runs process_chunk_batchtools.py to process one chunk.

suppressPackageStartupMessages({
  library(batchtools)
  library(jsonlite)
})

# batchtools getJobTable() may return done/running/error/expired as POSIXct (timestamps), not logical.
# sum() on POSIXct fails; count non-NA for timestamps, sum(..., na.rm=TRUE) for logical.
count_status <- function(x) {
  if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) sum(!is.na(x)) else sum(as.logical(x), na.rm = TRUE)
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript submit_chunks_batchtools.R <config_json> <chunks_json>")
}

config_file <- args[1]
chunks_file <- args[2]

# Verify files exist and are readable
if (!file.exists(config_file)) {
  stop(sprintf("ERROR: Config file does not exist: %s", config_file))
}
if (!file.exists(chunks_file)) {
  stop(sprintf("ERROR: Chunks file does not exist: %s", chunks_file))
}

# Read configuration with error handling
tryCatch({
  config <- fromJSON(config_file)
}, error = function(e) {
  cat(sprintf("ERROR: Failed to parse config JSON file: %s\n", config_file))
  cat(sprintf("Error message: %s\n", as.character(e)))
  cat("First 500 characters of file:\n")
  tryCatch({
    file_content <- readLines(config_file, n = 20, warn = FALSE)
    cat(paste(file_content, collapse = "\n"))
    cat("\n")
  }, error = function(e2) {
    cat("Could not read file content\n")
  })
  stop(e)
})

tryCatch({
  chunks <- fromJSON(chunks_file)
}, error = function(e) {
  cat(sprintf("ERROR: Failed to parse chunks JSON file: %s\n", chunks_file))
  cat(sprintf("Error message: %s\n", as.character(e)))
  stop(e)
})

# Extract configuration
reg_dir <- config$reg_dir
queue <- config$queue
cores_per_job <- config$cores_per_job
memory_per_job <- config$memory_per_job
walltime <- config$walltime
conda_env <- config$conda_env
python_script <- config$python_script
chunk_dir <- config$chunk_dir
alignment_job_types_json <- config$alignment_job_types_json
per_chunk_processes <- config$per_chunk_processes
job_name_prefix <- config$job_name_prefix
max_parallel <- ifelse(is.null(config$max_parallel), nrow(chunks), config$max_parallel)

# Create batchtools registry
# Remove existing directory if it exists (following InPAS pattern)
if (dir.exists(reg_dir)) {
  unlink(reg_dir, recursive = TRUE, force = TRUE)
}
# Create registry (following InPAS pattern: conf.file = NA, seed = 1)
reg <- makeRegistry(
  file.dir = reg_dir,
  conf.file = NA,
  work.dir = getwd(),
  seed = 1
)

# Configure LSF cluster functions
# Python passes absolute path to template_file, or NULL/empty if using default
template_file <- if (is.null(config$template_file) || config$template_file == "") {
  # Default: look for lsf_custom.tmpl in current directory (where R script runs)
  # Python should have set this, but fallback here for safety
  default_template <- "lsf_custom.tmpl"
  if (file.exists(default_template)) {
    normalizePath(default_template)  # Use absolute path
  } else {
    "lsf-simple"  # Fallback to built-in
  }
} else {
  # Use template file path from Python (should be absolute path)
  config$template_file
}

# Configure LSF cluster functions (following InPAS pattern with latency settings)
# InPAS uses: makeClusterFunctionsLSF(template = template_file, scheduler.latency = 1, fs.latency = 65)
if (template_file == "lsf-simple" || file.exists(template_file)) {
  cat(sprintf("Configuring LSF cluster functions with template: %s\n", template_file))
  reg$cluster.functions <- makeClusterFunctionsLSF(
    template = template_file,
    scheduler.latency = 1,  # Latency for scheduler operations (seconds) - from InPAS
    fs.latency = 65         # Latency for filesystem operations (seconds) - from InPAS
  )
  cat(sprintf("LSF cluster functions configured successfully.\n"))
  
  # Verify cluster functions are set correctly
  if (is.null(reg$cluster.functions)) {
    stop("ERROR: Cluster functions are NULL after configuration!")
  }
  cat(sprintf("Cluster functions type: %s\n", class(reg$cluster.functions)[1]))
  
  # Test LSF connectivity
  cat("Testing LSF connectivity...\n")
  test_result <- tryCatch({
    system("bqueues > /dev/null 2>&1", intern = FALSE)
  }, error = function(e) {
    return(1)
  })
  if (test_result != 0) {
    cat("WARNING: Cannot execute 'bqueues'. LSF may not be available or not in PATH.\n")
    cat("  This may prevent jobs from being submitted to LSF.\n")
  } else {
    cat("LSF connectivity test passed.\n")
  }
} else {
  stop(sprintf("Template file '%s' does not exist. Use --batchtools_template to specify a valid template file.", template_file))
}

# Define job function
process_chunk_job <- function(chunk_file, chunk_idx, chunk_dir, 
                              alignment_job_types_json, per_chunk_processes, 
                              python_script, conda_env) {
  # Build command
  cmd <- paste(
    "python3", python_script,
    shQuote(chunk_file),
    chunk_idx,
    shQuote(chunk_dir),
    shQuote(alignment_job_types_json),
    per_chunk_processes
  )
  
  # Add conda activation if specified
  if (!is.null(conda_env) && conda_env != "") {
    conda_init <- 'eval "$(conda shell.bash hook)" 2>/dev/null || source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true'
    conda_activate <- paste('conda activate', conda_env)
    cmd <- paste(conda_init, "&&", conda_activate, "&&", cmd)
  }
  
  # Execute command
  system(cmd)
}

# Submit jobs
cat(sprintf("Submitting %d chunk jobs to LSF cluster...\n", nrow(chunks)))
if (max_parallel < nrow(chunks)) {
  cat(sprintf("NOTE: max_parallel=%d limits concurrent RUNNING jobs. Jobs will be submitted in batches of %d.\n", max_parallel, max_parallel))
  cat(sprintf("      The script will wait for each batch to complete before submitting the next batch.\n"))
  cat(sprintf("      This ensures only %d jobs are running at any time.\n", max_parallel))
} else {
  cat(sprintf("NOTE: No max_parallel limit. All %d jobs will be submitted at once.\n", nrow(chunks)))
}

# Prepare job resources
# Convert walltime from "HH:MM" format to seconds (batchtools expects seconds)
walltime_seconds <- if (grepl(":", walltime)) {
  parts <- strsplit(walltime, ":")[[1]]
  as.numeric(parts[1]) * 3600 + as.numeric(parts[2]) * 60
} else {
  # Assume already in seconds or minutes - try to parse
  as.numeric(walltime) * 60  # Assume minutes if no colon
}

# Map resources to template variable names
# Template uses: resources$ncpus, resources$mpp, resources$queue, resources$walltime
resources <- list(
  ncpus = cores_per_job,        # Number of CPUs (maps to -n)
  mpp = memory_per_job,          # Memory per processor (maps to rusage[mem=...])
  queue = queue,                 # Queue name (maps to -q)
  walltime = walltime_seconds,   # Walltime in seconds (template converts to minutes)
  job.name = paste0(job_name_prefix, "_chunk")
)

# Submit jobs
# batchtools batchMap() supports more.args parameter for constant arguments
ids <- batchMap(
  fun = process_chunk_job,
  chunk_file = chunks$chunk_file,
  chunk_idx = chunks$chunk_idx,
  more.args = list(
    chunk_dir = chunk_dir,
    alignment_job_types_json = alignment_job_types_json,
    per_chunk_processes = per_chunk_processes,
    python_script = python_script,
    conda_env = conda_env
  ),
  reg = reg
)

# Submit jobs with max_parallel constraint
# If max_parallel is set, submit jobs in batches and wait for completion before submitting more
if (max_parallel < length(ids)) {
  # Submit in batches, waiting for jobs to complete before submitting next batch
  num_batches <- ceiling(length(ids) / max_parallel)
  all_submitted_ids <- c()
  remaining_ids <- ids
  
  for (batch_idx in 1:num_batches) {
    # Get next batch of jobs to submit
    batch_size <- min(max_parallel, length(remaining_ids))
    batch_ids <- remaining_ids[1:batch_size]
    remaining_ids <- remaining_ids[-(1:batch_size)]
    
    cat(sprintf("Submitting batch %d/%d: %d jobs (max %d concurrent)...\n", 
                batch_idx, num_batches, length(batch_ids), max_parallel))
    
    # Submit jobs and check for errors
    tryCatch({
      submitJobs(batch_ids, resources = resources, reg = reg)
      
      # Verify jobs were actually submitted by checking status
      Sys.sleep(2)  # Give LSF a moment to register jobs
      job_table <- getJobTable(ids = batch_ids, reg = reg)
      
      # Check if jobs have batch.id (LSF job IDs) - this confirms they were submitted
      submitted_count <- sum(!is.na(job_table$batch.id))
      cat(sprintf("  Verified: %d/%d jobs have LSF batch IDs\n", submitted_count, length(batch_ids)))
      
      if (submitted_count == 0) {
        cat("  WARNING: No jobs have LSF batch IDs! Jobs may not have been submitted to LSF.\n")
        cat("  Check batchtools logs in:", reg_dir, "\n")
        cat("  Try: bjobs -u $USER\n")
      } else {
        # Print LSF job IDs for verification
        lsf_job_ids <- job_table$batch.id[!is.na(job_table$batch.id)]
        cat(sprintf("  LSF job IDs: %s\n", paste(lsf_job_ids, collapse=", ")))
        cat(sprintf("  Verify with: bjobs %s\n", paste(lsf_job_ids, collapse=" ")))
      }
      
      all_submitted_ids <- c(all_submitted_ids, batch_ids)
    }, error = function(e) {
      cat(sprintf("  ERROR submitting batch %d: %s\n", batch_idx, as.character(e)))
      cat("  Check batchtools configuration and LSF template.\n")
      stop(e)
    })
    
    # Wait for this batch to complete before submitting next batch (if any remain)
    # Use manual polling (avoid waitForJobs() which can trigger sum() on POSIXct in some batchtools versions)
    if (length(remaining_ids) > 0) {
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
        Sys.sleep(30)
      }
      job_table <- getJobTable(ids = batch_ids, reg = reg)
      batch_done <- count_status(job_table$done)
      batch_error <- count_status(job_table$error)
      batch_expired <- count_status(job_table$expired)
      cat(sprintf("Batch %d completed (%d done, %d error, %d expired). Proceeding to next batch.\n",
                  batch_idx, batch_done, batch_error, batch_expired))
    }
  }
  submitted_ids <- all_submitted_ids
} else {
  # Submit all jobs at once (no limit)
  cat(sprintf("Submitting all %d jobs at once (no max_parallel limit)...\n", length(ids)))
  
  # Submit jobs and check for errors
  tryCatch({
    submitJobs(ids, resources = resources, reg = reg)
    
    cat(sprintf("Waiting for all jobs to complete...\n"))
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
      Sys.sleep(30)
    }
    # Verify jobs were actually submitted by checking status
    Sys.sleep(2)  # Give LSF a moment to register jobs
    job_table <- getJobTable(ids = ids, reg = reg)
    
    # Check if jobs have batch.id (LSF job IDs) - this confirms they were submitted
    submitted_count <- sum(!is.na(job_table$batch.id))
    cat(sprintf("Verified: %d/%d jobs have LSF batch IDs\n", submitted_count, length(ids)))
    
    if (submitted_count == 0) {
      cat("WARNING: No jobs have LSF batch IDs! Jobs may not have been submitted to LSF.\n")
      cat("Check batchtools logs in:", reg_dir, "\n")
      cat("Try: bjobs -u $USER\n")
      cat("Check template file:", template_file, "\n")
    } else {
      # Print first few LSF job IDs for verification
      lsf_job_ids <- job_table$batch.id[!is.na(job_table$batch.id)]
      cat(sprintf("LSF job IDs (first 5): %s\n", paste(head(lsf_job_ids, 5), collapse=", ")))
      cat(sprintf("Verify with: bjobs %s\n", paste(head(lsf_job_ids, 5), collapse=" ")))
    }
    
    submitted_ids <- ids
  }, error = function(e) {
    cat(sprintf("ERROR submitting jobs: %s\n", as.character(e)))
    cat("Check batchtools configuration and LSF template.\n")
    stop(e)
  })
}

# Get final status and LSF job IDs
final_job_table <- getJobTable(ids = submitted_ids, reg = reg)
lsf_job_ids <- final_job_table$batch.id[!is.na(final_job_table$batch.id)]

cat(sprintf("\n=== Submission Summary ===\n"))
cat(sprintf("Batchtools job IDs: %s\n", paste(submitted_ids, collapse = ", ")))
cat(sprintf("LSF job IDs: %s\n", ifelse(length(lsf_job_ids) > 0, paste(lsf_job_ids, collapse = ", "), "NONE (jobs not submitted to LSF!)")))
cat(sprintf("Registry directory: %s\n", reg_dir))
cat(sprintf("\n=== Verification Commands ===\n"))
if (length(lsf_job_ids) > 0) {
  cat(sprintf("Check LSF jobs: bjobs %s\n", paste(head(lsf_job_ids, 10), collapse=" ")))
  cat(sprintf("Check all your jobs: bjobs -u $USER\n"))
  cat(sprintf("Monitor jobs by name: bjobs -J %s*\n", job_name_prefix))
} else {
  cat("WARNING: No LSF job IDs found! Jobs were not submitted to LSF.\n")
  cat("Troubleshooting steps:\n")
  cat("  1. Check batchtools logs in:", reg_dir, "\n")
  cat("  2. Verify LSF template:", template_file, "\n")
  cat("  3. Check if batchtools can access LSF: bqueues\n")
  cat("  4. Check registry status: loadRegistry('", reg_dir, "'); getStatus(reg)\n", sep="")
}
cat(sprintf("Check queue limits: bqueues -l %s\n", queue))
cat(sprintf("Check status in R: loadRegistry('%s'); getStatus(reg)\n", reg_dir))
cat(sprintf("\nNOTE: If jobs are not visible in bjobs, check:\n"))
cat(sprintf("  1. LSF queue limits: bqueues -l %s (look for 'MAX' or 'JOB_LIMIT')\n", queue))
cat(sprintf("  2. Your job limits: bqueues -l %s | grep -i limit\n", queue))
cat(sprintf("  3. All your jobs: bjobs -u $USER (some may be PEND)\n"))
cat(sprintf("  4. Batchtools logs: %s/logs/\n", reg_dir))

# Write job IDs to file for Python to read (all IDs, not just submitted)
writeLines(as.character(ids), file.path(reg_dir, "job_ids.txt"))

# Write registry path for Python
writeLines(reg_dir, file.path(reg_dir, "registry_path.txt"))

