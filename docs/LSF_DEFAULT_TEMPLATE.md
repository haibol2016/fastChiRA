# Default LSF Template (lsf-simple) Content

This document shows the content of the default LSF template used by batchtools when you specify `template = "lsf-simple"`.

## Default Template Content

The `lsf-simple.tmpl` template in batchtools typically contains:

```bash
#!/bin/bash
#BSUB <%= resources %>
#BSUB -J <%= job.name %>
#BSUB -o <%= log.file %>
#BSUB -e <%= log.file %>
<%= rscript %>
```

## Template Variables

Batchtools replaces these variables when submitting jobs:

- `<%= resources %>` - Resource requirements string (queue, cores, memory, walltime)
  - Example: `-q long -n 8 -R "rusage[mem=2GB]" -W 240:00`
  
- `<%= job.name %>` - Job name
  - Example: `chira_bt_output_1234567890_chunk`
  
- `<%= log.file %>` - Log file path (both stdout and stderr)
  - Example: `/path/to/registry/logs/job_1.log`
  
- `<%= rscript %>` - R script command to execute
  - Example: `Rscript /path/to/script.R`

## How to View the Actual Template on Your System

If you have batchtools installed, you can view the template by running:

```bash
# Method 1: Using Rscript
Rscript -e "cat(readLines(system.file('templates', 'lsf-simple.tmpl', package='batchtools')), sep='\n')"

# Method 2: Find and cat the file directly
Rscript -e "cat(system.file('templates', 'lsf-simple.tmpl', package='batchtools'), '\n')"
# Then: cat <path_from_above>

# Method 3: Use the provided script
Rscript print_lsf_template.R
```

## What Gets Generated

When batchtools submits a job, it generates a script like this:

```bash
#!/bin/bash
#BSUB -q long -n 8 -R "rusage[mem=2GB]" -W 240:00
#BSUB -J chira_bt_output_1234567890_chunk
#BSUB -o /path/to/registry/logs/job_1.log
#BSUB -e /path/to/registry/logs/job_1.log
Rscript /path/to/process_chunk_batchtools.R args...
```

## Custom Template

If you need to customize the template (e.g., add module loads), see `lsf_custom.tmpl` for an example.

