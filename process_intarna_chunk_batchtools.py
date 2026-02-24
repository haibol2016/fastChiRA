#!/usr/bin/env python
"""
Standalone script to process a single IntaRNA chunk for batchtools job submission.
Called by submit_intarna_batchtools.R: one LSF job per chunk.
Reads query.fa and target.fa from chunk_dir; runs IntaRNA once with --outPairwise
(one result per same-index pair); writes result.csv.
"""
import os
import sys
import subprocess
import shlex

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

def process_chunk_standalone(chunk_dir, intarna_params, ncpus, conda_env=None):
    """
    Process one IntaRNA chunk: run IntaRNA once with -q query.fa -t target.fa --outPairwise
    (pairwise mode: one result per same-index query/target pair). Writes result.csv in chunk_dir.
    No subprocess timeout (batchtools job walltime limits runtime).
    """
    query_fa = os.path.join(chunk_dir, "query.fa")
    target_fa = os.path.join(chunk_dir, "target.fa")
    result_csv = os.path.join(chunk_dir, "result.csv")

    if not os.path.exists(query_fa) or not os.path.exists(target_fa):
        raise FileNotFoundError(f"query.fa or target.fa not found in {chunk_dir}")

    # If no sequences (empty FASTA), write header-only result
    with open(query_fa, 'r') as f:
        n_query = f.read().count('>')
    if n_query == 0:
        with open(result_csv, 'w') as f:
            f.write("id1;id2;start1;subseq1;start2;subseq2;hybridDPfull;E\n")
        return

    # Single IntaRNA run: pairwise mode (one result per query-target pair by index)
    cmd = (intarna_params.split() if isinstance(intarna_params, str) else list(intarna_params))
    cmd += ["-q", "query.fa", "-t", "target.fa", "--out", "result.csv"]

    if conda_env:
        conda_init = 'eval "$(conda shell.bash hook)" 2>/dev/null || source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true'
        conda_activate = f"conda activate {shlex.quote(conda_env)}"
        intarna_cmd = "IntaRNA " + " ".join(shlex.quote(a) for a in cmd)
        full_cmd = f"{conda_init} && {conda_activate} && {intarna_cmd}"
        result = subprocess.run(
            ["bash", "-c", full_cmd],
            capture_output=True,
            text=True,
            cwd=chunk_dir,
        )
    else:
        result = subprocess.run(
            ["IntaRNA"] + cmd,
            capture_output=True,
            text=True,
            cwd=chunk_dir,
        )

    if result.returncode != 0:
        raise RuntimeError(
            f"IntaRNA failed (exit {result.returncode}) in {chunk_dir}: "
            f"{result.stderr[:500] if result.stderr else result.stdout[:500]}"
        )


def main():
    # R passes: chunk_dir, intarna_params, ncpus; conda_env only when set (so 4 or 5 args total)
    if len(sys.argv) < 4:
        print(
            "Usage: process_intarna_chunk_batchtools.py <chunk_dir> <intarna_params> <ncpus> [conda_env]",
            file=sys.stderr,
        )
        sys.exit(1)

    chunk_dir = os.path.abspath(sys.argv[1])
    intarna_params = sys.argv[2]
    ncpus = int(sys.argv[3])
    conda_env = sys.argv[4].strip() if len(sys.argv) > 4 and sys.argv[4].strip() else os.environ.get("CONDA_DEFAULT_ENV", "") or ""

    if not os.path.isdir(chunk_dir):
        print(f"ERROR: chunk_dir is not a directory: {chunk_dir}", file=sys.stderr)
        sys.exit(1)

    try:
        process_chunk_standalone(chunk_dir, intarna_params, ncpus, conda_env=conda_env or None)
        print(f"SUCCESS: IntaRNA chunk completed: {chunk_dir}", file=sys.stderr)
        sys.exit(0)
    except Exception as e:
        print(f"FAILED: {chunk_dir}: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
