#!/usr/bin/env bash
#PBS -V
#PBS -d .
#PBS -t 1
#PBS -q condo
#PBS -j oe
#PBS -o /dev/null
#PBS -N run.snakemake
#PBS -l nodes=1:ppn=3
#PBS -l walltime=8:00:00


# An example bash script demonstrating how to run the entire snakemake pipeline
# This script creates two separate log files in the output dir:
# 	1) log - the basic snakemake log of completed rules
# 	2) qlog - a more detailed log of the progress of each rule and any errors

# Before running this snakemake pipeline, remember to complete the config file
# with the required input info.
# Also, make sure that this script is executed from the directory that it lives in!

out_path="out"
mkdir -p "$out_path/logs"

# clear leftover log files
echo ""> "$out_path/logs/log"
echo ""> "$out_path/logs/qlog"

# try to find and activate the snakemake conda env if we need it
if ! command -v 'snakemake' &>/dev/null && \
	command -v 'conda' &>/dev/null && \
   [ "$CONDA_DEFAULT_ENV" != "snakemake" ] && \
   conda info --envs | grep "$CONDA_ROOT/snakemake" &>/dev/null; then
        echo "Snakemake not detected. Attempting to switch to snakemake environment." >> "out/logs/log"
        eval "$(conda shell.bash hook)"
        conda activate snakemake
fi


# check: are we being executed from within qsub?
if [ "$ENVIRONMENT" = "BATCH" ]; then
    snakemake \
    --cluster "qsub -d . -V -q condo -l walltime=00:30:00 -l nodes=1:ppn={threads} -j oe -o \"$out_path/logs/qlog\"" \
    --latency-wait 60 \
    --use-conda \
    --conda-frontend conda \
    -k \
    -j 12 \
    -c 12 \
    "$@" >>"$out_path/logs/log" 2>>"$out_path/logs/qlog"
else
    snakemake \
    --latency-wait 60 \
    --use-conda \
    --conda-frontend conda \
    -k \
    -c 12 \
    "$@" 2>>"$out_path/logs/log" >>"$out_path/logs/qlog"
fi

exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
        slack "snakemake finished successfully" &>/dev/null
    else
        slack "snakemake simulate_gwas job failed" &>/dev/null
        slack "$(tail -n4 "$out_path/logs/log")" &>/dev/null
    fi
fi
exit "$exit_code"
