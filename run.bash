#!/usr/bin/env bash
#PBS -V
#PBS -d .
#PBS -t 1
#PBS -q condo
#PBS -o /dev/null
#PBS -e /dev/null
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

# you can specify a directory for all output here:
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
    --cluster "qsub -d . -V -q condo -l walltime=00:30:00 -l nodes=1:ppn={threads} -o /dev/null -e $out_path/logs/qlog" \
    --config out="$out_path" \
    --latency-wait 60 \
    --use-conda \
    --conda-frontend conda \
    -k \
    -j 12 \
    "$@" &>"$out_path/logs/log"
else
    snakemake \
    --config out="$out_path" \
    --latency-wait 60 \
    --conda-frontend conda \
    -k \
    -j 12 \
    "$@" 2>>"$out_path/logs/log" >>"$out_path/logs/qlog"
fi
