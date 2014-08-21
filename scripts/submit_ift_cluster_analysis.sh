#!/bin/sh

#SBATCH --job-name=ift_analysis
#SBATCH --partition=pre
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

#SBATCH --time=1-00:00:00

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8192

echo "Date:"
date '+%s'
echo "Using ACI / HCP / Slurm cluster."
echo "JobID = $SLURM_JOB_ID"
echo "Using $SLURM_NNODES nodes"
echo "Using $SLURM_NODELIST nodes."
echo "Number of cores per node: $SLURM_TASKS_PER_NODE"
echo "Submit directory: $SLURM_SUBMIT_DIR"
echo ""

# Executable
python /home/maldonis/model_analysis/scripts/ift_cluster_analysis.py $@

echo "Finished on:"
date '+%s'
