#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=30-00:00:00
#SBATCH --job-name=pp_run
#SBATCH --mem=120g
#SBATCH --partition=naimi
#SBATCH --array=1-5

module purge
module load R

# Simulation Run Script
Rscript --no-save --no-restore --verbose ./code/sim_run.R 10 100 "hello_world\!"  > ./misc/sim_run.Rout 2>&1