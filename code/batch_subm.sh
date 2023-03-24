#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=30-00:00:00
#SBATCH --job-name=het_sim_test
#SBATCH --mem=120g
#SBATCH --partition=naimi
#SBATCH --mail-type=ALL

module purge
module load R

# Simulation Run Script
Rscript --no-save --no-restore --verbose /home/anaimi/Heterogeneity-Sim/code/sim_run_grf.R 10000 1000  > /home/anaimi/Heterogeneity-Sim/misc/sim_run_grf.Rout 2>&1
