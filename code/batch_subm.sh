#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=30-00:00:00
#SBATCH --job-name=het_sim_test
#SBATCH --mem=120g
#SBATCH --partition=naimi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anaimi@emory.edu

module purge
module load R

# Simulation Run Script
for ii in {1..5} ; do
  Rscript --no-save --no-restore --verbose /home/anaimi/Heterogeneity-Sim/code/sim_run_grf.R 3 1000  > /home/anaimi/Heterogeneity-Sim/misc/sim_run_grf_3_{ii}.Rout 2>&1
  
  #Rscript --no-save --no-restore --verbose /home/anaimi/Heterogeneity-Sim/code/sim_run_grf.R 3 1000  > /home/anaimi/Heterogeneity-Sim/misc/sim_run_grf_3.Rout 2>&1
  
  #Rscript --no-save --no-restore --verbose /home/anaimi/Heterogeneity-Sim/code/sim_run_grf.R 9 1000  > /home/anaimi/Heterogeneity-Sim/misc/sim_run_grf_9.Rout 2>&1
  
  #Rscript --no-save --no-restore --verbose /home/anaimi/Heterogeneity-Sim/code/sim_run_grf.R 15 1000  > /home/anaimi/Heterogeneity-Sim/misc/sim_run_grf_15.Rout 2>&1
done