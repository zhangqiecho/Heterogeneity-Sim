#!/usr/bin/env bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qi.zhang2@emory.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=10:00:00
#SBATCH --job-name=cluster_sim2
#SBATCH --mem=30GB
#SBATCH --partition=naimi
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

start_time="$(date -u +%s)"
module purge
module load R/4.2.2

for ii in {3..5} 
do 
  Rscript --no-save --no-restore --verbose run_cluster_sim.R $ii 1000  > /home/qzha223/test/output/sim_run_${ii}_1000.Rout 2>&1
done
end_time="$(date -u +%s)"
echo "Job completed!"
echo -e "\n"
elapsed="$(($end_time-$start_time))"
echo "Total of $elapsed seconds elapsed for this job!"