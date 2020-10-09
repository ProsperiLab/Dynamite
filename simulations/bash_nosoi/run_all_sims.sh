#!/bin/bash
#SBATCH --job-name=nosoi_sims   #Job name	
#SBATCH --mail-type=NONE   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=brittany.rife@ufl.edu   # Where to send mail	
#SBATCH --time=1:00:00   # Walltime
#SBATCH --output=r_sim.%j.out   # Name output file 
#SBATCH --account=salemi
#SBATCH --qos=salemi
#SBTACH --ntasks=1
#SBATCH --cpus-per-task=1 ## Needs to be N+1
#SBATCH --mem=1g

# Start 1000 jobs

for i in $(seq 1000);

do

  sbatch nosoi_sim.sh $i

done

