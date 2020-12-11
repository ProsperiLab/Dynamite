#!/bin/bash
#SBATCH --job-name=seqgen_sims   #Job name	
#SBATCH --mail-type=NONE   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=brittany.rife@ufl.edu   # Where to send mail	
#SBATCH --time=00:10:00   # Walltime
#SBATCH --output=s_%j.out   # Name output file 
#SBATCH --account=salemi
#SBATCH --qos=salemi
#SBTACH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem=4g

#Record the time and compute node the job ran on
date; hostname; pwd
#Use modules to load the environment for R
module load R/4.0 #phyclust not available in 3.6

Rscript seqgen.R $1
date


