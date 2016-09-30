#!/bin/bash -l
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=ws
#SBATCH --mem=3G
#SBATCH --array=1-12
#SBATCH --time=10:00:00
#SBATCH --output=/home/scratch/users/stephan.rasp/results/array%a.out
#SBATCH --error=/home/scratch/users/stephan.rasp/results/array%a.err
#SBATCH --export=NONE
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=s.rasp@lmu.de 
#SBATCH --job-name=1-12coarse20mem6to24inc30min0height3
 
python py_loop.py --date $SLURM_ARRAY_TASK_ID --ana coarse --tstart 6 --tend 24 --nens 20 --tinc 30 --minmem 0 --height 3000
