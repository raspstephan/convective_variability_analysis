#!/bin/bash -l
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=ws
#SBATCH --mem=8G
#SBATCH --array=1-12
#SBATCH --time=10:00:00
#SBATCH --output=/home/scratch/users/stephan.rasp/results/array%a.out
#SBATCH --error=/home/scratch/users/stephan.rasp/results/array%a.err
#SBATCH --export=NONE
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=s.rasp@lmu.de 
#SBATCH --job-name=1-12vert50mem6to24inc60min0height3
 
python py_loop.py --date $SLURM_ARRAY_TASK_ID --ana vert --tstart 6 --tend 24 --nens 50 --tinc 60 --height 3000
