#!/bin/bash -l
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=ws
#SBATCH --mem=1G
#SBATCH --array=1-12
#SBATCH --time=10:00:00
#SBATCH --output=/home/scratch/users/stephan.rasp/results/array%a.out
#SBATCH --error=/home/scratch/users/stephan.rasp/results/array%a.err
#SBATCH --export=NONE
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=s.rasp@lmu.de 
#SBATCH --reservation=rebhuhn 
#SBATCH --job-name=1-12vert50mem6to24inc60watdr2
 
python py_loop.py --date $SLURM_ARRAY_TASK_ID --ana vert --tstart 6 --tend 24 --nens 50 --water True --dr 2 --tinc 60
