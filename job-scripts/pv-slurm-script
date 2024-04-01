#!/bin/bash

#SBATCH --job-name="paraview"
#SBATCH --mail-type=NONE         # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=alankardutta@iisc.ac.in    # Where to send mail.  Set this to your email address
#SBATCH -p debug
#SBATCH -t 2-20:00:00
#SBATCH -n 48 
#SBATCH --output=%x-%j.log

echo "Working Directory = $(pwd)"

cd $SLURM_SUBMIT_DIR
export PROG="pvserver --force-offscreen-rendering"
module load paraview 

srun $PROG
