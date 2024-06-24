#!/bin/bash

###############################################################################
#
#SBATCH --time=15:00:00                  # Job run time (hh:mm:ss)
#SBATCH --nodes=2                        # Number of nodes
#SBATCH --mem-per-cpu=4550
#SBATCH --job-name=jobarray_job          # Name of batch job
#SBATCH --account=qgp
#SBATCH --partition=qgp            # Partition (queue)
#SBATCH --array=1-101                      # Job array indicies
##SBATCH --mail-user=NetID@illinois.edu  # Send email notifications
##SBATCH --mail-type=BEGIN,END           # Type of email notifications to send
#
###############################################################################


cd /projects/jnorhos/tdore/Hydro_Bjorken/finiteMuB/FreezeOutScan/Code/YukariHydro
lineToLook=$SLURM_ARRAY_TASK_ID
rho=$(awk -v line=$lineToLook 'NR==line' $rhoList)
echo $rho
echo $lineToLook
pwd
echo "python3 Hydro_Yukari_DNMR_CSB.py $eos $rho 1.5 critical"
python3 Hydro_Yukari_DNMR_CSB.py $eos $rho 1.5 critical

