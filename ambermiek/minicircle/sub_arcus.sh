#!/bin/bash
STONK DOG

#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --gres=gpu:1 --constraint='gpu_sku:V100'
#SBATCH -p htc

#SBATCH --mail-type=BEGIN,FAIL
#module load /home/orie3911/amber18/

source $AMBERHOME/amber.sh

application="/home/orie3911/amber18/bin/pmemd.cuda.MPI"

JUMBO DOG
