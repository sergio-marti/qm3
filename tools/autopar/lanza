#!/bin/bash
#SBATCH --job-name="qm3.autopar"
#SBATCH --output=out
#SBATCH --error=err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --tasks-per-node=1

cd /storage/projects/vlc78/xexo/prueba
source /storage/projects/vlc78/xexo/qm3_git/rc
python3 autopar.py > log
