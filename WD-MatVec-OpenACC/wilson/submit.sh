#!/bin/bash
#SBATCH -N 1
#SBATCH -p GPU-shared
#SBATCH --ntasks-per-node 16
#SBATCH -t 5:00:00
#SBATCH --gres=gpu:p100:1
set -x

#cd /pylon5/phz3a8p/whytet24

#./2D-Wilson-accel
./testmvps
