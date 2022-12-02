#!/bin/bash
#SBATCH -N 1
#SBATCH --time=600
#SBATCH -o TG_16_p2.out
#SBATCH --mail-user stevenw@lanl.gov
#SBATCH --mail-type=ALL

srun  ./fierro ../meshes/mesh_taylor_16.geo
