#!/bin/bash
#SBATCH -N 1
#SBATCH --time=600
#SBATCH -o TG_3x3x1_p4.out
#SBATCH --mail-user stevenw@lanl.gov
#SBATCH --mail-type=ALL

srun  ./fierro ../meshes/mesh_taylor_3x3x1.geo
