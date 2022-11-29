#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=2000G
#SBATCH --time=600
#SBATCH -o TG_16_p2.out
#SBATCH --mail-user stevenw@lanl.gov
#SBATCH --mail-type=ALL

srun  ./fierro ../meshes/mesh_taylor_16.geo
