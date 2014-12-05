#!/bin/bash
#SBATCH -J qepps
#SBATCH -p normal
#SBATCH -N 2 #Number of nodes
#SBATCH -n 16 #Number of mpi tasks
#SBATCH -t 00:30:00
#SBATCH --mail-user=EMAIL@DOMAIN.EDU
#SBATCH --mail-type=end

ibrun ../qepps -lua ./graphene_phc-ppwg_uc.lua -st_pc_factor_mat_solver_package mumps -st_ksp_type preonly -st_pc_type lu 

