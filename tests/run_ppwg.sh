#!/bin/bash

ibrun -n 8 -o 0 ../qepps -lua ./ppwg.lua -st_pc_factor_mat_solver_package mumps -st_ksp_type preonly -st_pc_type lu "$@"
