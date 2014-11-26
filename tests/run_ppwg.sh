#!/bin/bash

MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

LUA_CFG=$MYDIR/ppwg.lua
QEPPS_BIN=$MYDIR/../qepps

ibrun -n 8 -o 0 $QEPPS_BIN -lua $LUA_CFG -pep_nev 6 -pep_tol 1e-16 \
 -st_ksp_type preonly -st_pc_type lu -st_pc_factor_mat_solver_package mumps -st_type sinvert -st_transform
 #-pep_basis monomial -pep_general -pep_type qarnoldi -st_type sinvert -st_transform "$@"
