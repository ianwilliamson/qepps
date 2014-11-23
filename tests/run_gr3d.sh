#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DATA_DIR="./comsol_gr3d"
$DIR/verifyBase.sh -lua $DIR/gr3d.lua -st_ksp_type preonly -st_pc_type lu -st_pc_factor_mat_solver_package mumps -st_type sinvert -st_transform
