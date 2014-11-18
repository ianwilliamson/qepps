#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DATA_DIR="./comsol_gr3d"
$DIR/verifyBase.sh \
	-params $DIR/common/freqs.txt \
	-E0 $DIR/$DATA_DIR/E.dat \
	-D0 $DIR/$DATA_DIR/D.dat \
	-K0 $DIR/$DATA_DIR/K.dat \
	"$@" \
	-st_ksp_type preonly -st_pc_type lu -st_pc_factor_mat_solver_package mumps
