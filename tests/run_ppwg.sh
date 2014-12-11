#!/bin/bash
# This script provides a launcher for the 2D parallel plate waveguide 
# model that can be run under an idev session on stampede
# 
# NOTE: Although this script will correctly resolve the location of the
# qepps exe and the lua configuration script, the lua script points to 
# data files at ./ppwg. This means that one should cd to the directory
# of this script (i.e. tests/) before executing this launcher. At some
# point in the future, referencing of the problem data from within
# LUA may be improved.

# Resolve the directory of this script
# taken from http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

QEPPS_EXE=$DIR/../qepps
CONFIG_LUA=$DIR/ppwg.lua
NUMPROC=8

ibrun -n $NUMPROC -o 0 $QEPPS_EXE -lua $CONFIG_LUA -st_pc_factor_mat_solver_package mumps -st_ksp_type preonly -st_pc_type lu "$@"
