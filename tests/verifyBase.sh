#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
QEPPS=$DIR/../qepps
ibrun -n 16 -o 0 $QEPPS -pep_nev 6 -pep_tol 1e-16 "$@"
