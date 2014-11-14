#!/bin/bash
QEPPS=../qepps
ibrun -n 8 -o 0 $QEPPS -info -lambda_tgt 1 -pep_nev 6 -pep_tol 1e-14 "$@" \
        -pep_basis monomial \
        -pep_general \
        -pep_type qarnoldi \
        -st_type sinvert \
        -st_transform
