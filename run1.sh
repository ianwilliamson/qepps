#!/bin/bash

ibrun -n 1 -o 0 ./qepps \
	-lambda_tgt 0 \
	-qep_nev 1 \
	-qep_ncv 2 \
	-qep_tol 1e-6 \
	-st_type sinvert \
	-qep_monitor \
	-qep_type linear \
	-qep_linear_cform 1 \
	-params $WORK/data_comsol/freqs.txt \
	-E0 $WORK/data_comsol/E.dat \
	-D0 $WORK/data_comsol/D.dat \
	-K0 $WORK/data_comsol/K.dat

#-qep_linear_explicitmatrix -qep_nev 4 -qep_ncv 20 -qep_tol 1e-9 -qep_terse > ex17_1.tmp 2>&1;
