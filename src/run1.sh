#!/bin/bash

ibrun -n 1 -o 0 \
	./qepps -lambda_tgt 55+22i -qep_nev 1 -qep_tol 1e-6 \
	-st_type sinvert \
	-qep_monitor -qep_type linear -qep_linear_cform 1 \
	-qep_st_ksp_type gmres -qep_st_pc_type bjacobi -qep_eps_type krylovschur \
	-params $WORK/data_comsol/freqs.txt \
	-E0 $WORK/data_comsol/E.dat \
	-D0 $WORK/data_comsol/D.dat \
	-K0 $WORK/data_comsol/K.dat


