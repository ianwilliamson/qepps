#!/bin/bash

DATA_DIR="./1E3"
./verifyBase.sh \
	-params ./common/freqs.txt \
	-E0 $DATA_DIR/E.dat \
	-D0 $DATA_DIR/D.dat \
	-K0 $DATA_DIR/K.dat
