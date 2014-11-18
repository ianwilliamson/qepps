#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DATA_DIR="./ppwg"
$DIR/verifyBase.sh \
	-params $DIR/common/freqs.txt \
	-E0 $DIR/$DATA_DIR/E.dat \
	-D0 $DIR/$DATA_DIR/D.dat \
	-K0 $DIR/$DATA_DIR/K.dat \
	-pep_basis monomial -pep_general -pep_type qarnoldi -st_type sinvert -st_transform \
	"$@"
