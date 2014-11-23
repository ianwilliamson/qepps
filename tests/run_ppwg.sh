#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DATA_DIR="./ppwg"
$DIR/verifyBase.sh -lua $DIR/ppwg.lua -pep_basis monomial -pep_general -pep_type qarnoldi -st_type sinvert -st_transform "$@"
