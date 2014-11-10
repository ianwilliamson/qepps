EXE=qepps
CFLAGS=-g -O0
PETSC_DIR=$(TACC_PETSC_DIR)
SLEPC_DIR=$(TACC_SLEPC_DIR)
SLEPC_LIB=$(TACC_SLEPC_LIB)
SLEPC_INC=$(TACC_SLEPC_INC)

include $(SLEPC_DIR)/conf/slepc_common

SRC_FILES = \
	load.c
OBJ_FILES = $(SRC_FILES:%.c=%.o)

$(EXE): $(EXE).o $(OBJ_FILES)
	-${CLINKER} $(EXE).o $(OBJ_FILES) -o ../$(EXE) ${SLEPC_LIB}
	${RM} $(EXE).o $(OBJ_FILES)

%.o: %.c
	echo ${PETSC_COMPILE}
	${PETSC_COMPILE} $< -o $@
