CFLAGS=
PETSC_DIR=$(TACC_PETSC_DIR)
SLEPC_DIR=$(TACC_SLEPC_DIR)
SLEPC_LIB=$(TACC_SLEPC_LIB)
SLEPC_INC=$(TACC_SLEPC_INC)

include $(SLEPC_DIR)/conf/slepc_common

qepps: qepps.o chkopts
	-${CLINKER} qepps.o -o qepps ${SLEPC_LIB}
	${RM} qepps.o

