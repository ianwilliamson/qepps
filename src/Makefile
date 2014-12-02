
PETSC_DIR=$(TACC_PETSC_DIR)
SLEPC_DIR=$(TACC_SLEPC_DIR)
SLEPC_LIB=$(TACC_SLEPC_LIB)
SLEPC_INC=$(TACC_SLEPC_INC)

LOCAL_LUA_DIR=./lua
LOCAL_LUA_INC=$(LOCAL_LUA_DIR)/src
LOCAL_LUA_LIB=$(LOCAL_LUA_DIR)/src
LOCAL_LUA_BIN=$(LOCAL_LUA_DIR)/src

GRVY_INC=$(TACC_GRVY_INC)
GRVY_LIB=$(TACC_GRVY_LIB)

WARN=-Wall -Wextra -Wshadow
LUA_INC=-I$(LOCAL_LUA_INC)
LUA_LIB=-L$(LOCAL_LUA_LIB) -llua -lm
CFLAGS=-O2 $(LUA_INC) $(WARN) -I$(GRVY_INC)

include $(SLEPC_DIR)/conf/slepc_common

SRC_FILES=sweeper.c lcomplex.c config.c
OBJ_FILES=$(SRC_FILES:%.c=%.o)

all: qepps

qepps: qepps.o $(OBJ_FILES)
	-${CLINKER} qepps.o $(OBJ_FILES) -o ../qepps ${SLEPC_LIB} $(LUA_LIB) -L$(GRVY_LIB) -lgrvy
	${RM} *.o

