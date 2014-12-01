
PETSC_DIR=$(TACC_PETSC_DIR)
SLEPC_DIR=$(TACC_SLEPC_DIR)
SLEPC_LIB=$(TACC_SLEPC_LIB)
SLEPC_INC=$(TACC_SLEPC_INC)

LUA_DIR=./lua
LUAINC_DIR=$(LUA_DIR)/src
LUALIB_DIR=$(LUA_DIR)/src
LUABIN_DIR=$(LUA_DIR)/src

WARN=-Wall -Wextra -Wshadow
LUA_INCS=-I$(LUAINC_DIR)
LUA_LIBS=-L$(LUALIB_DIR) -llua -lm
CFLAGS=-O2 $(LUA_INCS) $(WARN)

include $(SLEPC_DIR)/conf/slepc_common

SRC_FILES=sweeper.c lcomplex.c config.c
OBJ_FILES=$(SRC_FILES:%.c=%.o)

all: qepps

qepps: qepps.o $(OBJ_FILES)
	-${CLINKER} qepps.o $(OBJ_FILES) -o ../qepps ${SLEPC_LIB} $(LUA_LIBS)
	${RM} *.o

