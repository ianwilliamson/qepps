
PETSC_DIR=$(TACC_PETSC_DIR)
SLEPC_DIR=$(TACC_SLEPC_DIR)
SLEPC_LIB=$(TACC_SLEPC_LIB)
SLEPC_INC=$(TACC_SLEPC_INC)

LUA_DIR=./lua
LUAINC_DIR=$(LUA_DIR)/src
LUALIB_DIR=$(LUA_DIR)/src
LUABIN_DIR=$(LUA_DIR)/src

#CC=icc
CFLAGS=-g -O0 $(LUA_INCS) $(WARN) #-fPIC
WARN=
LUA_INCS=-I$(LUAINC_DIR)
LUA_LIBS=-L$(LUALIB_DIR) -llua -lm

include $(SLEPC_DIR)/conf/slepc_common

LUA_SRC=lcomplex.c config.c
LUA_OBJ=$(LUA_SRC:%.c=%.o)

PETSC_SRC=sweeper.c load.c assemble.c
PETSC_OBJ=$(PETSC_SRC:%.c=%.o)

all: qepps config_reader

qepps: qepps.o $(PETSC_OBJ) $(LUA_OBJ)
	-${CLINKER} qepps.o $(PETSC_OBJ) $(LUA_OBJ) -o ../qepps ${SLEPC_LIB} $(LUA_LIBS)
	${RM} qepps.o
	
config_reader: config_reader.o $(LUA_OBJ)
	-${CLINKER} config_reader.o $(LUA_OBJ) -o ../config_reader ${SLEPC_LIB} $(LUA_LIBS)
	${RM} config_reader.o 

#%.o: %.c
#	${PETSC_COMPILE} $< -o $@ $(LUA_LIBS)

#config:	$(LUA_OBJ) lua
#	${PETSC_COMPILE} -o ../$@ $@.o $(LUA_OBJ) $(LIBS)

#lua: 
#	$(MAKE) -C lua

#MAKESO=$(CC) -shared -Wl,-E
#COMPLEX_NAME=complex
#COMPLEX_LIB=l$(COMPLEX_NAME)
#(COMPLEX_NAME).so: $(COMPLEX_LIB).o
#	$(MAKESO) -o $@ $(COMPLEX_LIB).o
