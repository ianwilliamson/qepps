#include <slepcpep.h>
#include <petscbag.h>

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>

#include "types.h"
#include "config.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main (int argc, char **argv) 
{
  int i;
  if(argc >= 1)
  {
    lua_State *L = openConfigLUA(argv[1]);
    ParameterSet *parameters = parseConfigParametersLUA(L);
    
    printf("Found %d parameters...\n",parameters->num);
    //for (i=0; i < parameters->num; i++)
    //  printf(" %E \n",parameters->param[i]);
      
    free(parameters);

  }
}
