
#define Complex	double complex

static Complex Pget(lua_State *L, int i);

static int pushcomplex(lua_State *L, Complex z);

static int Leq(lua_State *L);

static int Ltostring(lua_State *L);

LUALIB_API int luaopen_complex(lua_State *L);


