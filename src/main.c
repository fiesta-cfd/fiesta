#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"

void error(lua_State *L, const char *fmt, ...){
    va_list argp;
    va_start(argp, fmt);
    vfprintf(stderr, fmt, argp);
    va_end(argp);
    lua_close(L);
    exit(EXIT_FAILURE);
}

int getglobbool(lua_State *L, const char *var) {
    int isnum, result;
    lua_getglobal(L, var);
    result = (int)lua_toboolean(L, -1);
    lua_pop(L,1);
    return result;
}

int getglobint(lua_State *L, const char *var) {
    int isnum, result;
    lua_getglobal(L, var);
    result = (int)lua_tointegerx(L, -1, &isnum);
    if (!isnum)
        error(L, "'%s' should be an integer\n", var);
    lua_pop(L,1);
    return result;
}

double getglobdbl(lua_State *L, const char *var) {
    int isnum;
    double result;
    lua_getglobal(L, var);
    result = (double)lua_tonumberx(L, -1, &isnum);
    if (!isnum)
        error(L, "'%s' should be a double\n", var);
    lua_pop(L,1);
    return result;
}

int main(){
    char buff[256];
    int w,h;
    double gamma;
    int ceq;
    lua_State *L = luaL_newstate(); //Opens Lua
    //luaL_openlibs(L);               //opens the standard libraries

    if (luaL_loadfile(L,"input.lua") || lua_pcall(L,0,0,0))
        error(L, "Cannot run config file: %s\n", lua_tostring(L, -1));
    w     = getglobint (L, "width" );
    h     = getglobint (L, "height");
    gamma = getglobdbl (L, "gamma" );
    ceq   = getglobbool(L, "ceq"   );

    printf("Width = %d, Height = %d, Gamma = %4.2f, CEQ = %d\n", w, h, gamma,ceq);

    lua_close(L);
    return 0;
}
