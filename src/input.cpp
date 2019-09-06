#include "input.hpp"
#include "Kokkos_Core.hpp"
#include "lsdebug.hpp"

// Lua error function
void error(lua_State *L, const char *fmt, ...){
    va_list argp;
    va_start(argp, fmt);
    vfprintf(stderr, fmt, argp);
    va_end(argp);
    lua_close(L);
    exit(EXIT_FAILURE);
}

//Lua get boolean value
int getglobbool(lua_State *L, const char *var) {
    int isnum, result;
    lua_getglobal(L, var);
    result = (int)lua_toboolean(L, -1);
    lua_pop(L,1);
    return result;
}

//Lua get integer value
int getglobint(lua_State *L, const char *var) {
    int isnum, result;
    lua_getglobal(L, var);
    result = (int)lua_tointegerx(L, -1, &isnum);
    if (!isnum)
        error(L, "'%s' should be an integer\n", var);
    lua_pop(L,1);
    return result;
}

//Lua get double value
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

struct inputConfig executeConfiguration(char * fname){

    struct inputConfig myConfig;
    /* Create Lua State */
    lua_State *L = luaL_newstate(); //Opens Lua
    luaL_openlibs(L);               //opens the standard libraries

    /* Open and run Lua configuration file */
    if (luaL_loadfile(L,fname) || lua_pcall(L,0,0,0))
        error(L, "Cannot run config file: %s\n", lua_tostring(L, -1));

    /* get global variables from Lua results */
    myConfig.glbl_nci    = getglobint (L, "ni" );
    myConfig.glbl_ncj    = getglobint (L, "nj" );
    myConfig.glbl_nck    = getglobint (L, "nk" );
    myConfig.nv          = getglobint (L, "nv" );
    myConfig.dx          = getglobdbl (L, "dx" );
    myConfig.dy          = getglobdbl (L, "dy" );
    myConfig.dz          = getglobdbl (L, "dz" );
    myConfig.gamma       = getglobdbl (L, "gamma" );
    myConfig.xProcs      = getglobint (L, "procsx");
    myConfig.yProcs      = getglobint (L, "procsy");
    myConfig.zProcs      = getglobint (L, "procsz");

    snprintf(myConfig.inputFname,32,"%s",fname);

    int isnum;
    double z;

    lua_getglobal(L,"f");
    lua_pushnumber(L,2.0);
    lua_pushnumber(L,3.0);
    if (lua_pcall(L,2,1,0) != LUA_OK)
        error(L, "error running function 'f': %s\n",lua_tostring(L, -1));
    z = lua_tonumberx(L,-1,&isnum);
    if (!isnum)
        error(L, "function 'f' should return a number");
    lua_pop(L,1);
    printf("Function is %f\n",z);

    /* Done with Lua */
    lua_close(L);

    /* calculate number of nodes from number of cells */
    myConfig.glbl_ni = myConfig.glbl_nci + 1;
    myConfig.glbl_nj = myConfig.glbl_ncj + 1;
    myConfig.glbl_nk = myConfig.glbl_nck + 1;

    return myConfig;
}

int loadInitialConditions(struct inputConfig cf, const Kokkos::View<double****> v){
    
    lua_State *L = luaL_newstate(); //Opens Lua
    luaL_openlibs(L);               //opens the standard libraries

    /* Open and run Lua configuration file */
    if (luaL_loadfile(L,cf.inputFname) || lua_pcall(L,0,0,0))
        error(L, "Cannot run config file: %s\n", lua_tostring(L, -1));
    
    typename Kokkos::View<double****>::HostMirror hostV = Kokkos::create_mirror_view(v);
    for (int k=0; k<cf.nck; ++k){
        for (int j=0; j<cf.ncj; ++j){
            for (int i=0; i<cf.nci; ++i){
                hostV(i,j,k,0) = 1.25;
            }
        }
    }

    Kokkos::deep_copy(v,hostV);
    lua_close(L);

    return 0;
}
