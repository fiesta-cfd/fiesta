#include "luaReader.hpp"
#include "input.hpp"
#include "Kokkos_Core.hpp"
#include "debug.hpp"
#include "fiesta.hpp"
#include "particle.hpp"
#include "string.h"
#include <getopt.h>
#include <iostream>
#include <string>
#include <regex>

// Construct object and open lua file
luaReader::luaReader(std::string fname){
  L = luaL_newstate();
  luaL_openlibs(L);

  if (luaL_loadfile(L, fname.c_str()) || lua_pcall(L, 0, 0, 0))
    error(L, "Cannot run config file: %s\n", lua_tostring(L, -1));
};

// Get Double Precision Numbers
double luaReader::getDouble(std::string key){
  return getglobdbl(L,key.c_str(),0,0.0);
}

double luaReader::getDouble(std::string key, double val){
  return getglobdbl(L,key.c_str(),1,val);
}

// Get Integer Numbers
int luaReader::getInt(std::string key){
  return getglobint(L,key.c_str(),0,0);
}

int luaReader::getInt(std::string key, int val){
  return getglobint(L,key.c_str(),1,val);
}

// Get Strings
std::string luaReader::getString(std::string key){
  return getglobstr(L,key.c_str(),0,"0");
}

std::string luaReader::getString(std::string key, std::string val){
  return getglobstr(L,key.c_str(),1,val);
}

// Get Boolean Values
bool luaReader::getBool(std::string key){
  return getglobbool(L,key.c_str(),0,0);
}

bool luaReader::getBool(std::string key, bool val){
  return getglobbool(L,key.c_str(),1,1);
}

// Get an array of doubles from a lua table
void luaReader::getDoubles(std::string key, int n, double *out){
  int isnum;
  lua_getglobal(L, key.c_str());
  for (int i=0; i < n; ++i) {
    lua_pushnumber(L, i + 1);
    lua_gettable(L, -2);
    out[i] = (double)lua_tonumberx(L, -1, &isnum);
    lua_pop(L, 1);
  }
}

// Get an array of ints from a lua table
void luaReader::getInts(std::string key, int n, int *out){
  int isnum;
  lua_getglobal(L, key.c_str());
  for (int i=0; i < n; ++i) {
    lua_pushnumber(L, i + 1);
    lua_gettable(L, -2);
    out[i] = (double)lua_tointegerx(L, -1, &isnum);
    lua_pop(L, 1);
  }
}

// Get a vector of strings from a lua table
void luaReader::getStrings(std::string key, int n, std::vector<std::string> &out){
  lua_getglobal(L, key.c_str());
  if (lua_istable(L,-1))
      for (int i=0; i<n; ++i) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, -2);
     
        const char *name_c = (char *)malloc(32 * sizeof(char));
     
        if (lua_isstring(L,-1)){
          name_c = lua_tostring(L, -1);
          out.push_back(std::string(name_c));
          lua_pop(L, 1);
        }else{
          error(L, "Error Reading Input File: Could not read value in array '%s', a string was expected.\n", key.c_str());
        }
        free(&name_c);
      }
  else
    error(L, "Error Reading Input File: A string array was expected at '%s'\n", key.c_str());
}

// Call lua function from c (takes integer arguments and returns a double)
double luaReader::call(std::string f, int n, ...){
  lua_getglobal(L, f.c_str());
  int isnum;

  va_list ap;
  va_start(ap, n);
  int v;
  double z;
  for(int i=0; i<n; ++i) {
    v = va_arg(ap, int);
    lua_pushnumber(L, v);
  }
  va_end(ap);

  if (lua_pcall(L, n, 1, 0) != LUA_OK)
    error(L, "error running function '%s': %s\n",f.c_str(),lua_tostring(L, -1));
  z = (double)lua_tonumberx(L, -1, &isnum);
  if (!isnum)
    error(L, "function '%s' should return a number\n",f.c_str());
  lua_pop(L, 1);

  return z;
}

// Close Script
void luaReader::close(){
  lua_close(L);
}


// Private Methods

// Lua error function
void luaReader::error(lua_State *L, const char *fmt, ...) {
  va_list argp;
  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);
  lua_close(L);
  exit(EXIT_FAILURE);
}

// Lua get boolean value
int luaReader::getglobbool(lua_State *L, const char *var, int defaultable, int default_value) {
  int result;
  int isnil;
  lua_getglobal(L, var);

  if (lua_isnoneornil(L,-1)){
    if (defaultable)
      return default_value;
    else
      error(L,"Error Reading Input File: Could not read value for '%s', a boolean value was expected.\n",var);
  }else{
    if ( lua_isstring(L,-1) ) {
      const char *iresult;
      iresult = lua_tostring(L,-1);

      std::string str(iresult);
      std::regex rgxtrue(R"(\.?(true|on|enable|enabled|1)\.?)",std::regex_constants::icase);
      std::regex rgxfalse(R"(\.?(false|off|disable|disabled|0)\.?)",std::regex_constants::icase);

      if ( std::regex_match(str,rgxtrue) ){
        result = 1;
      } else if ( std::regex_match(str,rgxfalse) ) {
        result = 0;
      } else {
        error(L,"Error Reading Input File: Unknown value for '%s', a boolean value was expected.\n",var);
      }
    } else {
      result = lua_toboolean(L,-1);
    }
  }
  lua_pop(L, 1);
  return result;
}

// Lua get integer value
int luaReader::getglobint(lua_State *L, const char *var, int defaultable, int default_value) {
  int isnum, result;
  lua_getglobal(L, var);
  if (lua_isnoneornil(L, -1)){
    if (defaultable)
      return default_value;
    else
      error(L, "Error Reading Input File: Could not read value for '%s', an integer was expected.\n", var);
  }
  result = (int)lua_tointegerx(L, -1, &isnum);
  if (!isnum)
    error(L, "Error Reading Input File: Could not read value for '%s', an integer was expected.\n", var);
  lua_pop(L, 1);
  return result;
}

// Lua get double value
double luaReader::getglobdbl(lua_State *L, const char *var, int defaultable, double default_value) {
  int isnum;
  double result;
  lua_getglobal(L, var);
  if (lua_isnoneornil(L, -1)){
    if (defaultable)
      return default_value;
    else
      error(L, "Error Reading Input File: Could not read value for '%s', an integer was expected.\n", var);
  }
  result = (double)lua_tonumberx(L, -1, &isnum);
  if (!isnum)
    error(L, "Error Reading Input File: Could not read value for '%s', an integer was expected.\n", var);
  lua_pop(L, 1);
  return result;
}

// Lua get string
std::string luaReader::getglobstr(lua_State *L, const char *var, int defaultable, std::string default_value) {
  // size_t len;
  const char *result;
  result = (char *)malloc(32 * sizeof(char));
  lua_getglobal(L, var);
  if (lua_isnoneornil(L, -1)){
    if (defaultable)
      return default_value;
    else
      error(L, "Error Reading Input File: Could not read value for '%s', a string was expected.\n", var);
  }

  if (lua_isstring(L,-1)){
    result = lua_tostring(L, -1);
  }else{
    error(L, "Error Reading Input File: Could not read value for '%s', a string was expected.\n", var);
  }

  return std::string(result);
}
