/*
  Copyright 2019-2021 The University of New Mexico

  This file is part of FIESTA.
  
  FIESTA is free software: you can redistribute it and/or modify it under the
  terms of the GNU Lesser General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option) any
  later version.
  
  FIESTA is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
  details.
  
  You should have received a copy of the GNU Lesser General Public License
  along with FIESTA.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "luaReader.hpp"
#include "input.hpp"
#include "Kokkos_Core.hpp"
#include "debug.hpp"
#include "kokkosTypes.hpp"
#include "string.h"
#include <getopt.h>
#include <iostream>
#include <string>
#include <regex>

using namespace std;

// Construct object and open lua file
luaReader::luaReader(std::string fname){
  L = luaL_newstate();
  luaL_openlibs(L);

  if (luaL_loadfile(L, fname.c_str()) || lua_pcall(L, 0, 0, 0))
    error(L, "Cannot run config file: %s\n", lua_tostring(L, -1));
};

// Int
template<>
void luaReader::get<int>(string key, int &out){
   out = getInt(key);
}

template<>
void luaReader::get<int>(string key, int &out, int val){
  if (undefined(key)) out = val;
  else out = getInt(key);
}

// Double
template<>
void luaReader::get<double>(string key, double &out){
  out = getDouble(key);
}

template<>
void luaReader::get<double>(string key, double &out, double val){
  if (undefined(key)) out = val;
  else out = getDouble(key);
}

// string
template<>
void luaReader::get<string>(string key,string &out){
  out = getString(key);
}

template<>
void luaReader::get<string>(string key, string &out, const char *val){
  if (undefined(key)) out = string(val);
  else out = getString(key);
}

template<>
void luaReader::get<string>(string key, string &out, string val){
  if (undefined(key)) out = val;
  else out = getString(key);
}

// get array of doubles
template<>
void luaReader::getArray(string key, double *out, int n){
  int isnum;
  lua_getglobal(L, key.c_str());
  for (int i=0; i < n; ++i) {
    lua_pushnumber(L, i + 1);
    lua_gettable(L, -2);
    out[i] = (double)lua_tonumberx(L, -1, &isnum);
    lua_pop(L, 1);
  }
  lua_pop(L,1);
}

// Get array of ints
template<>
void luaReader::getArray(string key, int *out, int n){
  int isnum;
  lua_getglobal(L, key.c_str());
  for (int i=0; i < n; ++i) {
    lua_pushnumber(L, i + 1);
    lua_gettable(L, -2);
    out[i] = (double)lua_tointegerx(L, -1, &isnum);
    lua_pop(L, 1);
  }
  lua_pop(L,1);
}

// Get a vector of strings from a lua table
template<>
void luaReader::getArray(string key, vector<string>& out, int n){
  lua_getglobal(L, key.c_str());
  if (lua_istable(L,-1))
      for (int i=0; i<n; ++i) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, -2);
     
        const char *name_c;
     
        if (lua_isstring(L,-1)){
          name_c = lua_tostring(L, -1);
          out.push_back(std::string(name_c));
          lua_pop(L, 1);
        }else{
          error(L, "Error Reading Input File: Could not read value in array '%s', a string was expected.\n", key.c_str());
        }
      }
  else
    error(L, "Error Reading Input File: A string array was expected at '%s'\n", key.c_str());
}

// Call lua function from c (takes integer arguments and returns a double)
double luaReader::call(string f, int n, ...){
  lua_getglobal(L, f.c_str());
  int isnum;

  va_list ap;
  va_start(ap, n);
  double v;
  double z;
  for(int i=0; i<n; ++i) {
    v = va_arg(ap, double);
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

// check if key is defined in lua file
bool luaReader::undefined(string key){
  lua_getglobal(L, key.c_str());
  if (lua_isnoneornil(L,-1)){
    lua_pop(L, 1);
    return true;
  }
  return false;
}

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
// bool luaReader::getBool(std::string key){
int luaReader::getInt(string key){
  int result,isnum;
  lua_getglobal(L, key.c_str());
  if ( lua_isstring(L,-1) ) {
    const char *iresult;
    iresult = lua_tostring(L,-1);

    std::string str(iresult);
    std::regex rgxtrue(R"(\.?(true|on|enable|enabled)\.?)",std::regex_constants::icase);
    std::regex rgxfalse(R"(\.?(false|off|disable|disabled)\.?)",std::regex_constants::icase);

    if ( std::regex_match(str,rgxtrue) ){
      result = 1;
    } else if ( std::regex_match(str,rgxfalse) ) {
      result = 0;
    } else {
      result = (int)lua_tointegerx(L, -1, &isnum);
      if (!isnum)
        error(L, "Error Reading Input File: Could not read value for '%s', an integer or boolean value was expected.\n", key.c_str());
    }
  } else {
    result = lua_toboolean(L,-1);
  }
  lua_pop(L, 1);

  return result;
}

// Lua get double value
double luaReader::getDouble(string key){
  int isnum;
  double result;
  lua_getglobal(L, key.c_str());
  result = (double)lua_tonumberx(L, -1, &isnum);
  if (!isnum)
    error(L, "Error Reading Input File: Could not read value for '%s', an integer was expected.\n", key.c_str());
  lua_pop(L, 1);

  return result;
}

// Lua get string
std::string luaReader::getString(string key){
  const char *result;
  lua_getglobal(L, key.c_str());
  if (lua_isstring(L,-1)){
    result = lua_tostring(L, -1);
  }else{
    error(L, "Error Reading Input File: Could not read value for '%s', a string was expected.\n", key.c_str());
  }
  lua_pop(L, 1);

  return std::string(result);
}
