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
#include "block.hpp"
#include "fmt/core.h"
#include "log2.hpp"
#include <typeinfo>

using namespace std;
using fmt::format;

// Construct object and open lua file
luaReader::luaReader(std::string f_, std::string r_):filename(f_),root(r_){
  L = luaL_newstate();
  lua_newtable(L);
  lua_setglobal(L,root.c_str());

  lua_getglobal(L,root.c_str());
  lua_newtable(L);
  lua_setfield(L,-2,"ceq");
  lua_newtable(L);
  lua_setfield(L,-2,"noise");
  lua_newtable(L);
  lua_setfield(L,-2,"restart");
  lua_newtable(L);
  lua_setfield(L,-2,"bc");
  lua_newtable(L);
  lua_setfield(L,-2,"time");
  lua_newtable(L);
  lua_setfield(L,-2,"viscosity");
  lua_newtable(L);
  lua_setfield(L,-2,"buoyancy");
  lua_newtable(L);
  lua_setfield(L,-2,"grid");
  lua_newtable(L);
  lua_setfield(L,-2,"mpi");
  lua_newtable(L);
  lua_setfield(L,-2,"progress");
  lua_newtable(L);
  lua_setfield(L,-2,"status");

  luaL_openlibs(L);

  if (luaL_loadfile(L, filename.c_str()) || lua_pcall(L, 0, 0, 0))
    error(L, "Cannot run config file: %s\n", lua_tostring(L, -1));
};

// get required
template<class T>
void luaReader::get(std::initializer_list<string> keys, T& n){
  bool found=true;
  int top=lua_gettop(L);

  lua_getglobal(L,root.c_str());

  string fullkey=root;
  for (auto key : keys) fullkey=fmt::format("{}.{}",fullkey,key);

  for (auto key : keys){
    lua_getfield(L,-1,key.c_str());
    if(lua_isnoneornil(L,-1)){
      found=false;
      break;
    }
  }
  if (found){
    getValue<T>(n);
    Fiesta::Log::debug("Found {}={}",fullkey,n);
  }else{
    Fiesta::Log::error("Could not find required parameter '{}' in '{}'",fullkey,filename);
    exit(EXIT_FAILURE);
  }
  lua_settop(L,top);
}

// get with default
template<class T>
void luaReader::get(std::initializer_list<string> keys, T& n, T d){
  bool found=true;

  string fullkey=root;
  for (auto key : keys) fullkey=fmt::format("{}.{}",fullkey,key);

  int top=lua_gettop(L);

  lua_getglobal(L,root.c_str());

  for (auto key : keys){
    lua_getfield(L,-1,key.c_str());
    if(lua_isnoneornil(L,-1)){
      found=false;
      break;
    }
  }
  if (found){
    getValue<T>(n);
    Fiesta::Log::debug("Found {}={}",fullkey,n);
  }else{
    Fiesta::Log::debugWarning("Could not find '{}' setting default ({})",fullkey,d);
    n=d;
  }
  lua_settop(L,top);
}

// get array
template<class T>
void luaReader::get(std::initializer_list<string> keys, vector<T>& v, int n){
  bool found=true;
  int top=lua_gettop(L);

  lua_getglobal(L,root.c_str());

  string fullkey=root;
  for (auto key : keys) fullkey=fmt::format("{}.{}",fullkey,key);

  for (auto key : keys){
    lua_getfield(L,-1,key.c_str());
    if(lua_isnoneornil(L,-1)){
      found=false;
      break;
    }
  }
  if (found){
    getArray<T>(v,n);
    Fiesta::Log::debug("Found {}",fullkey);
  }else{
    Fiesta::Log::error("Could not find required parameter '{}' in '{}'",fullkey,filename);
    exit(EXIT_FAILURE);
  }
  lua_settop(L,top);
}

void luaReader::getSpeciesData(struct inputConfig& cf){
  int isnum;

  lua_getglobal(L,root.c_str());
  lua_getfield(L, -1, "species");
  if(lua_isnoneornil(L,-1)){
    Fiesta::Log::error("Could not find {}.species in '{}'",root,filename);
    exit(EXIT_FAILURE);
  }

  if (lua_istable(L,-1)){
    cf.ns = lua_rawlen(L,-1);
    for (int i=0; i<cf.ns; ++i){
      lua_pushnumber(L,i+1);
      lua_gettable(L,-2);
      
      lua_pushstring(L,"name");
      lua_gettable(L,-2);
      cf.speciesName.push_back(lua_tostring(L,-1));
      lua_pop(L,1);

      lua_pushstring(L,"gamma");
      lua_gettable(L,-2);
      cf.gamma.push_back(lua_tonumberx(L,-1,&isnum));
      lua_pop(L,1);

      lua_pushstring(L,"M");
      lua_gettable(L,-2);
      cf.M.push_back(lua_tonumberx(L,-1,&isnum));
      lua_pop(L,1);

      lua_pushstring(L,"mu");
      lua_gettable(L,-2);
      cf.mu.push_back(lua_tonumberx(L,-1,&isnum));
      lua_pop(L,1);

      lua_pop(L,1);
    }
  }else{
      error(L, "Error Reading Input File: Could not read blocks.\n");
  }
  lua_pop(L,1);
  lua_pop(L,1);
}

void luaReader::getIOBlock(struct inputConfig& cf, rk_func* f, int ndim, vector<blockWriter<float> >& blocks){
  int isnum;
  size_t numElems;
  size_t numBlocks;
  size_t frq,avg;

  lua_getglobal(L,root.c_str());
  lua_getfield(L, -1, "ioviews");
  if(lua_isnoneornil(L,-1)){
    Fiesta::Log::error("Could not find {}.ioviews in '{}'",root,filename);
    exit(EXIT_FAILURE);
  }

  if (lua_istable(L,-1)){
    numBlocks=lua_rawlen(L,-1);
    for (int i=0; i<numBlocks; ++i){
      string myname,mypath;
      vector<size_t> start,limit,stride;
      lua_pushnumber(L,i+1);
      lua_gettable(L,-2);
      
      lua_getfield(L,-1,"name");
      myname.assign(lua_tostring(L,-1));
      lua_pop(L,1);

      lua_getfield(L,-1,"path");
      if (!lua_isnoneornil(L,-1))
        mypath.assign(lua_tostring(L,-1));
      else
        mypath.assign("./");
      lua_pop(L,1);

      lua_getfield(L,-1,"frequency");
      frq = lua_tointegerx(L,-1,&isnum);
      lua_pop(L,1);

      lua_getfield(L,-1,"average");
      if (!lua_isnoneornil(L,-1))
        avg = lua_tointegerx(L,-1,&isnum);
      else
        avg = 1;
      lua_pop(L,1);

      lua_getfield(L,-1,"start");
      if (lua_istable(L,-1)){
        numElems = lua_rawlen(L,-1);
        for(int j=0;j<numElems; ++j){
          lua_pushnumber(L,j+1);
          lua_gettable(L,-2);
          start.push_back((size_t)lua_tonumberx(L, -1, &isnum)-1);
          lua_pop(L,1);
        }
      }else{
        for (int j=0; j<ndim; ++j)
          start.push_back(0);
      }
      lua_pop(L,1);

      lua_getfield(L,-1,"limit");
      if (lua_istable(L,-1)){
        numElems = lua_rawlen(L,-1);
        for(int j=0;j<numElems; ++j){
          lua_pushnumber(L,j+1);
          lua_gettable(L,-2);
          limit.push_back((size_t)lua_tonumberx(L, -1, &isnum)-1);
          lua_pop(L,1);
        }
      }else{
          limit.push_back(cf.glbl_nci-1);
          limit.push_back(cf.glbl_ncj-1);
          if (ndim > 2) limit.push_back(cf.glbl_nck-1);
      }
      lua_pop(L,1);

      lua_getfield(L,-1,"stride");
      if (lua_istable(L,-1)){
        numElems = lua_rawlen(L,-1);
        for(int j=0;j<numElems; ++j){
          lua_pushnumber(L,j+1);
          lua_gettable(L,-2);
          stride.push_back((size_t)lua_tonumberx(L, -1, &isnum));
          lua_pop(L,1);
        }
      }else{
        for (int j=0; j<ndim; ++j)
          stride.push_back(1);
      }
      lua_pop(L,1);

      blocks.push_back(blockWriter<float>(cf,f,myname,mypath,avg,frq,start,limit,stride));
      lua_pop(L,1);
    }
  }else{
      error(L, "Error Reading Input File: Could not read blocks.\n");
  }
}

// Call lua function from c (takes integer arguments and returns a double)
double luaReader::call(string f, int n, ...){
  lua_getglobal(L,root.c_str());
  lua_getfield(L, -1, f.c_str());
  if(lua_isnoneornil(L,-1)){
    Fiesta::Log::error("Could not find {}.{} in '{}'",root,f,filename);
    exit(EXIT_FAILURE);
  }
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
  lua_pop(L, 1);

  return z;
}

// Close Script
void luaReader::close(){
  lua_close(L);
}

// Private Methods
template<>
void luaReader::getValue(double& n){
  int isnum;
  n=(double)lua_tonumberx(L,-1,&isnum);
}

template<>
void luaReader::getValue(int& n){
  int isnum;
  n=(int)lua_tointegerx(L,-1,&isnum);
}
 
template<>
void luaReader::getValue(bool& n){
  n=(bool)lua_toboolean(L,-1);
}
 
template<>
void luaReader::getValue(string& n){
  if (lua_isstring(L,-1)){
    n = lua_tostring(L, -1);
  }
}

template<>
void luaReader::getArray(vector<double>& out, int n){
  int isnum;
  //lua_getglobal(L, key.c_str());
  for (int i=0; i < n; ++i) {
    lua_pushnumber(L, i + 1);
    lua_gettable(L, -2);
    out.push_back((double)lua_tonumberx(L, -1, &isnum));
    lua_pop(L, 1);
  }
  //lua_pop(L,1);
}

template<>
void luaReader::getArray(vector<int>& out, int n){
  int isnum;
  //lua_getglobal(L, key.c_str());
  for (int i=0; i < n; ++i) {
    lua_pushnumber(L, i + 1);
    lua_gettable(L, -2);
    out.push_back((size_t)lua_tointegerx(L, -1, &isnum));
    lua_pop(L, 1);
  }
  //lua_pop(L,1);
}

template<>
void luaReader::getArray(vector<size_t>& out, int n){
  int isnum;
  //lua_getglobal(L, key.c_str());
  for (int i=0; i < n; ++i) {
    lua_pushnumber(L, i + 1);
    lua_gettable(L, -2);
    out.push_back((size_t)lua_tointegerx(L, -1, &isnum));
    lua_pop(L, 1);
  }
  //lua_pop(L,1);
}

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

// explicit instantiation
template void luaReader::get<string>(std::initializer_list<string>,string&);
template void luaReader::get<int>(std::initializer_list<string>,int&);
template void luaReader::get<double>(std::initializer_list<string>,double&);
template void luaReader::get<bool>(std::initializer_list<string>,bool&);
// explicit instantiation
template void luaReader::get<string>(std::initializer_list<string>,string&,string);
template void luaReader::get<int>(std::initializer_list<string>,int&,int);
template void luaReader::get<double>(std::initializer_list<string>,double&,double);
template void luaReader::get<bool>(std::initializer_list<string>,bool&,bool);
// explicit instantiation
template void luaReader::get<int>(std::initializer_list<string>,vector<int>&,int);
template void luaReader::get<double>(std::initializer_list<string>,vector<double>&,int);
template void luaReader::get<size_t>(std::initializer_list<string>,vector<size_t>&,int);
