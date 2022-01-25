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

#ifndef LUAREADER_H
#define LUAREADER_H

#include <string>
#include "lua.hpp"
#include <vector>
#include "block.hpp"
#include "input.hpp"
#include "rkfunction.hpp"

class luaReader {

public:
  luaReader(std::string fname,std::string root);
  void close();

  template <class T>
  void getArray(std::string,T*,int);
 
  template <class T>
  void getArray(std::vector<T>&,int);

  void getSpeciesData(struct inputConfig&);
  void getIOBlock(struct inputConfig&, std::unique_ptr<class rk_func>&, int, vector<blockWriter<float>>&);

  template <class T>
  void get(std::initializer_list<string> keys, T& n);

  template <class T, class S>
  void get(std::initializer_list<string> keys, T& n, S d);

  template <class T>
  void get(std::initializer_list<string> keys, std::vector<T>& v, int n);

  template <class T>
  void getValue(T& n);

  FSCAL call(std::string, int ,...);

private:
  lua_State *L;

  string root;
  string filename;

  void error(lua_State *, const char *, ...);

  bool undefined(std::string key);

  bool getBool(std::string);
  int getInt(std::string);
  FSCAL getDouble(std::string);
  std::string getString(std::string);

};

#endif
