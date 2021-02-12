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

#include <string>
#include "lua.hpp"
#include <vector>

class luaReader {

public:
  luaReader(std::string fname);
  void close();

  template <class T>
  void get(std::string,T&);

  template <class T>
  void get(std::string,T&,T);

  template <class T>
  void get(std::string,T&,const char *);

  template <class T>
  void getArray(std::string,T*,int);
 
  template <class T>
  void getArray(std::string,std::vector<T>&,int);

  double call(std::string, int ,...);

private:
  lua_State *L;

  void error(lua_State *, const char *, ...);

  bool undefined(std::string key);

  bool getBool(std::string);
  int getInt(std::string);
  double getDouble(std::string);
  std::string getString(std::string);

};
