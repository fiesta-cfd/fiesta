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
