#include <string>
#include "lua.hpp"

class luaReader {

public:
  luaReader(std::string fname);
  void close();

  double getDouble(std::string,double);
  double getDouble(std::string);
  int getInt(std::string,int);
  int getInt(std::string);
  std::string getString(std::string);
  std::string getString(std::string,std::string);
  bool getBool(std::string);
  bool getBool(std::string,bool);

  void getDoubles(std::string,int,double*);
  void getInts(std::string,int,int*);
  void getStrings(std::string,int,std::vector<std::string>&);

  double call(std::string, int ,...);

private:
  lua_State *L;

  void error(lua_State *, const char *, ...);
  int getglobbool(lua_State *, const char *, int, int);
  int getglobint(lua_State *, const char *, int, int);
  double getglobdbl(lua_State *, const char *, int, double);
  std::string getglobstr(lua_State *, const char *, int, std::string);
};
