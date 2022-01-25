#include <string>
#include <vector>
#include <filesystem>
#include <regex>
#include <fmt/core.h>
#include <getopt.h>
#include <unistd.h>

#include "lua.hpp"
#include "log2.hpp"

void error(lua_State *L, const char *fmt, ...) {
  va_list argp;
  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);
  lua_close(L);
  exit(EXIT_FAILURE);
}

int main(int argc, char **argv){
  int verbosity = 4;
  bool colorFlag = false;
  static struct option long_options[] = {
      {"verbosity", optional_argument, NULL, 'v'},
      {"color", optional_argument, NULL, 'c'},
      {NULL, 0, NULL, 0}};

  std::string copt;
  int c = 1;
  int opt_index;

  while ((c = getopt_long(argc, argv, "Vv:ctn:", long_options, &opt_index)) != -1) {
    switch (c) {
      case 'c':
        if (optarg)
          copt = std::string(optarg);
        else
          copt = "on";
        break;
      case 'v':
        if (optarg)
          verbosity = atoi(optarg);
        else
          verbosity = 3;
        break;
    }
  }
  colorFlag = false;
  if (copt.compare("on") == 0)
    colorFlag = true;
  else if (copt.compare("auto") == 0 && isatty(fileno(stdout)))
      colorFlag = true;

  Log::Logger(verbosity,colorFlag,0);

  std::string filename{"fiesta.lua"};
  std::string root{"fiesta"};

  lua_State *L;
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
  lua_setfield(L,-2,"hdf5");
  lua_newtable(L);
  lua_setfield(L,-2,"eos");
  lua_newtable(L);
  lua_setfield(L,-2,"progress");
  lua_newtable(L);
  lua_setfield(L,-2,"status");

  luaL_openlibs(L);

  if (luaL_loadfile(L, filename.c_str()) || lua_pcall(L, 0, 0, 0))
    error(L, "Cannot run config file: %s\n", lua_tostring(L, -1));
  int isnum;
  size_t numElems;
  size_t numBlocks;
  size_t frq,avg;
  bool defaultable=false;
  int dcount=0;

  lua_getglobal(L,root.c_str());
  lua_getfield(L, -1, "ioviews");
  if(lua_isnoneornil(L,-1)){
    Log::error("Could not find {}.ioviews in '{}'",root,filename);
    exit(EXIT_FAILURE);
  }

  std::string rpath;
  int top=lua_gettop(L);
  lua_getglobal(L,root.c_str());
  lua_getfield(L,-1,"restart");
  lua_getfield(L,-1,"path");
  if (lua_isstring(L,-1)){
    rpath = lua_tostring(L, -1);
  }else{
    rpath = "./";
  }

  if (std::filesystem::exists(rpath)){
    std::string rpattern = fmt::format(".*/restart-(\\d+|auto)\\.(h5|xmf)");
    dcount = 0;
    for(auto const& de : std::filesystem::directory_iterator{rpath}){
      if (std::regex_match(de.path().c_str(),std::regex(rpattern))){
        Log::warning("Deleting '{}'.",de.path().c_str());
        std::filesystem::remove(de);
        dcount+=1;
      }
    }
    Log::message("Deleted {} restart files.",dcount);
    if (std::filesystem::is_empty(rpath)){
      std::filesystem::remove(rpath);
      Log::warning("Removing empty restart directory '{}'.",rpath);
    }else{
      Log::error("Restart directory '{}' is not empty.  Skipping.",rpath);
    }
  }

  lua_settop(L,top);
  Log::debug("Found Restart Path: '{}'",rpath);

  if (lua_istable(L,-1)){
    numBlocks=lua_rawlen(L,-1);
    for (size_t i=0; i<numBlocks; ++i){
      std::string myname,mypath;
      std::vector<size_t> start,limit,stride;
      lua_pushnumber(L,i+1);
      lua_gettable(L,-2);
      
      lua_getfield(L,-1,"name");
      myname.assign(lua_tostring(L,-1));
      lua_pop(L,1);

      lua_getfield(L,-1,"path");
      if (!lua_isnoneornil(L,-1)){
        mypath.assign(lua_tostring(L,-1));
      }else{
        mypath.assign("./");
      }

      std::string mypattern = fmt::format(".*/{}-\\d+\\.(h5|xmf)",myname);

      dcount = 0;
      if (std::filesystem::exists(mypath)){
        for(auto const& de : std::filesystem::directory_iterator{mypath}){
          if (std::regex_match(de.path().c_str(),std::regex(mypattern))){
            Log::warning("Deleting '{}'.",de.path().c_str());
            std::filesystem::remove(de);
            dcount += 1;
          }
        }
        Log::message("Deleted {} '{}' files.",dcount,myname);
        if (std::filesystem::is_empty(mypath)){
          std::filesystem::remove(mypath);
          Log::warning("Removing empty '{}' directory '{}'.",myname,mypath);
        }else{
          Log::error("'{}' directory '{}' is not empty.  Skipping.",myname,mypath);
        }
      }

      lua_pop(L,1);

      lua_pop(L,1);
    }
  }else{
      error(L, "Error Reading Input File: Could not read blocks.\n");
  }
}
