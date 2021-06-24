project ( lua C )
cmake_minimum_required ( VERSION 3.12 )

include_directories ( src ${CMAKE_CURRENT_BINARY_DIR} )

set ( LUA_SRC src/lapi.c src/lcode.c src/lctype.c src/ldebug.c src/ldo.c src/ldump.c src/lfunc.c src/lgc.c src/llex.c
  src/lmem.c src/lobject.c src/lopcodes.c src/lparser.c src/lstate.c src/lstring.c src/ltable.c
  src/ltm.c src/lundump.c src/lvm.c src/lzio.c src/lauxlib.c src/lbaselib.c src/lbitlib.c src/lcorolib.c src/ldblib.c src/liolib.c
  src/lmathlib.c src/loslib.c src/lstrlib.c src/ltablib.c src/lutf8lib.c src/loadlib.c src/linit.c )

add_library(lua STATIC ${LUA_SRC})
set_target_properties(lua PROPERTIES COMPILE_DEFINITIONS "LUA_USE_LINUX")

target_include_directories(lua INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src)
add_library(Lua::lua ALIAS lua)
