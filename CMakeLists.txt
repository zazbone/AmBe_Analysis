cmake_minimum_required(VERSION 3.4)

project(RootDictExample LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_PLATFORM_INDEPENDENT_CODE ON)

find_package(ROOT 6.20 CONFIG REQUIRED)
# If you want to support <6.20, add this line:
# include("${ROOT_DIR}/modules/RootNewMacros.cmake")
# However, it was moved and included by default in 6.201

root_generate_dictionary(G__AmBe include/AmBe.hpp LINKDEF AmBeLinkDef.h)

#Sources
file(GLOB AmBe_SOURCES src/*.cpp)
file(GLOB AmBe_HEADERS include/*.hpp)

add_library(AmBe SHARED ${AmBe_SOURCES} ${AmBe_HEADERS} G__AmBe.cxx)
target_include_directories(AmBe PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(AmBe PUBLIC ROOT::Core)

## Alternative to add the dictionary to an existing target:
# add_library(DictExample SHARED DictExample.cxx DictExample.h)
# target_include_directories(DictExample PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
# target_link_libraries(DictExample PUBLIC ROOT::Core)
# root_generate_dictionary(G__DictExample DictExample.h MODULE DictExample LINKDEF DictLinkDef.h)