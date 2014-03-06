cmake_minimum_required(VERSION 2.8)
find_package(VTK REQUIRED)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++0x -fopenmp -D_GLIBCXX_PARALLEL")
#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++0x")
set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")

include_directories(${VTK_INCLUDE_DIRS} . ../../src ../../eigen3)

add_executable (sphdem sphdem.cpp)
TARGET_LINK_LIBRARIES(sphdem ${VTK_LIBRARIES})
