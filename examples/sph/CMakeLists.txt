cmake_minimum_required(VERSION 2.8)
find_package(VTK REQUIRED)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++0x")
set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")

include_directories(${VTK_INCLUDE_DIRS} . ../../src ../../eigen3)

add_executable (sph sph.cpp)
TARGET_LINK_LIBRARIES(sph ${VTK_LIBRARIES})
