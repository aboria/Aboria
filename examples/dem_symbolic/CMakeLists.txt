cmake_minimum_required(VERSION 2.8)
find_package(Boost 1.50.0 REQUIRED)
find_package(VTK)
if (VTK_FOUND)
    add_definitions(-DHAVE_VTK)
endif(VTK_FOUND)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++0x -fopenmp -D_GLIBCXX_PARALLEL")
#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++0x")
#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++0x -ftemplate-backtrace-limit=0")
set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")

include_directories(${VTK_INCLUDE_DIRS} ${Boost_INCLUDE_DIR} ../../src)

add_executable(dem_symbolic dem_symbolic.cpp)
target_link_libraries(dem_symbolic ${VTK_LIBRARIES} ${Boost_LIBRARIES})

