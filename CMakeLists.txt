cmake_minimum_required(VERSION 2.8)
project(Aboria)
find_package(Boost 1.50.0 REQUIRED)
find_package(VTK)
if (VTK_FOUND)
    add_definitions(-DHAVE_VTK)
endif(VTK_FOUND)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cxxtest/build_tools/cmake/" ${CMAKE_MODULE_PATH})
find_package(CxxTest)

#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++0x -fopenmp -D_GLIBCXX_PARALLEL")
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++0x")
#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++0x -ftemplate-backtrace-limit=0")
set(CMAKE_CPP_FLAGS_DEBUG  "${CMAKE_CPP_FLAGS_DEBUG} -DBOOST_DEBUG_PYTHON" )
set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")

include_directories(${CXXTEST_INCLUDES} ${VTK_INCLUDE_DIRS} ${Boost_INCLUDE_DIR}  ${PYTHON_INCLUDE_DIRS} src)

enable_testing()
add_subdirectory(tests)


