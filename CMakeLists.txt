cmake_minimum_required(VERSION 2.8)
project(Aboria)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cxxtest/build_tools/cmake/"
                      "${CMAKE_SOURCE_DIR}/cmake"
                        ${CMAKE_MODULE_PATH})
find_package(Boost 1.50.0 REQUIRED)

option(Aboria_USE_VTK "Use VTK library" OFF)
if (Aboria_USE_VTK)
    find_package(VTK REQUIRED)
    add_definitions(-DHAVE_VTK)
endif()

option(Aboria_USE_THRUST "Use CUDA Thrust library" OFF)
if (Aboria_USE_THRUST)
    find_package(Thrust REQUIRED)
    add_definitions(-DHAVE_THRUST)
endif()

option(Aboria_USE_EIGEN "Use Eigen Linear Algebra library" OFF)
if (Aboria_USE_EIGEN)
    find_package(Eigen3 REQUIRED)
    add_definitions(-DHAVE_EIGEN)
endif()

#export compiler flags for code completion engines
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 )

set(CXXTEST_ROOT ${CMAKE_SOURCE_DIR})
find_package(CxxTest)

#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++0x -fopenmp -D_GLIBCXX_PARALLEL")
#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++0x ")
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++0x -ftemplate-backtrace-limit=0")
set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")

if (Aboria_COVERAGE)
    set(COVERAGE_FLAGS --coverage)
    #set(COVERAGE_FLAGS -fprofile-arcs -ftest-coverage)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COVERAGE_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${COVERAGE_FLAGS}")
    message("-- Building with coverage reporting")
endif()

include_directories(src)
include_directories(SYSTEM ${CXXTEST_INCLUDES} ${EIGEN3_INCLUDE_DIR} ${VTK_INCLUDE_DIRS} ${Boost_INCLUDE_DIR}  ${PYTHON_INCLUDE_DIRS})

enable_testing()
if (CXXTEST_FOUND)
    add_subdirectory(tests)
endif()
#add_subdirectory(examples)


add_custom_target(generate_doc 
    bjam
    WORKING_DIRECTORY doc
    COMMENT Generating documentation...
    VERBATIM)
