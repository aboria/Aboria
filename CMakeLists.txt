# The name of our project is "PDE_BD". CMakeLists files in this project can 
# refer to the root source directory of the project as ${PDE_BD_SOURCE_DIR} and 
# to the root binary directory of the project as ${PDE_BD_BINARY_DIR}. 
project (Aboria) 
cmake_minimum_required (VERSION 2.8) 


FIND_PACKAGE(PythonInterp REQUIRED)
FIND_PACKAGE(PythonLibs REQUIRED)
FIND_PACKAGE(Boost 1.50.0 COMPONENTS python date_time system filesystem timer REQUIRED)
find_package(VTK REQUIRED)


set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++0x")

set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")


set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${Aboria_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${Aboria_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${Aboria_BINARY_DIR}/bin)

# Appends the cmake/modules path inside the MAKE_MODULE_PATH variable which stores the
# directories of additional CMake modules (ie. MacroOutOfSourceBuild.cmake):
set(CMAKE_MODULE_PATH ${Aboria_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH})

# The macro below forces the build directory to be different from source directory:
include(MacroOutOfSourceBuild)

enable_testing()
set(Aboria_INCLUDE_DIRECTORIES  ${Aboria_SOURCE_DIR}/src ${Aboria_SOURCE_DIR}/eigen3 ${Boost_INCLUDE_DIR} ${VTK_INCLUDE_DIRS})
add_subdirectory (src)
#add_subdirectory (python)
add_subdirectory (testing)
