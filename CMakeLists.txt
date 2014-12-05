cmake_minimum_required(VERSION 2.8)
find_package(Boost 1.50.0 COMPONENTS filesystem system date_time python REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(VTK REQUIRED)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cxxtest/build_tools/cmake/" ${CMAKE_MODULE_PATH})
find_package(CxxTest)

#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++0x -fopenmp -D_GLIBCXX_PARALLEL")
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++0x")
set(CMAKE_CPP_FLAGS_DEBUG  "${CMAKE_CPP_FLAGS_DEBUG} -DBOOST_DEBUG_PYTHON" )
set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")

include_directories(${CXXTEST_INCLUDES} ${VTK_INCLUDE_DIRS} ${Boost_INCLUDE_DIR}  ${PYTHON_INCLUDE_DIRS} src)

set(TESTING_DIR tests)
enable_testing()
cxx_test(constructors ${TESTING_DIR}/constructors.h)
target_link_libraries(constructors ${VTK_LIBRARIES} ${Boost_LIBRARIES}) 
cxx_test(particle_container ${TESTING_DIR}/particle_container.h)
target_link_libraries(particle_container ${VTK_LIBRARIES} ${Boost_LIBRARIES}) 
cxx_test(neighbours ${TESTING_DIR}/neighbours.h)
target_link_libraries(neighbours ${VTK_LIBRARIES} ${Boost_LIBRARIES}) 