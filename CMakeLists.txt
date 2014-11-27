cmake_minimum_required(VERSION 2.8)
find_package(Boost 1.50.0 COMPONENTS filesystem system date_time python REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(VTK REQUIRED)
find_package(CxxTest)

#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++0x -fopenmp -D_GLIBCXX_PARALLEL")
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++0x")
set(CMAKE_CPP_FLAGS_DEBUG  "${CMAKE_CPP_FLAGS_DEBUG} -DBOOST_DEBUG_PYTHON" )
set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")

include_directories(${CXXTEST_INCLUDE_DIR} ${VTK_INCLUDE_DIRS} ${Boost_INCLUDE_DIR}  ${PYTHON_INCLUDE_DIRS} src)

set(TESTING_DIR cxxtest)
enable_testing()
CXXTEST_ADD_TEST(unit_tests unit_tests.cpp
                      ${TESTING_DIR}/constructors.h
                      ${TESTING_DIR}/particle_container.h
                      )
target_link_libraries(unit_tests ${PYTHON_LIBRARIES} ${VTK_LIBRARIES} ${Boost_LIBRARIES}) # as needed
