cmake_minimum_required(VERSION 2.8)
project(Aboria)

set(Aboria_LOG_LEVEL 1 CACHE STRING "Logging level (1 = least, 3 = most)")
add_definitions(-DABORIA_LOG_LEVEL=${Aboria_LOG_LEVEL})

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cxxtest/build_tools/cmake/"
                      "${CMAKE_SOURCE_DIR}/cmake"
                        ${CMAKE_MODULE_PATH})




                    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wreturn-type -Wno-deprecated -std=c++11")
list(APPEND CMAKE_CXX_FLAGS "-Wno-deprecated -std=c++11 -ftemplate-backtrace-limit=0")
#set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} "-Wno-deprecated -std=c++11"")

option(Aboria_USE_THRUST "Use CUDA Thrust library" OFF)
if (Aboria_USE_THRUST)
    set(Aboria_USE_OPENMP ON CACHE BOOL "Use OpenMP for shared memory parallism" FORCE)
    find_package(Thrust REQUIRED)
    find_package(CUDA REQUIRED)

    add_definitions(-D__aboria_have_thrust__)

    option(Aboria_THRUST_TEST_NVCC "Turn on tests using nvcc" OFF)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        option(Aboria_THRUST_TEST_CUDA_CLANG "Turn on tests using cuda clang" OFF)
    endif()
    #http://stackoverflow.com/questions/34996295/trying-to-get-cuda-7-5-to-work-with-gcc-5-x

    add_definitions(-D_MWAITXINTRIN_H_INCLUDED)
    set(
    CUDA_NVCC_FLAGS
    ${CUDA_NVCC_FLAGS};
    --expt-relaxed-constexpr
    #--keep
    #-lineinfo
    --expt-extended-lambda
    -ftemplate-backtrace-limit=0
    #-gencode arch=compute_30,code=sm_30 
    -gencode arch=compute_35,code=sm_35 
    #-gencode arch=compute_50,code=sm_50 
    #-gencode arch=compute_52,code=sm_52 
    #-gencode arch=compute_50,code=compute_50
    )
#   if (CMAKE_BUILD_TYPE STREQUAL Debug)
#       set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-G)
#   endif()
endif()


option(Aboria_USE_OPENMP "Use OpenMP for shared memory parallism" OFF)
if (Aboria_USE_OPENMP)
    find_package(OpenMP REQUIRED)
    add_definitions(-DHAVE_OPENMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


find_package(Boost 1.50.0 REQUIRED serialization)
list(APPEND Aboria_LIBRARIES "${Boost_LIBRARIES}")

option(Aboria_USE_VTK "Use VTK library" OFF)
if (Aboria_USE_VTK)
    find_package(VTK REQUIRED)
    add_definitions(-DHAVE_VTK)
endif()

option(Aboria_USE_GPERFTOOLS "Use Google Profiling tools" OFF)
if (Aboria_USE_GPERFTOOLS)
    find_package(Gperftools REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
    set(GPERFTOOLS_LIBRARIES "-Wl,--no-as-needed ${GPERFTOOLS_LIBRARIES} -Wl,--as-needed") 
    add_definitions(-DHAVE_GPERFTOOLS)
endif()



option(Aboria_USE_EIGEN "Use Eigen Linear Algebra library" OFF)
if (Aboria_USE_EIGEN)
    find_package(Eigen3 REQUIRED)
    option(Aboria_USE_INTEL_MKL "Use Intel MKL for all supported Eigen ops. Remeber to set BLA_VENDOR appropriatly)" OFF)
    add_definitions(-DHAVE_EIGEN)
    if (Aboria_USE_INTEL_MKL)
        set(BLA_VENDOR Intel10_64lp)
        find_package(BLAS REQUIRED)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}")
        list(APPEND Aboria_LIBRARIES "${BLAS_LIBRARIES}")
        add_definitions(-DEIGEN_USE_MKL_ALL)
    endif()
endif()

#export compiler flags for code completion engines
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 )

set(CXXTEST_ROOT ${CMAKE_SOURCE_DIR})
find_package(CxxTest)

option(Aboria_COVERAGE "Compile with coverage info" OFF)
if (Aboria_COVERAGE)
    #set(COVERAGE_FLAGS --coverage)
    set(COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COVERAGE_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${COVERAGE_FLAGS}")
    message("-- Building with coverage reporting")
endif()

include_directories(src)
include_directories(third-party)
include_directories(SYSTEM ${CXXTEST_INCLUDES} ${EIGEN3_INCLUDE_DIR} ${VTK_INCLUDE_DIRS} ${Boost_INCLUDE_DIR}  ${PYTHON_INCLUDE_DIRS} ${THRUST_INCLUDE_DIR})

enable_testing()
if (CXXTEST_FOUND)
    add_subdirectory(tests)
endif()

option(Aboria_BUILD_DOCUMENTATION "Build Aboria documentation" OFF)
if (Aboria_BUILD_DOCUMENTATION)
    add_subdirectory(doc)
endif()


