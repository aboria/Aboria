/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Aboria.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef DOC_GETTING_STARTED_H_
#define DOC_GETTING_STARTED_H_

#include <cxxtest/TestSuite.h>

//[getting_started
/*`


[section Installation and Getting Started] [section Getting prerequisites]

This software is tested on Ubuntu 14.04LTS with the GCC compiler (version
5.4.1), and Clang compiler (version 3.8.0). See our
[@https://travis-ci.org/martinjrobins/Aboria Travis CI page] for more details.

You will need to install a C++ compiler with C++14 support, and the
[@https://cmake.org/ CMake] build software prior to using Aboria, which you can
do on a Debian-based OS using

``
sudo apt-get install build-essential cmake
``

The only required dependency is the [@http://www.boost.org Boost] library.
Optional dependencies are [@http://www.vtk.org/ The Visualization Toolkit],
[@http://eigen.tuxfamily.org Eigen] (version >= 3.3~beta1),
[@http://www.h2lib.org H2Lib], [@http://thrust.github.io Trust] and
[@http://www.openmp.org OpenMP], all of which add extra functionality. To
install all these dependencies in a Debian-based OS you can type

``
sudo apt-get install libboost-dev libvtk5-dev libeigen3-dev libthrust-dev
``
[note replace `libvtk5-dev` with `libvtk6-dev` if necessary]

[note If you wish to use H2Lib you will need to download and compile the source
manually]

[endsect]

[section Compiling a program using Aboria]

Aboria is a header-only library, so at a minimum you will need to add the
`Aboria/src` directory to your include path and include the `Aboria.h` header,
e.g.

``
#include <Aboria.h>
``

If you wish to use any of the optional dependencies you will need to install,
include and/or link the required library as normal, and define one or more of
the following compiler definitions to "turn on" this functionality within Aboria

[table Optional dependencies compiler definitions [[Library Name] [Compiler
definition]] [[VTK] [HAVE_VTK]] [[Eigen] [HAVE_EIGEN]] [[H2Lib] [HAVE_H2LIB]]
[[Thrust] [HAVE_THRUST]]
]

If you are familiar with compiling C++ projects, this might be all the
information you need to incorporate Aboria into your own build system. If you
wish for more details, the following provides a step-by-step guide on compiling
your first Aboria program, using the popular CMake build system.

First clone the Aboria GitHub repository like so

[teletype]
```
$ git clone https://github.com/martinjrobins/Aboria
```
[c++]

Then copy and paste the code below into a C++ source file named
`getting_started.cpp`.

*/

#include "Aboria.h"

using namespace Aboria;

//<-
class DocGettingStartedTest : public CxxTest::TestSuite {
public:
  ABORIA_VARIABLE(cu_velocity, vdouble2, "velocity")

  void test_getting_started(void) {
    //->
    //=int main() {
    /*
     * Create a 2d particle container type with one
     * additional variable "velocity", represented
     * by a 2d double vector
     */
    ABORIA_VARIABLE(velocity, vdouble2, "velocity")
    typedef Particles<std::tuple<velocity>, 2> container_t;
    typedef typename container_t::position position;

    /*
     * create a particle set with size N
     */
    const int N = 100;
    container_t particles(N);

    std::uniform_real_distribution<double> uni(0, 1);
    std::default_random_engine gen;
    for (int i = 0; i < N; ++i) {
      /*
       * set a random position, and initialise velocity
       */
      get<position>(particles)[i] = vdouble2(uni(gen), uni(gen));
      get<velocity>(particles)[i] = vdouble2(0, 0);
    }

    /*
     * write particle container to a vtk
     * unstructured grid file
     */
//<-
#ifdef HAVE_VTK
    //->
    vtkWriteGrid("aboria", 0, particles.get_grid(true));
//<-
#endif
    //->

    //=}
    /*`

    Now copy and paste the CMake config file below into another file called
    `CMakeLists.txt`. Note that this assumes that the Aboria repository is in
    your current directory, which it will be if you just cloned the repo using
    the git command above.


    [teletype]
    ```
    cmake_minimum_required(VERSION 2.8)

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

    # Boost
    find_package(Boost 1.50.0 REQUIRED serialization)
    list(APPEND LIBRARIES ${Boost_LIBRARIES})
    list(APPEND INCLUDES ${Boost_INCLUDE_DIRS})

    # VTK
    find_package(VTK REQUIRED)
    if (VTK_FOUND)
        add_definitions(-DHAVE_VTK)
    endif(VTK_FOUND)
    list(APPEND LIBRARIES ${VTK_LIBRARIES})
    list(APPEND INCLUDES ${VTK_INCLUDE_DIRS})

    # Aboria
    set(Aboria_LOG_LEVEL 1 CACHE STRING "Logging level (1 = least, 3 = most)")
    add_definitions(-DABORIA_LOG_LEVEL=${Aboria_LOG_LEVEL})
    list(APPEND INCLUDES Aboria/src)
    list(APPEND INCLUDES Aboria/third-party)

    include_directories(src ${INCLUDES})

    set(SOURCE
        getting_started.cpp
        )

    add_executable(getting_started ${SOURCE})
    target_link_libraries(getting_started ${LIBRARIES})
    ```

    If you wish to use [@http://www.vtk.org/ The Visualisation Toolkit] or
    [@http://eigen.tuxfamily.org Eigen] with Aboria then you will need to define
    the `HAVE_VTK` and `HAVE_EIGEN` compiler definitions, as done using the
    above `CMakeLists.txt`


    Finally, configure and compile `getting_started.cpp`, then run it like so

    ```
    $ cmake .
    $ make
    $ ./getting_started
    ```
    [c++]

    You now should have a file `aboria00000.vtu` in your directory, which you
    can open and view using any application that supports the VTK unstructured
    grid format (such as Paraview)

    [endsect]

    [section Compiling with Eigen]

    Aboria interfaces with the [@http://eigen.tuxfamily.org Eigen] library to
    provide matrix-free linear operators, linear solvers and preconditioners.
    See [link aboria.evaluating_and_solving_kernel_op] for more details.

    Assuming you are using CMake as per the instructions above, you can include
    Eigen in the build process by inserting the following into your
    `CMakeLists.txt`

    [teletype]
    ```
    # optional if you already have a `FindEigen3.cmake` on your system
    set(CMAKE_MODULE_PATH  "${CMAKE_SOURCE_DIR}/Aboria/cmake"
                           ${CMAKE_MODULE_PATH})
    # end optional

    find_package(Eigen3 REQUIRED)
    list(APPEND INCLUDES ${EIGEN3_INCLUDE_DIR})
    add_definitions(-DHAVE_EIGEN)
    ```

    The first line sets the `CMAKE_MODULE_PATH` to search in the Aboria `cmake`
    directory for additional `.cmake` files. This is only necessary if you don't
    already have a `FindEigen3.cmake` module or config file on your system and
    wish to use the one bundled with Aboria.

    The most important line for using the Eigen functionality within Aboria is
    the `add_definitions` instruction, which adds the `HAVE_EIGEN` compiler
    definition. This definition is used to "turn on" all Eigen functionality
    within Aboria.

    [endsect] [section Compiling with VTK]

    Aboria uses [@http://www.vtk.org/ The Visualisation Toolkit] to allow it to
    write out or read in particle data as a `vtkUnstructuredGrid`. The example
    `CMakeLists.txt` above shows how you can add VTK to the build process. The
    important lines are:

    ```
    find_package(VTK REQUIRED)
    add_definitions(-DHAVE_VTK)
    list(APPEND LIBRARIES ${VTK_LIBRARIES})
    list(APPEND INCLUDES ${VTK_INCLUDE_DIRS})
    ```

    [endsect] [section Compiling with H2Lib]

    Aboria uses [@http://www.h2lib.org H2Lib] for storing, evaluating and
    solving hierarchical matrices. You can use the bundled `FindH2Lib.cmake` to
    find H2Lib on your system. Note that you will need to download and compile
    H2Lib first, according to their instructions. You will also need to enable
    OpenMP, see below for more details on how to do this.

    ```
    set(CMAKE_MODULE_PATH  "${CMAKE_SOURCE_DIR}/Aboria/cmake"
                           ${CMAKE_MODULE_PATH})
    set(H2Lib_ROOT $ENV{HOME}/git/H2Lib)
    find_package(H2Lib REQUIRED)
    list(APPEND LIBRARIES ${H2Lib_LIBRARIES})
    list(APPEND INCLUDES ${H2Lib_INCLUDE_DIRS})
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}
    ${H2Lib_LINKER_FLAGS}") add_definitions(-DHAVE_H2LIB)
    ```
    [endsect] [section Compiling with OpenMP]

    Aboria can be run using multiple cores using OpenMP. The majority of
    parallel regions in Aboria are implemented using the OpenMP backend of
    Thrust, not raw OpenMP, so to get the best parallel performance Thrust
    should also be used.

    To add OpenMP and Thrust to the build process, add the following to your
    `CMakeLists.txt`

    ```
    # OpenMP
    find_package(OpenMP REQUIRED)
    add_definitions(-DHAVE_OPENMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

    # Thrust
    find_package(Thrust REQUIRED)
    add_definitions(-DHAVE_THRUST)
    ```

    [endsect] [section Compiling with CUDA]

    The parallel STL-like algorithms provided with Thrust are used to provide
    the majority of parallism in Aboria. It is possible to use the native CUDA
    backend of Thrust to run code on the GPU.

    [caution This mode is experimental, and has not been thoroughly tested. It
    is also quite slow, so any pull requests to fix this are most welcome!]

    There are a few restrictions on the code that you can write while using the
    CUDA backend, and these are detailed in
    [link aboria.parallelism_in_aboria.cuda] (basically if you are already
    familiar with Thrust you should be comfortable using the CUDA backend for
    Aboria). For now, just copy the following into a `getting_started.cu` file
    */

//<-
#ifdef HAVE_THRUST
    //->

    //=#include "Aboria.h"
    //=using namespace Aboria;
    //=ABORIA_VARIABLE(cu_velocity, vdouble2, "velocity")
    //=
    //=int main() {

    typedef Particles<std::tuple<cu_velocity>, 2, thrust::device_vector,
                      CellListOrdered>
        cu_container_t;
    typedef typename cu_container_t::position cu_position;

    /*
     * create a particle set with size N
     */
    cu_container_t cu_particles(N);

    /*
     * set a random position
     */
    thrust::tabulate(get<cu_position>(cu_particles).begin(),
                     get<cu_position>(cu_particles).end(),
                     [] __device__(const int i) {
                       thrust::default_random_engine gen;
                       thrust::uniform_real_distribution<float> uni(0, 1);
                       gen.discard(i);
                       return vdouble2(uni(gen), uni(gen));
                     });

    /*
     * init velocity
     */
    thrust::fill(get<cu_velocity>(cu_particles).begin(),
                 get<cu_velocity>(cu_particles).end(), vdouble2(0, 0));

    /*
     * write particle container to a vtk
     * unstructured grid file
     */
//<-
#ifdef HAVE_VTK
    //->
    vtkWriteGrid("aboria", 0, particles.get_grid(true));
//<-
#endif
    //->

    //=}
//<-
#endif
    //->

    /*`

    To enable the CUDA backend and compile the above source with the
    CUDA compiler `nvcc`, add the following to your `CMakeLists.txt`

    ```
    find_package(CUDA REQUIRED)
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}  --expt-relaxed-constexpr
                                             --expt-extended-lambda
                                             -std=c++14")

    set(SOURCE_CU
        getting_started.cu
        )

    cuda_add_executable(getting_started_cu ${SOURCE_CU})
    target_link_libraries(getting_started_cu ${LIBRARIES})
    ```

    [endsect] [section Putting it all together]

    For completness, here is a possible `CMakeLists.txt` file combining all the
    options shown above

    ```
    cmake_minimum_required(VERSION 2.8)

    set(CMAKE_MODULE_PATH  "${CMAKE_SOURCE_DIR}/Aboria/cmake"
                           ${CMAKE_MODULE_PATH})

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

    # Boost
    find_package(Boost 1.50.0 REQUIRED serialization)
    list(APPEND LIBRARIES ${Boost_LIBRARIES})
    list(APPEND INCLUDES ${Boost_INCLUDE_DIRS})

    # VTK
    find_package(VTK REQUIRED)
    if (VTK_FOUND)
        add_definitions(-DHAVE_VTK)
    endif(VTK_FOUND)
    list(APPEND LIBRARIES ${VTK_LIBRARIES})
    list(APPEND INCLUDES ${VTK_INCLUDE_DIRS})

    # H2Lib
    set(H2Lib_ROOT $ENV{HOME}/git/H2Lib)
    find_package(H2Lib REQUIRED)
    list(APPEND LIBRARIES ${H2Lib_LIBRARIES})
    list(APPEND INCLUDES ${H2Lib_INCLUDE_DIRS})
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}
    ${H2Lib_LINKER_FLAGS}") add_definitions(-DHAVE_H2LIB)

    # OpenMP
    find_package(OpenMP REQUIRED)
    add_definitions(-DHAVE_OPENMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

    # CUDA
    find_package(CUDA REQUIRED)
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}  --expt-relaxed-constexpr
                                             --expt-extended-lambda
                                             -std=c++14")

    # Thrust
    find_package(Thrust REQUIRED)
    add_definitions(-DHAVE_THRUST)

    # Aboria
    set(Aboria_LOG_LEVEL 1 CACHE STRING "Logging level (1 = least, 3 = most)")
    add_definitions(-DABORIA_LOG_LEVEL=${Aboria_LOG_LEVEL})
    list(APPEND INCLUDES Aboria/src)
    list(APPEND INCLUDES Aboria/third-party)

    include_directories(src ${INCLUDES})

    set(SOURCE
        getting_started.cpp
        )

    add_executable(getting_started ${SOURCE})
    target_link_libraries(getting_started ${LIBRARIES})

    set(SOURCE_CU
        getting_started.cu
        )

    cuda_add_executable(getting_started_cu ${SOURCE_CU})
    target_link_libraries(getting_started_cu ${LIBRARIES})
    ```

    [endsect]

    [endsect]
    */
    //]
  }
};

#endif /* DOC_GETTING_STARTED_H_ */
