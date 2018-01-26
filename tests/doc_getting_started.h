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

      
    [section Installation and Getting Started]
    [section Getting prerequisites]

    This software is tested on Ubuntu 14.04LTS with the GCC compiler (version 
    5.4.1), and Clang compiler (version 3.8.0). See our 
    [@https://travis-ci.org/martinjrobins/Aboria Travis CI page] for more details.

    You will need to install a C++ compiler with C++11 support, and the 
    [@https://cmake.org/ CMake] build software prior to using Aboria, which you can 
    do on a Debian-based OS using

    ``
    sudo apt-get install build-essential cmake
    ``

    The only required dependency is the [@http://www.boost.org Boost] library.  
    Optional dependencies are [@http://www.vtk.org/ The Visualization Toolkit] and 
    [@http://eigen.tuxfamily.org Eigen] (version >= 3.3~beta1), both of which add 
    extra functionality.  To install all these dependencies in a Debian-based OS you 
    can type

    ``
    sudo apt-get install libboost-dev libvtk5-dev libeigen3-dev
    ``

    [endsect]

    [section Compiling a program using Aboria]

    Aboria is a header-only library, so at a minimum you will need to add the 
    `Aboria/src` directory to your include path and include the `Aboria.h` header, e.g.

    ``
    #include <Aboria.h>
    ``

    The following provides a step-by-step guide on compiling your first Aboria program

    First clone the Aboria GitHub repository like so

    [teletype]
    ```
    $ git clone https://github.com/martinjrobins/Aboria
    ```
    [c++]

    Then copy and paste the code below into a C++ source file named `getting_started.cpp`.

    */

#include "Aboria.h"

using namespace Aboria;


//<-
class DocGettingStartedTest : public CxxTest::TestSuite {
public:

void test_getting_started(void) {
//->
//=int main() {
    /*
     * Create a 2d particle container type with one 
     * additional variable "velocity", represented 
     * by a 2d double vector
     */
    ABORIA_VARIABLE(velocity,vdouble2,"velocity")
    typedef Particles<std::tuple<velocity>,2> container_type;
    typedef typename container_type::position position;

    /*
     * create a particle set with size N
     */
    const int N = 100;
    container_type particles(N);

    std::uniform_real_distribution<double> uni(0,1);
    std::default_random_engine gen;
    for (int i = 0; i < N; ++i) {
        /*
         * set a random position, and initialise velocity
         */
        get<position>(particles)[i] = vdouble2(uni(gen),uni(gen));
        get<velocity>(particles)[i] = vdouble2(0,0);
    }

    /*
     * write particle container to a vtk
     * unstructured grid file
     */
//<-
#ifdef HAVE_VTK
//->
    vtkWriteGrid("aboria",0,particles.get_grid(true));
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

    # Boost
    find_package(Boost 1.50.0 REQUIRED)
    list(APPEND LIBRARIES ${Boost_LIBRARIES})
    list(APPEND INCLUDES ${Boost_INCLUDE_DIRS})

    # VTK 
    find_package(VTK REQUIRED)
    if (VTK_FOUND)
        add_definitions(-DHAVE_VTK)
    endif(VTK_FOUND)
    list(APPEND LIBRARIES ${VTK_LIBRARIES})
    list(APPEND INCLUDES ${VTK_INCLUDE_DIRS})

    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++11")

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
    [@http://eigen.tuxfamily.org Eigen] with Aboria then you will need to define the 
    `HAVE_VTK` and `HAVE_EIGEN` compiler definitions, as done using the above `CMakeLists.txt`


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
    [endsect]
    */
    //]
    }
};

#endif /* DOC_GETTING_STARTED_H_ */
