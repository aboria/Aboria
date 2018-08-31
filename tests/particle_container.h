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

#ifndef PARTICLE_CONTAINER_H_
#define PARTICLE_CONTAINER_H_

#include <cxxtest/TestSuite.h>

#include "Level1.h"

using namespace Aboria;

class ParticleContainerTest : public CxxTest::TestSuite {
public:
  template <template <typename, typename> class V,
            template <typename> class SearchMethod>
  void helper_add_particle1(void) {
    typedef Particles<std::tuple<>, 3, V, SearchMethod> Test_type;
    typedef typename Test_type::position position;
    Test_type test;
    typename Test_type::value_type p;
    get<position>(p) = vdouble3::Constant(0);
    test.push_back(p);
    TS_ASSERT_EQUALS(test.size(), 1);
  }

  template <template <typename, typename> class V,
            template <typename> class SearchMethod>
  void helper_add_particle2(void) {
    ABORIA_VARIABLE(scalar, double, "scalar")
    typedef std::tuple<scalar> variables_type;
    typedef Particles<variables_type, 3, V, SearchMethod> Test_type;
    typedef typename Test_type::position position;
    Test_type test;
    typename Test_type::value_type p;
    get<position>(p) = vdouble3::Constant(0);
    test.push_back(p);
    TS_ASSERT_EQUALS(test.size(), 1);
  }

  template <template <typename, typename> class V,
            template <typename> class SearchMethod>
  void helper_add_particle2_dimensions(void) {
    ABORIA_VARIABLE(scalar, double, "scalar")
    typedef std::tuple<scalar> variables_type;
    typedef Particles<variables_type, 6, V, SearchMethod> Test_type;
    Test_type test;
    typename Test_type::value_type p;
    typedef Vector<double, 6> double6;
    typedef position_d<6> position;
    get<position>(p) = double6::Constant(2.0);
    test.push_back(p);
    TS_ASSERT_EQUALS(test.size(), 1);
    TS_ASSERT(
        ((double6)get<position>(test[0]) == double6::Constant(2.0)).all());
  }

  template <template <typename, typename> class V,
            template <typename> class SearchMethod>
  void helper_add_delete_particle(void) {
    ABORIA_VARIABLE(scalar, double, "scalar")
    typedef std::tuple<scalar> variables_type;
    typedef Particles<variables_type, 3, V, SearchMethod> Test_type;
    Test_type test;
    typename Test_type::value_type p;
    test.push_back(p);
    test.pop_back();
    TS_ASSERT_EQUALS(test.size(), 0);
    test.push_back(p);
    test.push_back(p);
    test.push_back(p);
    TS_ASSERT_EQUALS(test.size(), 3);
    test.erase(test.begin());
    TS_ASSERT_EQUALS(test.size(), 2);
    test.erase(test.begin(), test.end());
    TS_ASSERT_EQUALS(test.size(), 0);

    typename Test_type::reference p_ref = test[0];
    get<id>(p_ref) = 101;
    typename Test_type::value_type p_value = test[0];
    TS_ASSERT_EQUALS(get<id>(p_value), 101);
  }

  void test_documentation(void) {
#if not defined(__CUDACC__)
    //[particle_container
    /*`
    [section Particle Container]
    [section Creating Particles]

    The main particles data-structure, or container, is called [classref
    Aboria::Particles]. It is templated using a tuple of variable types,
    explained below.  For example, the following creates a set of particles
    which each have (along with the standard variables such as position, id etc)
    a data package consisting of one double variable type named scalar.
    */
    //=using namespace Aboria;

    ABORIA_VARIABLE(scalar, double, "my scalar")
    typedef Particles<std::tuple<scalar>> MyParticles;
    MyParticles particles;

    /*`

    You can set the dimension of the container by using an optional unsigned
    integer template argument (defaults to 3). For example, if you wanted a
    container of particles in 2D space, you would use
    */

    typedef Particles<std::tuple<scalar>, 2> MyParticles2;
    //<-
    MyParticles2 particles_not_used2;
    //->

    /*`
    If you wanted each particle to have a potential variable held as a `double`,
    as well as a velocity variable held as a [classref Aboria::vdouble3] vector
    class, then you would write the following
    */

    ABORIA_VARIABLE(potential, double, "potential energy")
    ABORIA_VARIABLE(velocity, vdouble3, "velocity")
    typedef Particles<std::tuple<potential, velocity>> MyParticles3;
    //<-
    MyParticles3 particles_not_used3;
    //->

    /*`

    Note that there is a special case for boolean variables, which must be
    represented by an integer, rather than a boolean. This is due to the STL
    specialisation of a boolean STL vector, which conflicts with the internal
    design of Aboria. For example, here we can use an 8-bit unsigned integer to
    stand in for the boolean `flag` variable.
    */

    ABORIA_VARIABLE(flag, uint8_t, "my flag variable")
    //<-
    flag flag_not_used;
    (void)flag_not_used;
    //->

    /*`
    You can give the `MyParticles` constructor a single `int` argument to
    initialise the container with `n` particles:
    */

    const int n = 100;
    MyParticles particles2(n);

    /*`
    To create new particles simply use the `value_type` of the container type.
    For example, to create a new particle you could write
    */

    MyParticles::value_type p;

    /*`
    Each `value_type` is a tuple of values, of the types specified by each
    variable. You can retrieve or set these value using the [funcref
    Aboria::get] function, which is templated on the variable type. For example,
    say you wanted to set the `scalar` variable for particle `p`:
    */

    get<scalar>(p) = 1.0;

    /*`
    You can print the value back out, again using the [funcref Aboria::get]
    function
    */

    std::cout << "the scalar variable equals " << get<scalar>(p) << std::endl;

    /*`
    The `value_type` of the [classref Aboria::Particles Particles] container
    also has, a position, a unique id and a boolean flag indicating if this
    particle is alive or not. The position type is dependent on the dimension,
    so the best way is to get the type from the container type, i.e.
    */

    typedef MyParticles::position position;
    get<position>(p) = vdouble3(0, 0, 0);

    /*`
    Getting the id or alive flag from a `value_type` is much simpler
    */

    std::cout << "the particle id is " << get<id>(p) << std::endl;
    std::cout << "the particle alive flag is " << get<alive>(p) << std::endl;

    /*`
    Once you are happy with your particle, you can add it to the container using
    the [memberref Aboria::Particles::push_back] member function
    */

    particles.push_back(p);

    /*`
    [endsect]

    [section Multidimensional Data Types]

    Aboria provides an internal vector type [classref Aboria::Vector] for types
    representing a vector of dimension `d`. [classref Aboria::Vector] is
    templated on the type of each element and the number of dimensions:
    */

    Vector<double, 3> dim3vector;

    /*`
    There are a number of predefined `double`, `int`, and `bool` vector types,
    up to dimension 7, and typedefed by the pattern v<type><dim>. E.g. [classref
    Aboria::vdouble3], [classref Aboria::vdouble6], [classref Aboria::vint2],
    [classref Aboria::vbool5]...

    [endsect]

    [section Working with particles within the container]

    You can use the indexing operator [memberref
    Aboria::Particles::operator\[\]] to simply loop through the container
    */

    for (size_t i = 0; i < particles.size(); i++) {
      std::cout << "Accessing particle with id = " << get<id>(particles[i])
                << "\n";
    }

    /*`
    Note that the index operator [memberref Aboria::Particles::operator\[\]]
    returns a [classref Aboria::Particles::reference], which is defined as a
    tuple containing references to each of the variables. This is different from
    a reference to [classref Aboria::Particles::value_type].

    Or you can use the normal STL [memberref Aboria::Particles::begin()] and
    [memberref Aboria::Particles::end()] functions that return random access
    iterators to the beginning and end of the container.
    */

    for (auto i = particles.begin(); i != particles.end(); i++) {
      std::cout << "Accessing particle with id = " << get<id>(*i) << "\n";
    }

    /*`
    Or
    */

    for (auto i : particles) {
      std::cout << "Accessing particle with id = " << get<id>(i) << "\n";
    }

    /*`
    Or you can use the STL algorithm `for_each`. If you are using a GCC
    compiler, you can turn on the parallel mode to enable this loop to be run in
    parallel
    */

    std::for_each(particles.begin(), particles.end(), [](auto i) {
      std::cout << "Accessing particle with id = " << get<id>(i) << "\n";
    });

    /*`
    [endsect]

    [section Internal Data for Variables]

    Each variable is held internally by a STL vector `std::vector`. If you wish
    to directly access this vector, then you can use the normal [funcref
    Aboria::get] functions to get it.
    */

    std::vector<size_t> &ids = get<id>(particles);
    std::vector<double> &scalars = get<scalar>(particles);
    //<-
    (void)ids;
    (void)scalars;
    //->

    /*`
    [endsect]

    [section Particle's `value_type` versus `reference`]

    When you index an individual particle using the bracket operator [memberref
    Aboria::Particles::operator\[\]], it returns a [classref Aboria::getter_type
    getter_type], which is essentially a tuple of references to the variables
    for that particle. This [classref Aboria::getter_type getter_type] is
    `typedef`-ed to [classref Aboria::Particles::reference], and acts as the
    reference type for the container. Similarly, the [classref
    Aboria::Particles::value_type value_type] for the continer is also a
    [classref Aboria::getter_type], but instead holds a tuple of values instead
    of references.

    Reading the above paragraph, you will note the fundamental difference from
    normal STL containers, in that [classref Aboria::Particles::value_type
    value_type]& is *not the same* as [classref Aboria::Particles::value_type
    reference]. This can be relevant when writing functors for STL algorithms,
    where you will need to be sure if you need a [classref
    Aboria::Particles::value_type value_type]& or a [classref
    Aboria::Particles::value_type reference].

    For example, the `std::sort` algorithm internally stores a `value_type` of
    an element which is used in the comparison, so the functor needs to be
    equivalent to the following

    ``
    bool cmp(const value_type& a, const value_type& b)
    ``

    However, the `std::transform` algorithm can use a `unaryop` functor
    equivalent to

    ``
    Ret fun(const reference a)
    ``

    Which is more efficient than `value_type&`, since dereferencing the iterator
    will result in a `reference`.

    [note Fortunatelly, c++14 makes all this a lot easier, since you can just
    use the `auto` keyword and let the compiler deduce the correct type!]

    [endsect]

    [section Important differences from STL containers]

    The [classref Aboria::Particles] data structure acts fairly typically like a
    normal STL random-access container, with a few important differences.  It
    has methods like [memberref Aboria::Particles::push_back push_back],
    [memberref Aboria::Particles::clear clear],
    [memberref Aboria::Particles::size size],
    [memberref Aboria::Particles::erase erase].  It provides
    subtypes like [classref Aboria::Particles::value_type value_type],
    [classref Aboria::Particles::reference reference],
    [classref Aboria::Particles::const_reference const_reference],
    [classref Aboria::Particles::iterator iterator],
    [classref Aboria::Particles::const_iterator const_iterator].
    All of the normal algorithms in the standard library *should*
    work with this container, if you find any that don't please let us know and
    we will try to fix this.

    The main differences between [classref Aboria::Particles] and normal STL
    containers are:

    1. The difference between
    [classref Aboria::Particles::value_type value_type]& and
    [classref Aboria::Particles::reference reference] mentioned described
    earlier.

    2. Additional member functions are available to suit the specific purpose of
    this container, for example the
    [memberref Aboria::Particles::push_back push_back]
    function can take a vector data-type for the particle position, and the
    [memberref Aboria::Particles::get_query get_query] function for neighbour
    searching.

    3.  When using the neighbourhood searching capabilities of the container,
    the order of the particles in the particle container might change due to
    internal sorting for neighbourhood searching efficiency. So do not assume
    that the particle ordering is fixed. For example, the [memberref
    Aboria::Particles::push_back push_back] member function can reorder the
    particles if neighbourhood searching is turned on.

    [endsect]

    [section Conversion to VTK formats]

    It is possible to convert the [classref Aboria::Particles] data structure to
    a [@https://www.vtk.org/doc/nightly/html/classvtkUnstructuredGrid.html  VTK
    unstructured grid] class using the [memberref Aboria::Particles::get_grid]
    function. This function will write out each particle as a 3D point in the
    unstructured grid. By default all the particle's variables are converted
    into VTK data arrays and added to the grid, [*except] for those with names
    starting with the character "_".

    In order to write out the resultant grid to a file using the VTK data
    format, Aboria provides a useful helper function [funcref
    Aboria::vtkWriteGrid], which can write out the grid along with any constant
    fields (e.g. a timestamp) that you may need. For example, the following code
    writes out the entire contents of the the particle set to the file
    `doc00001.vtu`, along with a constant field named "time" containing the
    value 1.0.

    ```
    vtkWriteGrid("doc", 0, particles.get_grid(true), {{"time", 1.0}});
    ```

    [endsect]

    [endsect]
     */
    //]
#endif
  }

  void test_vtk_output(void) {
#ifdef HAVE_VTK
    ABORIA_VARIABLE(scalar, double, "my scalar")
    ABORIA_VARIABLE(scalar2, double, "_should not be in vtk file")
    typedef Particles<std::tuple<scalar, scalar2>> MyParticles;
    MyParticles particles(5);
    vtkWriteGrid("test", 0, particles.get_grid(true));
    vtkWriteGrid("test_with_field", 0, particles.get_grid(true),
                 {{"time", 1.0 / 2.0}});
#endif
  }

  void test_std_vector_CellList(void) {
    helper_add_particle1<std::vector, CellList>();
    helper_add_particle2<std::vector, CellList>();
    helper_add_particle2_dimensions<std::vector, CellList>();
    helper_add_delete_particle<std::vector, CellList>();
  }

  void test_std_vector_CellListOrdered(void) {
    helper_add_particle1<std::vector, CellListOrdered>();
    helper_add_particle2<std::vector, CellListOrdered>();
    helper_add_particle2_dimensions<std::vector, CellListOrdered>();
    helper_add_delete_particle<std::vector, CellListOrdered>();
  }

  void test_thrust_vector_CellListOrdered(void) {
#if defined(__CUDACC__)
    helper_add_particle1<thrust::device_vector, CellListOrdered>();
    helper_add_particle2<thrust::device_vector, CellListOrdered>();
    helper_add_particle2_dimensions<thrust::device_vector, CellListOrdered>();
    helper_add_delete_particle<thrust::device_vector, CellListOrdered>();
#endif
  }
};

#endif /* PARTICLE_CONTAINER_H_ */
