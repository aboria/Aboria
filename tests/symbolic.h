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

#ifndef SYMBOLICTEST_H_
#define SYMBOLICTEST_H_

#include <cxxtest/TestSuite.h>

#include "Aboria.h"

using namespace Aboria;

class SymbolicTest : public CxxTest::TestSuite {
public:
  void test_documentation(void) {
    ABORIA_VARIABLE(flag, uint8_t, "my flag");
    size_t N = 100;
    typedef Particles<std::tuple<flag>> Particles_t;
    typedef typename Particles_t::position position;
    Particles_t particles(N);

    std::uniform_real_distribution<double> uniform(0, 1);
    for (size_t i = 0; i < N; ++i) {
      auto &gen = get<generator>(particles)[i];
      get<position>(particles)[i] =
          vdouble3(uniform(gen), uniform(gen), uniform(gen));
    }

    particles.init_neighbour_search(vdouble3::Constant(0),
                                    vdouble3::Constant(1),
                                    vdouble3::Constant(false));

    //[symbolic
    /*`

    [section Symbolic Expressions]

    [section Setup]
    To start using symbolic expressions, you first need to define a set of
    symbols to represent your variables, as well as labels to represent your
    particle set(s).

    A symbol representing each particle's `position` is defined as
    */

    Symbol<position> p;

    /*`
    A symbol representing the variable `alive` is defined as
    */

    Symbol<alive> _alive;

    /*`
    A label representing the particle set `particles` with type `MyParticles` is
    defined as
    */

    Label<0, Particles_t> a(particles);
    Label<1, Particles_t> b(particles);

    /*`
    The first template argument is the *id* of the label, and the second is the
    type of the particle set label refers to. Note that the type of each label
    must be unique, as this type is used to determine which particle you mean
    when using the label within expressions. You can use the *id* template
    argument to ensure that each label type is unique.

    Labels refer to a specific particle set. For example, given a bivariate
    neighbour expression involving two particles from `particles`, the label `a`
    defined above would refer to the first particle, and `b` would refer to the
    second. Note that the result values of a bivariate expression will form a
    matrix of values, with particle `a` corresponding to the row of the matrix,
    and particle `b` corresponding to the column.

    [endsect]

    [section Constant Expressions]

    Now we have defined our labels and symbols, we can create a simple
    expression to set the position of all particles to `double3(0,0,1)`
    */

    p[a] = vdouble3(0, 0, 1);

    /*`
    If we want to add `1` to every particle position, we can use an increment
    expression
    */

    p[a] += 1;

    /*`
    A special case involves the `alive` variable flag. If we use our `_alive`
    symbol to set all the alive flag's to `false` like so
    */

    _alive[a] = false;

    /*`
    Then after the expression is complete Aboria will call the [memberref
    Aboria::Particles::delete_particles delete_particles] member function to
    delete all the particles with `get<alive>(particle) == false` (i.e. all of
    them). After this expression the particle container will be empty.

    [endsect]


    [section Univariate Expressions]

    Single-particle dependent expressions can easily be built by using a single
    label on the RHS of the expression. For example we can also add a constant
    value `1` to each particle position like so
    */

    p[a] = p[a] + 1;

    /*`
    Or we can use a flag variable `flag`, along with its symbol `f`, so define
    two sets of particles within the same container and only change the particle
    position if `get<flag>(particle) == true`. Recall that we can't use boolean
    flags so we can instead use a `uint8_t` datatype. For example

    ``
    ABORIA_VARIABLE(flag,uint8_t,"my flag");

    // create particle container and set flags and positions here
    ``

    */

    Symbol<flag> f;
    p[a] = if_else(f[a], p[a] + 1, p[a]);

    /*`
    We can even selectively delete particle using our `f` and `_alive` symbols.
    */

    _alive[a] = if_else(f[a], true, false);

    /*`
    [endsect]

    [section Bivariate Expressions]

    So far we have used constant or univariate expressions on the RHS of an
    assignment operator. This makes sense because we can only assign an
    expression to a single particle if that expression depends on a constant or
    that particle's variables. However, we might also want our expression to
    depend on the other particles in the container, or another container. In
    this case we wish to use a bivariate expression. For example, the following
    equation defined a sum over all other particles (with label `b`), using the
    function `w`, which depends on the separation of particle `a` and `b`.

    $$
    p_a = \sum_b w(p_b-p_a)
    $$

    Normally for a bivariate expression we wish to accumulate the contribution
    from other particles in order to arrive at a result. This is most often done
    with a summation. In Aboria, we can define an accumulator object, which
    takes a single template argument which is the function or functor to
    accumulate with. For example, the following defines a summation accumulator
    using `std::plus`
    */

    Accumulate<std::plus<vdouble3>> sum;

    /*`
    Once we have this summation accumulator, we can write an expression similar
    to the equation above like so, using $w(p_b) = (p_b - p_a)$.Note that the
    first argument to `sum` is set to `true`, indicating the summation is over
*all* particles in label `b`.
*/

    p[a] = sum(b, p[b] - p[a]);

    /*`
     We might also want to sum the inverse distance between particle pairs,
     like so
     */

    p[a] = sum(b, 1.0 / (p[b] - p[a]));

    /*`
    Unfortunately if `a` and `b` refer to the same particle container this will
    result a divide-by-zero at runtime, when  `a` and `b` are labels for the
    same particle. Therefore we can restrict the evaluation of the sum by
    setting the second argument to ensure that the `id` of particles `a` and `b`
    are non-identical. Recall that `id` is a built-in variable that contains a
    unique id for each particle in the container.
    */

    Symbol<id> _id;
    p[a] = sum(b, if_else(_id[a] != _id[b], 1.0, 0.0) / (p[b] - p[a]));

    /*`
    So the first argument to `sum` is the label to sum over, the second is the
    conditional expression that must evaluate to `true` to be included in the
    summation, and the third is an expression that provides the value to be
    included in the summation.

    There is a special case for the conditional expression, when you want to sum
    over all particles within a certain radius. This can be expressed using the
    Aboria::AccumulateWithinDistance summation
    */

    AccumulateWithinDistance<std::plus<vdouble3>> sum_within_two(2);
    auto dx = create_dx(a, b);
    p[a] = sum_within_two(b, dx / pow(norm(dx), 2));

    /*`
    where `dx` is a [classref Aboria::Dx Dx] symbol representing the
    *shortest* vector from `b` to `a`. Note that this might be different
    from `p[a]-p[b]` for periodic domains. The symbolic function
    [functionref Aboria::norm norm] returns the 2-norm, or magnitude, of the
    vector `dx`, and the symbolic function [functionref Aboria::pow pow]
    returns the first argument to the power of the second (in much the same
    way as std::pow, but lazily evaluated).

    In this case Aboria will perform a summation over neighbours closer
    than a radius of 2, and will use the neighbourhood searching facility
    described in [link aboria.neighbourhood_searching] to find these
    neighbouring particles. Note that you need to call [memberref
    Aboria::Particles::init_neighbour_search] before any bivariate
    expressions using neighbourhood searching.

        As another, more complete example, we can use the
    neighbourhood-searching expressions to count the number of particles
    within a distance of `2` of each individual particle, storing the result
    in a variable called `count`.
    */

    ABORIA_VARIABLE(count, int, "number of surrounding particles")
    typedef Particles<std::tuple<count>> MyParticles;
    typedef MyParticles::position position_t;

    // add some particles
    MyParticles particles2(N);

    std::uniform_real_distribution<double> uni(0, 1);
    for (size_t i = 0; i < N; ++i) {
      auto &gen = get<generator>(particles2)[i];
      get<position_t>(particles2)[i] = vdouble3(uni(gen), uni(gen), uni(gen));
    }

    // initialise neighbour searching
    particles2.init_neighbour_search(
        vdouble3::Constant(0), vdouble3::Constant(1), vbool3::Constant(false));

    // define symbols and labels, and sum
    Symbol<count> c;
    Label<0, MyParticles> i(particles2);
    Label<1, MyParticles> j(particles2);
    AccumulateWithinDistance<std::plus<int>> neighbour_sum(2);
    auto dx_ij = create_dx(i, j);

    // count neighbouring particles within a distance of 2
    c[i] = neighbour_sum(j, 1);

    /*`
    [endsect]
    [endsect]
     */

    //]
  }

  void helper_create_double_vector(void) {
    ABORIA_VARIABLE(scalar, double, "scalar")
    typedef Particles<std::tuple<scalar>> ParticlesType;
    ParticlesType particles;

    Symbol<scalar> s;
    Label<0, ParticlesType> a(particles);
  }

  void helper_create_default_vectors(void) {
    ABORIA_VARIABLE(scalar, double, "scalar")
    typedef Particles<std::tuple<scalar>> ParticlesType;
    typedef position_d<3> position;
    ParticlesType particles;

    Symbol<position> p;
    Symbol<id> id_;
    Symbol<alive> alive_;
  }

  void helper_transform(void) {
    ABORIA_VARIABLE(scalar, double, "scalar")
    typedef Particles<std::tuple<scalar>> ParticlesType;
    typedef position_d<3> position;
    ParticlesType particles;

    Symbol<position> p;
    Symbol<id> id_;
    Symbol<alive> alive_;
    Symbol<scalar> s;
    Label<0, ParticlesType> a(particles);

    s[a] = 0;

    ParticlesType::value_type particle;
    particles.push_back(particle);
    particles.push_back(particle);

    s[a] = 0;

    TS_ASSERT_EQUALS(get<scalar>(particles[0]), 0);
    TS_ASSERT_EQUALS(get<scalar>(particles[1]), 0);

    s[a] = 1;

    TS_ASSERT_EQUALS(get<scalar>(particles[0]), 1);
    TS_ASSERT_EQUALS(get<scalar>(particles[1]), 1);

    s[a] = s[a] + 1;

    TS_ASSERT_EQUALS(get<scalar>(particles[0]), 2);
    TS_ASSERT_EQUALS(get<scalar>(particles[1]), 2);

    s[a] = s[a] + 1;

    TS_ASSERT_EQUALS(get<scalar>(particles[0]), 3);
    TS_ASSERT_EQUALS(get<scalar>(particles[1]), 3);

    p[a] = vdouble3(1, 2, 3);

    TS_ASSERT_EQUALS(get<position>(particles[0])[0], 1);
    TS_ASSERT_EQUALS(get<position>(particles[0])[1], 2);
    TS_ASSERT_EQUALS(get<position>(particles[0])[2], 3);

    p[a] += vdouble3(1, 2, 3);

    TS_ASSERT_EQUALS(get<position>(particles[0])[0], 2);
    TS_ASSERT_EQUALS(get<position>(particles[0])[1], 4);
    TS_ASSERT_EQUALS(get<position>(particles[0])[2], 6);

    p[a] = p[a] * s[a];

    TS_ASSERT_EQUALS(get<position>(particles[0])[0], 6);
    TS_ASSERT_EQUALS(get<position>(particles[0])[1], 12);
    TS_ASSERT_EQUALS(get<position>(particles[0])[2], 18);

    p[a] = if_else(id_[a] == 0, vdouble3(0, 0, 0), vdouble3(3, 2, 1));

    TS_ASSERT_EQUALS(get<position>(particles[0])[0], 0);
    TS_ASSERT_EQUALS(get<position>(particles[0])[1], 0);
    TS_ASSERT_EQUALS(get<position>(particles[0])[2], 0);

    TS_ASSERT_EQUALS(get<position>(particles[1])[0], 3);
    TS_ASSERT_EQUALS(get<position>(particles[1])[1], 2);
    TS_ASSERT_EQUALS(get<position>(particles[1])[2], 1);
  }

  void helper_neighbours(void) {
    ABORIA_VARIABLE(scalar, double, "scalar")

    typedef Particles<std::tuple<scalar>> ParticlesType;
    typedef position_d<3> position;
    ParticlesType particles;

    vdouble3 min(-1, -1, -1);
    vdouble3 max(1, 1, 1);
    vdouble3 periodic(true, true, true);
    double diameter = 0.1;
    particles.init_neighbour_search(min, max, periodic);

    Symbol<position> p;
    Symbol<id> id_;
    Symbol<alive> alive_;
    Symbol<scalar> s;
    Label<0, ParticlesType> a(particles);
    Label<1, ParticlesType> b(particles);
    auto dx = create_dx(a, b);
    AccumulateWithinDistance<std::plus<double>> sum(diameter);

    particles.push_back(vdouble3(0, 0, 0));
    particles.push_back(vdouble3(diameter * 2, 0, 0));

    s[a] = sum(b, 1);

    TS_ASSERT_EQUALS(get<scalar>(particles[0]), 1);
    TS_ASSERT_EQUALS(get<scalar>(particles[1]), 1);

    p[a] =
        if_else(id_[a] == 0, vdouble3(diameter / 2.0, 0, 0), vdouble3(0, 0, 0));

    TS_ASSERT_EQUALS(get<position>(particles[0])[0], diameter / 2.0);
    TS_ASSERT_EQUALS(get<position>(particles[0])[1], 0);
    TS_ASSERT_EQUALS(get<position>(particles[0])[2], 0);

    TS_ASSERT_EQUALS(get<position>(particles[1])[0], 0);
    TS_ASSERT_EQUALS(get<position>(particles[1])[1], 0);
    TS_ASSERT_EQUALS(get<position>(particles[1])[2], 0);

    p[a] = 0.5 * (p[a] + s[a]);

    TS_ASSERT_EQUALS(get<position>(particles[0])[0],
                     0.5 * (diameter / 2.0 + 1));
    TS_ASSERT_EQUALS(get<position>(particles[0])[1], 0.5);
    TS_ASSERT_EQUALS(get<position>(particles[0])[2], 0.5);

    TS_ASSERT_EQUALS(get<position>(particles[1])[0], 0.5);
    TS_ASSERT_EQUALS(get<position>(particles[1])[1], 0.5);
    TS_ASSERT_EQUALS(get<position>(particles[1])[2], 0.5);

    s[a] = sum(b, 1);

    TS_ASSERT_EQUALS(get<scalar>(particles[0]), 2);
    TS_ASSERT_EQUALS(get<scalar>(particles[1]), 2);

    Accumulate<std::plus<double>> sum_all;
    static_assert(detail::is_const<decltype(sum_all(a, sum(b, 1)))>::value,
                  "result of is_const on expression "
                  "sum_all(a,sum(b,1)) is not true");

    const double sum_of_sum = eval(sum_all(a, sum(b, 1.0)));
    TS_ASSERT_EQUALS(sum_of_sum, 4);

    p[a] = if_else(id_[a] == 0, vdouble3(0, 0, 0),
                   vdouble3(diameter / 2.0, diameter / 2.0, diameter / 2.0));
    s[a] = if_else(id_[a] == 0, 0, 1);

    TS_ASSERT_EQUALS(get<scalar>(particles[0]), 0);
    TS_ASSERT_EQUALS(get<scalar>(particles[1]), 1);

    TS_ASSERT_EQUALS(get<position>(particles[0])[0], 0);
    TS_ASSERT_EQUALS(get<position>(particles[0])[1], 0);
    TS_ASSERT_EQUALS(get<position>(particles[0])[2], 0);

    TS_ASSERT_EQUALS(get<position>(particles[1])[0], diameter / 2.0);
    TS_ASSERT_EQUALS(get<position>(particles[1])[1], diameter / 2.0);
    TS_ASSERT_EQUALS(get<position>(particles[1])[2], diameter / 2.0);

    AccumulateWithinDistance<std::plus<vdouble3>> sumVect(diameter);
    Accumulate<std::plus<vdouble3>> sum_all_vect;

    const vdouble3 sum_of_sumv = eval(sum_all_vect(
        a, sumVect(b, vdouble3(0, 0, 0) + 0.5 * (s[a] / 2.0 + s[b] / 10.0))));

    for (int i = 0; i < 3; ++i) {
      TS_ASSERT_DELTA(sum_of_sumv[i], 0.6,
                      10 * std::numeric_limits<double>::epsilon());
    }

    p[a] = sumVect(b, vdouble3(0, 0, 0) + 0.5 * (s[a] / 2.0 + s[b] / 10.0));

    TS_ASSERT_EQUALS(get<position>(particles[0])[0], 0.05);
    TS_ASSERT_EQUALS(get<position>(particles[0])[1], 0.05);
    TS_ASSERT_EQUALS(get<position>(particles[0])[2], 0.05);

    TS_ASSERT_EQUALS(get<position>(particles[1])[0], 0.55);
    TS_ASSERT_EQUALS(get<position>(particles[1])[1], 0.55);
    TS_ASSERT_EQUALS(get<position>(particles[1])[2], 0.55);

    //
    // test inf norm range sum
    //
    AccumulateWithinDistance<std::plus<double>, -1> box_sum(diameter);
    get<position>(particles)[0] = vdouble3(0, 0, 0);
    get<position>(particles)[1] = 0.99 * vdouble3(diameter, diameter, diameter);
    particles.update_positions();

    s[a] = sum(b, 1);

    TS_ASSERT_EQUALS(get<scalar>(particles[0]), 1);
    TS_ASSERT_EQUALS(get<scalar>(particles[1]), 1);

    s[a] = box_sum(b, 1);

    TS_ASSERT_EQUALS(get<scalar>(particles[0]), 2);
    TS_ASSERT_EQUALS(get<scalar>(particles[1]), 2);
  }

  void helper_level0_expressions(void) {
    ABORIA_VARIABLE(scalar, double, "scalar")

    typedef Particles<std::tuple<scalar>> ParticlesType;
    typedef position_d<3> position;
    ParticlesType particles;

    Symbol<position> p;
    Symbol<id> id_;
    Symbol<alive> alive_;
    Symbol<scalar> s;
    Label<0, ParticlesType> a(particles);
    Label<1, ParticlesType> b(particles);
    auto dx = create_dx(a, b);

    Accumulate<std::plus<double>> sum;
    Accumulate<Aboria::max<double>> max;
    max.set_init(-1);
    Accumulate<Aboria::min<double>> min;
    min.set_init(1000);

    particles.push_back(vdouble3(0, 0, 0));
    particles.push_back(vdouble3(2, 0, 0));

    double result = eval(sum(a, if_else(p[a][0] < 1, 1, 0)));
    TS_ASSERT_EQUALS(result, 1);
    result = eval(sum(a, 1));
    TS_ASSERT_EQUALS(result, 2);
    result = eval(sum(a, p[a][0]));
    TS_ASSERT_EQUALS(result, 2);
    int result2 = eval(max(a, id_[a]));
    TS_ASSERT_EQUALS(result2, 1);
    result2 = eval(min(a, id_[a]));
    TS_ASSERT_EQUALS(result2, 0);
    particles.push_back(vdouble3(0, 0, 0));
    result2 = eval(max(a, id_[a]));
    TS_ASSERT_EQUALS(result2, 2);
  }

  void test_default() {
    helper_create_default_vectors();
    helper_create_double_vector();
    helper_transform();
    helper_neighbours();
    helper_level0_expressions();
  }
};

#endif /* SYMBOLICTEST_H_ */
