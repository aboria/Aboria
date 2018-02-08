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

#ifndef OPERATORSTEST_H_
#define OPERATORSTEST_H_

#include <cxxtest/TestSuite.h>

#include "Aboria.h"

using namespace Aboria;

class OperatorsTest : public CxxTest::TestSuite {
public:
  void test_documentation(void) {
#ifdef HAVE_EIGEN
    //[operators
    /*`
    [section Evaluating and Solving Kernel Operators with Eigen]

    Given that Aboria can describe non-linear operators, this naturally covers
    linear operators (i.e. matrices) as well. Consider the summation operator
    given by the kernel function $K(x_i,x_j)$, over a set of $N$ particles

    $$
    a_i = \sum_j^N b_j K(x_i,x_j) \text{ for } i=1..N
    $$

    This is a common enough operator that can be used in many areas. If
    $K(x_i,x_j) = 1/(x_j-x_i)$, then the operator might be calculating the
    force on a set of charged particles via a Coulomb force. If
    $K(x_i,x_j) = \sqrt{(x_j-x_i)^2 + c^2}$, then the operator might be used
    for function interpolation using the multiquadric basis function.

    One way to evaluate this operator is to use a matrix to store the values of
    $K(x_i,x_j)$ for each particle pair, leading to a matrix $\mathbf{K}$ with
    storage size $N^2$. Then the summation operator above is equivalent to the
    matrix-vector product

    $$
    \mathbf{a} = \mathbf{K} \mathbf{b}.
    $$

    However, $\mathbf{K}$ could be too large to fit in memory, or the values of
    $x_i$ and $x_j$ might  change too frequently for this to be useful. Or you
    may wish to take advantage of the fact that $K(x_i,x_j)$ is a continuous
    function and use method such as Chebyshev interpolation or Fast Multipole
    methods to efficiently calculate the action of the operator on the vector
    $\mathbf{b}$. For any of these reasons and many others, Aboria can help you
    out.

    Aboria provides functionality to describe linear operators arising from the
    evaluation of kernel functions at a set of points (i.e. particles in
    Aboria's terminology) in $N$ dimensional space. From now on we will refer to
    these types of operators as [*kernel operators]. To provide the concept and
    API of a matrix or linear operator, we will use the C++ linear algebra
    library
    [@eigen.tuxfamily.org Eigen]. Aboria provides functionality to wrap kernel
    operators as [classref Aboria::MatrixReplacement], so that Eigen can treat
    them as normal dense or sparse matrices.

    Using this wrapping, we can take advantage of Eigen to:

    # calculate the action of a kernel operator on a vector (i.e. a
    matrix-vector multiplication) # use Eigen's iterative solvers to solve a
    large linear system of equations arising from a kernel operator.

    [section Creating Dense Operators]

    The most general case involves a kernel operator $K$ that is non-zero
    for every possible particle pair. The kernel function can depend on the
    particle positions and/or the variables assigned to each particle. For
    example, say we had a particle set with particle positions $\mathbf{x}_i$
    for $i=1..N$, and a single variable $a_i$.  We wish to create a summation
    operator using the kernel function

    $$
    K(\mathbf{x}_i,a_i,\mathbf{x}_j,a_j) = \frac{a_i
    a_j}{||\mathbf{x}_j-\mathbf{x}_i|| + \epsilon}
    $$

    were $||.||$ refers to the 2-norm, or magnitude of a vector.

    First we need a particle set to apply the operator to. We will create a
    particle set containing $N=100$ particles with a single additional variable
    $a$.
    */

    const size_t N = 100;
    const double epsilon = 0.1;
    ABORIA_VARIABLE(a, double, "a");
    typedef Particles<std::tuple<a>> particle_type;
    typedef particle_type::position position;
    particle_type particles(N);
    std::default_random_engine gen;
    std::uniform_real_distribution<double> uniform(0, 1);
    for (size_t i = 0; i < N; ++i) {
      get<position>(particles)[i] =
          vdouble3(uniform(gen), uniform(gen), uniform(gen));
      get<a>(particles)[i] = uniform(gen);
    }
    /*`
    For convenience, we will also define a constant reference type to refer to
    each particle in the container
    */
    typedef particle_type::const_reference const_particle_reference;

    /*`
    We then create a dense kernel operator using the
    [funcref Aboria::create_dense_operator] function
    */
    auto K = create_dense_operator(
        particles, particles,
        [epsilon](const_particle_reference i, const_particle_reference j) {
          const auto dx = get<position>(j) - get<position>(i);
          return (get<a>(i) * get<a>(j)) / (dx.norm() + epsilon);
        });

    /*`
    Note that [funcref Aboria::create_dense_operator] takes three arguments. The
    first two are particle containers which give the two particle sets involved
    in the operator. The first container holds the particles indexed by $i$ in
    the kernel function, and the second holds the particles indexed by $j$. For
    a matrix representation, you might say that these refer to the rows and
    columns of the matrix.

    The third argument to [funcref Aboria::create_dense_operator] can be a
    function object, or C++ lambda expression. Basically any valid C++ object
    that can be called with two arguments, the first of type
    `const_particle_reference` (i.e. a constant reference to a particle in the
    set indexed by $i$), and the second of type `const_particle_reference` (i.e.
    a constant reference to a particle in the set indexed by $j$). Note that in
    this case the particle sets indexed by $i$ and $j$ is the same particle set
    `particles`. However, in other cases you may want $i$ and $j$ to index
    different particle sets, in which case the types for arguments 1 and 2 could
    be different.

    In the code above, we are using a lambda expression as our function object,
    and create one that returns the particular kernel function
    $K(\mathbf{x}_i,a_i,\mathbf{x}_j,a_j)$.

    Once we have created the operator `K`, we can use it within Eigen as if it
    were a normal matrix.  For example, to apply `K` to a vector `b`, we could
    write the following

    */
    Eigen::VectorXd b = Eigen::VectorXd::LinSpaced(N, 0, 1.0);
    Eigen::VectorXd c_1 = K * b;
    /*`
    This line of code calculates the following

    $$
    c_i = \sum_j^N b_j K(\mathbf{x}_i,a_i,\mathbf{x}_j,a_j) \text{ for } i=1..N
    $$

    Note that rather then storing all the values of K that are needed for this
    summation, Aboria will instead evaluate these values as they are needed.
    Therefore the memory requirements are only $\mathcal{O}(N)$, rather than
    $\mathcal{O}(N^2)$ for a traditional matrix. However, this requires
    evaluating the kernel function for each pair of particles at a cost of
    $\mathcal{O}(N^2)$. If you wish to calculate this operation approximately,
    then Aboria can perform the same operation using the Fast Multipole Method
    or Hierarchical Matrices, please see subsequent sections for more details on
    how to do this.

    [caution the [classref Aboria::MatrixReplacement] operator `K` cannot be
    used, for example, in multiplications or additions with other Eigen
    matrices. Thus far, it has only really been tested with matrix-vector
    multiplication and Eigen's iterative solvers]

    If we wish to perform the same operator, but using a traditional matrix, we
    can use `K`'s [memberref Aboria::MatrixReplacement::assemble] member
    function to fill in a normal Eigen matrix with the values of the kernel
    function $K(\mathbf{x}_i,a_i,\mathbf{x}_j,a_j)$. This might be useful if you
    wish to perform the same operation repeatedly, or in one of Eigen's direct
    solvers.

    */
    Eigen::MatrixXd K_eigen(N, N);
    K.assemble(K_eigen);

    Eigen::VectorXd c_2 = K_eigen * b;

    //<-
    TS_ASSERT(c_2.isApprox(c_1));
    //->
    /*`
    [endsect]
    [section Creating Sparse Operators]

    Is is common in particle-based methods that the kernel function $K$ be
    non-zero only for particle pairs separated by less than a certain radius. In
    this case we have a summation operation like so

    $$
    c_i = \sum_j^N b_j K_s(\mathbf{x}\_i,a_i,\mathbf{x}\_j,a_j)  \text{ for }
    i=1..N
    $$

    where $K_s$ is a truncated version of $K$ that is only non-zero for
    $||\mathbf{dx}\_{ij}||<r$, where $\mathbf{dx}\_{ij}$ is the shortest vector
    between particles $i$ and $j$.  Note that for non-periodic systems, this
    will be
    $\mathbf{dx}\_{ij}=\mathbf{x}_j-\mathbf{x}_i$.


    $$
    K_s(\mathbf{x}\_i,a_i,\mathbf{x}\_j,a_j) =
        \begin{cases}
            K(\mathbf{x}\_i,a_i,\mathbf{x}\_j,a_j), & \text{for }
    ||\mathbf{dx}\_{ij}||<r \\\\\ 0 & \text{otherwise}. \end{cases}
    $$

    Since the summation is only non-zero for $||\mathbf{dx}_{ij}||<r$, we wish
    to aim for better than $\mathcal{O}(N^2)$ time and combine the sum with a
    spatial search of radius $r$.

    Lets assume that we wish a similar kernel function as before

    $$
    K(\mathbf{x}\_i,a_i,\mathbf{x}\_j,a_j) = \frac{a_i
    a_j}{||\mathbf{dx}\_{ij}|| + \epsilon}
    $$

    We can create the operator `K_s` in Aboria like so (setting $r=0.1$ in this
    case)
    */
    const double r = 0.1;
    auto K_s = create_sparse_operator(
        particles, particles, r,
        [epsilon](const vdouble3 &dx, const_particle_reference i,
                  const_particle_reference j) {
          return (get<a>(i) * get<a>(j)) / (dx.norm() + epsilon);
        });

    /*`
    Sparse operators support periodic domains, and therefore we now need the
    `dx` argument in the kernel function object. This is a distance value that
    gives the shortest vector from particle `i` to particle `j`.

    When applied to a vector, this operator will use the neighbour search of the
    `particles` container to perform a neighbour search for all particle pairs
    where
    $||\mathbf{dx}\_{ij}||<r$.

    Before we can use this operator, we need to make sure that the neighbour
    search for `particles` is initialised. By default, the particle container
    was created using three spatial dimensions, so we need to set up a domain
    from $(0,0,0)$ to $(1,1,1)$ which is not periodic in all three directions.

    */

    vdouble3 min = vdouble3::Constant(0);
    vdouble3 max = vdouble3::Constant(1);
    vbool3 periodic = vbool3::Constant(false);
    particles.init_neighbour_search(min, max, periodic);

    /*`

    Once this is done, we can then apply the operator to the vector `b` from
    before

    */

    Eigen::VectorXd c_3 = K_s * b;

    /*`

    Once again, we can write out `K_s` to a traditional matrix. This time, we
    will write out the values of `K_s` to a sparse matrix, so we can still
    obtain an efficient operator

    */
    Eigen::SparseMatrix<double> K_s_eigen(N, N);
    K_s.assemble(K_s_eigen);

    Eigen::VectorXd c_4 = K_s_eigen * b;
    //<-
    TS_ASSERT(c_4.isApprox(c_3));
    //->
    /*`
    [endsect]


    [section Creating Chebyshev Operators]

    Lets assume that the kernel function of interest depends only on the
    particle positions

    $$
    c_i = \sum_j^N b_j K(\mathbf{x}_i,\mathbf{x}_j) \text{ for } i=1..N
    $$

    and that we can approximate $K(\mathbf{x}_i,\mathbf{x}_j)$ using
    interpolation. That is, if $w_l(x)$ denotes a set of interpolating
    functions, then

    $$
    K(\mathbf{x}_i,\mathbf{x}_j) \approx \sum_l \sum_m
    K(\mathbf{x}_l,\mathbf{x}_m) w_l(\mathbf{x}_i) w_m(\mathbf{x}_j)
    $$

    Using this idea, we can use chebyshev interpolation to interpolate
    $K(\mathbf{x}_i,\mathbf{x}_j)$ at the chebyshev nodes, leading to an
    operator that can be applied to a vector at a reduced cost. How reduced
    depends on the number of chebyshev nodes that are chosen. A small number of
    nodes means that the error in the approximation grows, while a larger number
    of nodes means increased computational cost.

    Note that this idea has been published in the following paper, which was
    used as a reference for Aboria's implementation:

    Fong, William, and Eric Darve. "The black-box fast multipole method."
    Journal of Computational Physics 228.23 (2009): 8712-8725.

    One restriction on the kernel function in this context is that it must be a
    smooth function of only position. With this in mind, we will define the
    following kernel function

    $$
    K(\mathbf{x}\_i,\mathbf{x}\_j) = \sqrt{||\mathbf{dx}\_{ij}|| + \epsilon}
    $$

    which is the multiquadric radial basis function.

    A kernel operator using this function and chebyshev interpolation can be
    created in Aboria using the [funcref Aboria::create_chebyshev_operator]
    function

    */
    const unsigned int n = 10;
    auto K_c = create_chebyshev_operator(
        particles, particles, n,
        [epsilon](const vdouble3 &i, const vdouble3 &j) {
          return std::sqrt((j - i).norm() + epsilon);
        });

    /*`

    where `n` sets the number of chebyshev nodes in each dimension. Here the
    function object given to [funcref Aboria::create_chebyshev_operator] must
    have two arguments representing the positions $\mathbf{x}\_i$ and
    $\mathbf{x}\_j$.

    [caution Once `K_c` is created, it assumes that the positions of the
    particles in `particles` do not change. If this is not true, then you can
    use the function of the underlying kernel class [classref
    Aboria::KernelChebyshev] to update the positions.]

    [note Aboria's neighbourhood searching does not need to be initialized in
    this case, as no neighbourhood queries are used.]

    As before, this operator can be applied to an Eigen vector

    */

    Eigen::VectorXd c_5 = K_c * b;

    /*`


    [endsect]

    [section Creating Fast Multipole Method Operators]

    One disadvantage of using chebyshev interpolation to construct an operator
    is that it is not valid if the kernel function is discontinuous. This is
    violated by many kernels which contain singularities when $\mathbf{x}_i =
    \mathbf{x}_j$. Also, often large numbers of chebyshev nodes are required,
    which does not reduce the computational cost sufficiently while retaining
    accuracy.

    Another powerful class of methods are based on the Fast Multipole Method,
    which has been also implemented in Aboria. These are valid for kernel
    functions with singularities and which boasts $\mathcal{O}(N)$ complexity.

    Once again, we will use chebyshev interpolation, but now will construct a
    tree data structure using one of Aboria's tree data structures (kd-tree or
    oct-tree). When the kernel operator is applied to a vector, the fast
    multiplole algorithm will be invoked, which leads to traversals up and down
    the tree, using chebyshev interpolation to compress the kernel function for
    clusters of particles that are ['well separated]. The interaction of
    particle clusters that are not well separated are calculated directly. The
    details of this are in the following article, which was used as a reference
    for Aboria's implementation:

    Fong, William, and Eric Darve. "The black-box fast multipole method."
    Journal of Computational Physics 228.23 (2009): 8712-8725.

    A kernel operator using the fast multipole method can be created in Aboria
    using the [funcref Aboria::create_fmm_operator] function

    */
    const unsigned int n2 = 4;
    auto K_f = create_fmm_operator<n2>(
        particles, particles,
        [epsilon](const vdouble3 &i, const vdouble3 &j) {
          return std::sqrt((j - i).norm() + epsilon);
        },
        [epsilon](const_particle_reference i, const_particle_reference j) {
          return std::sqrt(
              ((get<position>(j) - get<position>(i)).norm() + epsilon));
        });

    /*`

    where `n2` sets the number of chebyshev nodes in each dimension. Note that
    this is typically much less than what is required by the chebyshev operator
    discussed in the previous section. Here we need to pass two function objects
    to [funcref Aboria::create_fmm_operator], both of which represent the kernel
    applied in different situations. The first takes two positions that are
    guarenteed to be not equal, and which in FMM terminology are used to create
    the M2L (multipole to local) operators. The second function object takes two
    particle references from the row and column particle sets, and is used to
    create the P2P (particle to particle) operator. Note that the position of
    `i` and `j`, and/or the particles themselves, might be identical. Therefore
    this function object must handle any kernel singularities.

    [note the [funcref Aboria::create_fmm_operator] function creates a kernel
    operator in which  the fast multipole algorithm is applied completely "on
    the fly". The fmm can be significantly speeded up for static particle
    positions by storing a hierarchical matrix. Please see the next section for
    more details of how to do this in Aboria.]

    As before, this operator can be applied to an Eigen vector

    */

    Eigen::VectorXd c_6 = K_f * b;

    /*`




    [endsect]


    */
    //<-
#ifdef HAVE_H2LIB
    //->
    /*`
    [section Creating Hierarchical Matrix Operators]

    The fast multipole method (fmm) results in multiple tree traversals while
    applying linear operators to the multipole vectors. Rather than calculating
    these linear operators during the algorithm, it is also possible to
    pre-calculate them and store them in a hierarchical H2 matrix. The advantage
    of this is that applying the operator to a vector is much faster, the
    downside is that it assumes that the particle positions do not change.

    A kernel operator which does this can be created using the
    [funcref Aboria::create_h2_operator] function

    */
    auto K_h = create_h2_operator(
        particles, particles, n2,
        [epsilon](const vdouble3 &i, const vdouble3 &j) {
          return std::sqrt((j - i).norm() + epsilon);
        },
        [epsilon](const_particle_reference i, const_particle_reference j) {
          return std::sqrt(
              ((get<position>(j) - get<position>(i)).norm() + epsilon));
        });

    /*`

    As before, this operator can be applied to an Eigen vector

    */

    Eigen::VectorXd c_7 = K_h * b;

    /*`

    [endsect]
    */
    //<-
#endif // HAVE_H2LIB
    //->
    /*`

    [section Block Operators]

    It is common that you would like to compose operators in a tiled or block
    format, and Aboria provides a functionality to do this using the [funcref
    Aboria::create_block_operator].

    Let us assume that we wish to compose the two operators `K` and `K_s` from
    before, and want to perform the following combined operator

    $$
    \begin{align}
    e_i &= \sum_j^N d_j K_s(\mathbf{x}\_i,a_i,\mathbf{x}\_j,a_j)  \text{ for }
    i=1...N \\\\ e_{i+N} &= \sum_j^N d_{j+N}
    K(\mathbf{x}\_i,a_i,\mathbf{x}\_j,a_j) \text{ for } i=1...N \end{align}
    $$

    where $e_i$ and $d_j$ are elements of vectors $\mathbf{e}$ and $\mathbf{d}$
    of size $2N$.  Using matrix notation, and using $\mathbf{K}$ and
    $\mathbf{K}_s$ to represent the operators `K` and `K_s`, this is equivalent
    to

    $$
    \mathbf{e} = \begin{pmatrix}
      \mathbf{K}_s   & 0 \\\\\
      0 & \mathbf{K}
    \end{pmatrix} \mathbf{d}
    $$

    We first need operators representing the zero matrices in the upper right
    and lower left corners of the block operator. We create these in Aboria like
    so

    */
    auto Zero = create_zero_operator(particles, particles);
    /*`

    and then we can create the block operator `Full` like so
    */
    auto Full = create_block_operator<2, 2>(K_s, Zero, Zero, K);
    /*`
    Finally we can create vectors `e` and `d` and apply the block operator
    */

    Eigen::VectorXd d(2 * N);
    d.head(N) = Eigen::VectorXd::LinSpaced(N, 0, 1.0);
    d.tail(N) = Eigen::VectorXd::LinSpaced(N, 0, 1.0);

    Eigen::VectorXd e = Full * d;
    //<-
    TS_ASSERT(e.head(N).isApprox(c_3));
    TS_ASSERT(e.tail(N).isApprox(c_1));
    //->

    /*`

    [endsect]


    [section Iterative Solvers]

    The [classref Aboria::MatrixReplacement] class can multiply other Eigen
    vectors, and can be used in Eigen's iterative solvers. Both
    `Eigen::IdentityPreconditioner` and `Eigen::DiagonalPreconditioner`
    preconditioners are supported. Below is an example of how to use Eigen's
    GMRES iterative solver to solve the equation

    $$\mathbf{c} = \mathbf{K} \mathbf{h}$$

    for input vector $\mathbf{h}$.

    We can simply pass the dense operator `K` to Eigen's GMRES iterative solver
    like so
    */

    Eigen::GMRES<decltype(K), Eigen::DiagonalPreconditioner<double>>
        gmres_matrix_free;
    gmres_matrix_free.setMaxIterations(2 * N);
    gmres_matrix_free.set_restart(2 * N + 1);
    gmres_matrix_free.compute(K);
    Eigen::VectorXd h_1 = gmres_matrix_free.solve(c_1);
    /*`
    This will solve the equation in a matrix-free fashion. Alternatively, we can
    use the normal matrix `K_eigen` that we assembled previously to solve the
    equation
    */

    Eigen::GMRES<decltype(K_eigen), Eigen::DiagonalPreconditioner<double>>
        gmres_matrix;
    gmres_matrix.setMaxIterations(2 * N);
    gmres_matrix.set_restart(2 * N + 1);
    gmres_matrix.compute(K_eigen);
    Eigen::VectorXd h_2 = gmres_matrix.solve(c_1);

    //<-
    std::cout << "GMRES (matrix):    #iterations: " << gmres_matrix.iterations()
              << ", estimated error: " << gmres_matrix.error()
              << " actual error: " << (K_eigen * h_2 - c_1).norm() / c_1.norm()
              << std::endl;

    std::cout << "GMRES (matrix-free):    #iterations: "
              << gmres_matrix_free.iterations()
              << ", estimated error: " << gmres_matrix_free.error()
              << " actual error: " << (K * h_1 - c_1).norm() / c_1.norm()
              << std::endl;

    TS_ASSERT(h_1.isApprox(h_2, 1e-11));
    TS_ASSERT(h_1.isApprox(b, 1e-11));
//->

/*`
[endsect]
[endsect]
*/
//]

/*
 *
 * TODO: Will keep this for now until I can get this working, but Eigen doesn't
 * properly support this, so will have to think about it... probable make the
 * interleaved operator a tensor??
[section Interleaved Operators]

Rather than combining multiple operators via blocking, we can also opt to
interleave the operators together. This can often be more efficient in
terms of memory bandwidth. Currently Aboria can interleave operators as
long as they are of a similar type. For example, below will will show how
to interleave two different sparse operators.

We wish to calculate the following operator

$$
\begin{align}
m_{2i-1} &= \sum_j^N n_{2j-1} K^1_s(\mathbf{x}\_i,a_i,\mathbf{x}\_j,a_j)  \text{
for } i=1...N \\\\ m_{2i} &= \sum_j^N n_{2j}
K^2_s(\mathbf{x}\_i,a_i,\mathbf{x}\_j,a_j) \text{ for } i=1...N \end{align}
$$

where $K^1_s$ and $K^2_s$ are given by

$$
K^1_s(\mathbf{x}\_i,a_i,\mathbf{x}\_j,a_j) =
    \begin{cases}
        \frac{a_i a_j}{||\mathbf{dx}\_{ij}||}, & \text{for }
||\mathbf{dx}\_{ij}||<r \\\\\ 0 & \text{otherwise}. \end{cases}
$$
$$
K^2_s(\mathbf{x}\_i,a_i,\mathbf{x}\_j,a_j) =
    \begin{cases}
        \frac{a_i + a_j}{||\mathbf{dx}\_{ij}||}, & \text{for }
||\mathbf{dx}\_{ij}||<r \\\\\ 0 & \text{otherwise}. \end{cases}
$$

We can achieve this in Aboria by creating a sparse operator that returns a
two-dimensional vector type [classref Aboria::double2], instead of a `double`

        auto K_interleaved = create_sparse_operator(particles,particles,
                r,
                [](const_position_reference dx,
                   const_particle_reference i,
                   const_particle_reference j) {
                return double2(
                        (get<a>(i) * get<a>(j)) / dx.norm(),
                        (get<a>(i) + get<a>(j)) / dx.norm()
                        );
                });
Then we can create a new Eigen vector type using [classref Aboria::double2] as
the `Scalar` type, and create two vectors of this type, `m` and `n`.

        typedef Eigen::Matrix<double2,Eigen::Dynamic,1> double2_vector_type;
        double2_vector_type n =
double2_vector_type::LinSpaced(N,double2(0),double2(1.0));

Finally we can apply the interleaved operator to the vector `n`, and store
the result in `m`.


        double2_vector_type m = K_interleaved*n;

        for (int i=0; i<N; ++i) {
            TS_ASSERT_DELTA(m[i][0],c_3[i],std::numeric_limits<double>::epsilon());
        }

[endsect]
*/
#endif // HAVE_EIGEN
  }

  void test_dense_operator(void) {
#ifdef HAVE_EIGEN
    ABORIA_VARIABLE(scalar1, double, "scalar1")
    ABORIA_VARIABLE(scalar2, double, "scalar2")

    typedef Particles<std::tuple<scalar1, scalar2>> ParticlesType;
    typedef position_d<3> position;
    ParticlesType particles;

    double s_init1 = 1.0;
    double s_init2 = 2.0;
    double diameter = 0.1;

    ParticlesType::value_type p;
    get<position>(p) = vdouble3(0, 0, 0);
    get<scalar1>(p) = s_init1;
    get<scalar2>(p) = s_init2;
    particles.push_back(p);
    get<position>(p) = vdouble3(diameter * 0.9, 0, 0);
    particles.push_back(p);
    get<position>(p) = vdouble3(diameter * 1.8, 0, 0);
    particles.push_back(p);

    const size_t n = 3;

    //      3  3  3
    // A =  3  3  3
    //      3  3  3
    auto A = create_dense_operator(
        particles, particles,
        [](ParticlesType::const_reference a, ParticlesType::const_reference b) {
          return get<scalar1>(a) + get<scalar2>(b);
        });

    Eigen::VectorXd v(3);
    v << 1, 2, 3;
    Eigen::VectorXd ans(3);
    ans = A * v;
    for (size_t i = 0; i < n; i++) {
      TS_ASSERT_EQUALS(ans[i], (s_init1 + s_init2) * v.sum());
    }
    v << 0, 2, 1;
    ans = A * v;
    for (size_t i = 0; i < n; i++) {
      TS_ASSERT_EQUALS(ans[i], (s_init1 + s_init2) * v.sum());
    }

    Eigen::MatrixXd A_copy(n, n);
    Eigen::VectorXd ans_copy(n);
    A.assemble(A_copy);
    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        TS_ASSERT_DELTA(A_copy(i, j), A.coeff(i, j),
                        std::numeric_limits<double>::epsilon());
      }
    }

    ans_copy = A_copy * v;

    for (size_t i = 0; i < n; i++) {
      TS_ASSERT_EQUALS(ans[i], ans_copy[i]);
    }

    //        3   3   3
    //       -1  -1  -1
    // A2 =   3   3   3
    //       -1  -1  -1
    //        3   3   3
    typedef Eigen::Matrix<double, 2, 1> eigen_vector;
    auto A2 = create_dense_operator(
        particles, particles,
        [](ParticlesType::const_reference a, ParticlesType::const_reference b) {
          return eigen_vector(get<scalar1>(a) + get<scalar2>(b),
                              get<scalar1>(a) - get<scalar2>(b));
        });

    v << 1, 2, 3;
    ans = A2 * v;
    for (size_t i = 0; i < n; i++) {
      TS_ASSERT_EQUALS(ans[2 * i], (s_init1 + s_init2) * v.sum());
      TS_ASSERT_EQUALS(ans[2 * i + 1], (s_init1 - s_init2) * v.sum());
    }
    v << 0, 2, 1;
    ans = A2 * v;
    for (size_t i = 0; i < n; i++) {
      TS_ASSERT_EQUALS(ans[2 * i], (s_init1 + s_init2) * v.sum());
      TS_ASSERT_EQUALS(ans[2 * i + 1], (s_init1 - s_init2) * v.sum());
    }

    A_copy.resize(2 * n, n);
    A2.assemble(A_copy);
    for (size_t i = 0; i < 2 * n; i++) {
      for (size_t j = 0; j < n; j++) {
        TS_ASSERT_DELTA(A_copy(i, j), A2.coeff(i, j),
                        std::numeric_limits<double>::epsilon());
      }
    }

    ans_copy = A_copy * v;

    for (size_t i = 0; i < 2 * n; i++) {
      TS_ASSERT_EQUALS(ans[i], ans_copy[i]);
    }

#endif // HAVE_EIGEN
  }

  void test_sparse_operator(void) {
#ifdef HAVE_EIGEN
    ABORIA_VARIABLE(scalar1, double, "scalar1")
    ABORIA_VARIABLE(scalar2, double, "scalar2")

    typedef Particles<std::tuple<scalar1, scalar2>> ParticlesType;
    typedef position_d<3> position;
    ParticlesType particles;

    double diameter = 0.1;
    vdouble3 min = vdouble3::Constant(-1);
    vdouble3 max = vdouble3::Constant(1);
    vbool3 periodic = vbool3::Constant(false);

    double s_init1 = 1.0;
    double s_init2 = 2.0;
    ParticlesType::value_type p;
    get<position>(p) = vdouble3(0, 0, 0);
    get<scalar1>(p) = s_init1;
    get<scalar2>(p) = s_init2;
    particles.push_back(p);
    get<position>(p) = vdouble3(diameter * 0.9, 0, 0);
    particles.push_back(p);
    get<position>(p) = vdouble3(diameter * 1.8, 0, 0);
    particles.push_back(p);

    const size_t n = 3;

    particles.init_neighbour_search(min, max, periodic);

    //      3  3  0
    // C =  3  3  3
    //      0  3  3
    auto C = create_sparse_operator(particles, particles, diameter,
                                    [](const position::value_type &dx,
                                       ParticlesType::const_reference a,
                                       ParticlesType::const_reference b) {
                                      return get<scalar1>(a) + get<scalar2>(b);
                                    });

    Eigen::VectorXd v(n);
    v << 1, 2, 3;

    // 3  3  0   1   9
    // 3  3  3 * 2 = 18
    // 0  3  3   3   15
    Eigen::VectorXd ans = C * v;

    for (size_t i = 0; i < n; i++) {
      double sum = 0;
      for (size_t j = 0; j < n; j++) {
        if ((get<id>(particles[i]) == 0) && (get<id>(particles[j]) == 2)) {
          sum += 0;
        } else if ((get<id>(particles[i]) == 2) &&
                   (get<id>(particles[j]) == 0)) {
          sum += 0;
        } else {
          sum += (s_init1 + s_init2) * v[j];
        }
      }
      TS_ASSERT_EQUALS(ans[i], sum);
    }

    Eigen::MatrixXd C_copy(n, n);
    Eigen::VectorXd ans_copy(n);
    C.assemble(C_copy);
    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        TS_ASSERT_DELTA(C_copy(i, j), C.coeff(i, j),
                        std::numeric_limits<double>::epsilon());
      }
    }

    ans_copy = C_copy * v;

    for (size_t i = 0; i < n; i++) {
      TS_ASSERT_EQUALS(ans[i], ans_copy[i]);
    }

    Eigen::SparseMatrix<double> C_sparse(n, n);
    C.assemble(C_sparse);
    TS_ASSERT_EQUALS(C_sparse.nonZeros(), 7);
    for (int k = 0; k < C_sparse.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(C_sparse, k); it;
           ++it) {
        TS_ASSERT_EQUALS(it.value(), C.coeff(it.row(), it.col()));
      }
    }

    ans_copy = C_sparse * v;

    for (size_t i = 0; i < n; i++) {
      TS_ASSERT_EQUALS(ans[i], ans_copy[i]);
    }

    //       3  3  0
    //      -1 -1  0
    // C2 =  3  3  3
    //      -1 -1 -1
    //       0  3  3
    //       0 -1 -1
    typedef Eigen::Matrix<double, 2, 1> eigen_vector;
    auto C2 = create_sparse_operator(
        particles, particles, diameter,
        [](const position::value_type &dx, ParticlesType::const_reference a,
           ParticlesType::const_reference b) {
          return eigen_vector(get<scalar1>(a) + get<scalar2>(b),
                              get<scalar1>(a) - get<scalar2>(b));
        });

    v << 1, 2, 3;

    //  3  3  0   1    9
    // -1 -1  0 * 2 = -3
    //  3  3  3   3   18
    // -1 -1 -1       -7
    //  0  3  3       15
    //  0 -1 -1       -5
    ans = C2 * v;

    std::vector<double> check = {9, -3, 18, -7, 15, -5};
    for (size_t i = 0; i < n; i++) {
      TS_ASSERT_EQUALS(ans[i], check[i]);
    }

    C_copy.resize(2 * n, n);
    C2.assemble(C_copy);
    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        TS_ASSERT_DELTA(C_copy(i, j), C2.coeff(i, j),
                        std::numeric_limits<double>::epsilon());
      }
    }

    ans_copy = C_copy * v;

    for (size_t i = 0; i < n; i++) {
      TS_ASSERT_EQUALS(ans[i], ans_copy[i]);
    }

    C_sparse.resize(2 * n, n);
    C2.assemble(C_sparse);
    TS_ASSERT_EQUALS(C_sparse.nonZeros(), 14);
    for (int k = 0; k < C_sparse.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(C_sparse, k); it;
           ++it) {
        TS_ASSERT_EQUALS(it.value(), C2.coeff(it.row(), it.col()));
      }
    }

    ans_copy = C_sparse * v;

    for (size_t i = 0; i < n; i++) {
      TS_ASSERT_EQUALS(ans[i], ans_copy[i]);
    }

#endif // HAVE_EIGEN
  }

  void test_block_operator(void) {
#ifdef HAVE_EIGEN
    ABORIA_VARIABLE(scalar1, double, "scalar1")
    ABORIA_VARIABLE(scalar2, double, "scalar2")

    typedef Particles<std::tuple<scalar1, scalar2>> ParticlesType;
    typedef position_d<3> position;
    ParticlesType particles, augment;

    double diameter = 0.1;
    vdouble3 min = vdouble3::Constant(-1);
    vdouble3 max = vdouble3::Constant(1);
    vbool3 periodic = vbool3::Constant(false);

    ParticlesType::value_type p;
    get<position>(p) = vdouble3(0, 0, 0);
    get<scalar1>(p) = 1;
    get<scalar2>(p) = 0.1;
    particles.push_back(p);
    get<position>(p) = vdouble3(diameter * 0.9, 0, 0);
    get<scalar1>(p) = 2;
    get<scalar2>(p) = 0.2;
    particles.push_back(p);
    get<position>(p) = vdouble3(diameter * 1.8, 0, 0);
    get<scalar1>(p) = 3;
    get<scalar2>(p) = 0.3;
    particles.push_back(p);
    get<scalar1>(p) = 0;
    get<scalar2>(p) = 0;
    get<position>(p) = vdouble3(-diameter * 1.8, 0, 0);
    augment.push_back(p);

    const size_t n = 3;

    particles.init_neighbour_search(min, max, periodic);

    //      1  1  1
    // A =  2  2  2
    //      3  3  3
    auto A = create_dense_operator(
        particles, particles,
        [](ParticlesType::const_reference a, ParticlesType::const_reference b) {
          return get<scalar1>(a);
        });

    Eigen::VectorXd v(n);
    v << 1, 1, 1;
    Eigen::VectorXd ans(n);
    ans = A * v;
    TS_ASSERT_EQUALS(ans[0], 3);
    TS_ASSERT_EQUALS(ans[1], 6);
    TS_ASSERT_EQUALS(ans[2], 9);

    // B = 0.1
    //     0.2
    //     0.3
    auto B = create_dense_operator(
        particles, augment,
        [](ParticlesType::const_reference a, ParticlesType::const_reference b) {
          return get<scalar2>(a);
        });

    // C = 0.1 0.2 0.3
    auto C = create_dense_operator(
        augment, particles,
        [](ParticlesType::const_reference a, ParticlesType::const_reference b) {
          return get<scalar2>(b);
        });

    auto Zero = create_zero_operator(augment, augment);

    //         1   1   1   0.1
    //         2   2   2   0.2
    // Full =  3   3   3   0.3
    //         0.1 0.2 0.3 0
    auto Full = create_block_operator<2, 2>(A, B, C, Zero);
    std::cout << Full.rows() << " " << Full.cols() << std::endl;

    v.resize(n + 1);
    v << 1, 1, 1, 1;
    ans.resize(n + 1);

    //         1   1   1   0.1   1   3.1
    //         2   2   2   0.2   1   6.2
    // ans  =  3   3   3   0.3 * 1 = 9.3
    //         0.1 0.2 0.3 0     1   0.6
    ans = Full * v;
    TS_ASSERT_DELTA(ans[0], 3.1, std::numeric_limits<double>::epsilon());
    TS_ASSERT_DELTA(ans[1], 6.2, std::numeric_limits<double>::epsilon());
    TS_ASSERT_DELTA(ans[2], 9.3, std::numeric_limits<double>::epsilon());
    TS_ASSERT_DELTA(ans[3], 0.6, std::numeric_limits<double>::epsilon());

    Eigen::MatrixXd Full_copy(n + 1, n + 1);
    Eigen::VectorXd ans_copy(n + 1);
    Full.assemble(Full_copy);
    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        TS_ASSERT_EQUALS(Full_copy(i, j), Full.coeff(i, j));
      }
    }

    ans_copy = Full_copy * v;

    for (size_t i = 0; i < n; i++) {
      TS_ASSERT_EQUALS(ans[i], ans_copy[i]);
    }

    Eigen::SparseMatrix<double> Full_sparse(n + 1, n + 1);
    Full.assemble(Full_sparse);
    TS_ASSERT_EQUALS(Full_sparse.nonZeros(), 15);
    for (int k = 0; k < Full_sparse.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(Full_sparse, k); it;
           ++it) {
        TS_ASSERT_EQUALS(it.value(), Full.coeff(it.row(), it.col()));
      }
    }

    ans_copy = Full_sparse * v;

    for (size_t i = 0; i < n; i++) {
      TS_ASSERT_EQUALS(ans[i], ans_copy[i]);
    }

#endif // HAVE_EIGEN
  }
};

#endif /* OPERATORSTEST_H_ */
