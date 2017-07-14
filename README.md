[![TravisCI](https://travis-ci.org/martinjrobins/Aboria.svg?branch=master)](https://travis-ci.org/martinjrobins/Aboria)
[![Coverage](https://coveralls.io/repos/martinjrobins/Aboria/badge.svg?branch=master&service=github)](https://coveralls.io/github/martinjrobins/Aboria?branch=master)
<!---
[![AppVeyor](https://ci.appveyor.com/api/projects/status/6aimud6e8tvxfwgm?svg=true)](https://ci.appveyor.com/project/martinjrobins/aboria)
-->

Aboria implements an Embedded Domain Specific Language (eDSL) in C++ for 
specifying expressions over particles and their neighbours in N dimensional 
space, with the aim of providing a useful library for implementing 
particle-based numerical algorithms, for example Molecular Dynamics, Smoothed 
Particle Hydrodynamics or Radial Basis Functions. 

Aboria provides:
* a STL compatible container class to store a particle set containing
  a position and unique id for each particle, as well as any number of 
  user-defined variables with arbitrary types.
* the ability to embed each particle set within a hypercube N-dimensional
  domain with arbitrary periodicity
* flexible neighbourhood queries that return iterators, and can use any 
  L-norm distance measure (L1 = Manhattan distance, L2 = Euclidean distance, 
  ... , Linf = Chebyshev distance)
* an expression template API for forming non-linear operators over the 
  particles. This can be used, for example, to implement interaction forces
  in Molecular Dynamics.
* an API for forming linear kernel operators from C++ lambda functions, This
  can be used, for example, to implement Radial Basis Function kernels. These 
  can be wrapped as [Eigen](eigen.tuxfamily.org) matrices in order to solve 
  linear systems based on kernel operators.
    
    
Aboria is distributed under a BSD 3-Clause License, see LICENCE for more 
details. For documentation see the [Aboria 
website](https://martinjrobins.github.io/Aboria). If you are interested in 
contributing to Aboria, having trouble getting it working or just have a 
question, send me an email at <martin.robinson@cs.ox.ac.uk> or create a
GitHub issue or pull request.

