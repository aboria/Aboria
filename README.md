[![TravisCI](https://travis-ci.org/martinjrobins/Aboria.svg?branch=master)](https://travis-ci.org/martinjrobins/Aboria)
[![Coverage](https://coveralls.io/repos/martinjrobins/Aboria/badge.svg?branch=master&service=github)](https://coveralls.io/github/martinjrobins/Aboria?branch=master)
<!---
[![AppVeyor](https://ci.appveyor.com/api/projects/status/6aimud6e8tvxfwgm?svg=true)](https://ci.appveyor.com/project/martinjrobins/aboria)
-->

UPDATE (12/07/2016): The next release of Aboria is currently being documented 
and tested in the `version2` branch. It will be merged in the next few weeks. 
This release:
* reworks the internal storage of the container class to model a set of zipped 
  vectors (currently `std::vector`s)
* adds more compile-time checking of expression correctness
* updates the bucket-search neighbourhood searching algorithm to use STL 
  algorithms only
* adds matrix-free linear algebra capabilities via Eigen <eigen.tuxfamily.org>

Aboria implements an expressive Domain Specific Language (DSL) in C++ for 
specifying expressions over particles and their neighbours in 3D space. The 
library is header-only and based on expression templates for efficient and 
natural expression of mathematical expressions applied over the particles.

The particle data is contained in a STL compatible container. Each particle has 
a 3D position and user-defined data-package (for other variables such as 
velocity, density etc) and is optionally embedded within a cuboidal spatial 
domain (for neighbourhood searches) that can be periodic or not. Users can 
implement standard C++ and STL algorithms directly on the particle data, or use 
the expression template API to naturally form operations over the particle set.

The motivation behind Aboria is to provide a useful library for implementing 
particle-based numerical algorithms, for example Smoothed Particle Hydrodynamics 
or Molecular/Langevin Dynamics.

Aboria is distributed under a BSD 3-Clause License, see LICENCE for more 
details.

For documentation see the [Aboria 
website](https://martinjrobins.github.io/Aboria). If you are interested in 
contributing to Aboria, having trouble getting it working or just have a 
question, send me an email at <martin.robinson@cs.ox.ac.uk>.

