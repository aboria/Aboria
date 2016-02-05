[![TravisCI](https://travis-ci.org/martinjrobins/Aboria.svg?branch=master)](https://travis-ci.org/martinjrobins/Aboria)
<!---
[![AppVeyor](https://ci.appveyor.com/api/projects/status/6aimud6e8tvxfwgm?svg=true)](https://ci.appveyor.com/project/martinjrobins/aboria)
[![Coverity](https://scan.coverity.com/projects/5951/badge.svg)](https://scan.coverity.com/projects/5951)
[![Coverage](https://coveralls.io/repos/martinjrobins/Aboria/badge.svg?branch=master&service=github)](https://coveralls.io/github/martinjrobins/Aboria?branch=master)
-->

Aboria implements a STL container of particles in 3D space. The library is header-only.
The container supports random access of particles, as well as the normal STL algorithms.
Neighbourhood searches are possible, using a bucket search method (uniform bucket spacing).

Aboria is distributed under a BSD 3-Clause License, see LICENCE for more details

The motivation behind Aboria is to provide a useful library for implementing 
particle-based numerical algorithms, for example Smoothed Particle Hydrodynamics 
or Molecular Dynamics. Each particle has a 3D position and user-defined 
data-package (for other variables such as velocity, density etc) and is 
optionally embedded within a cuboidal spatial domain (for neighbourhood 
searches) that can be periodic or not. Each particle also has its own random 
number generator that is seeded via its own unique id.
