Aboria
=====

Aboria implements a STL container of particles in 3D space. The library is header-only.
The container supports random access of particles, as well as the normal STL algorithms.
Neighbourhood searches are possible, using a bucket search method (uniform bucket spacing)

This code is in alpha, and the API is currently in flux, use at your own risk. There is an example
of how to use the code in the example subdirectory, which implements a simple DEM algorithm, with 
linear spring contact forces. You will need cmake and VTK installed to compile the example.