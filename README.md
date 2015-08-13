Aboria
=====

[![TravisCI](https://travis-ci.org/martinjrobins/Aboria.svg?branch=master)](https://travis-ci.org/martinjrobins/Aboria)
<!---
[![AppVeyor](https://ci.appveyor.com/api/projects/status/6aimud6e8tvxfwgm?svg=true)](https://ci.appveyor.com/project/martinjrobins/aboria)
[![Coverity](https://scan.coverity.com/projects/5951/badge.svg)](https://scan.coverity.com/projects/5951)
[![Coverage](https://coveralls.io/repos/martinjrobins/Aboria/badge.svg?branch=master&service=github)](https://coveralls.io/github/martinjrobins/Aboria?branch=master)
-->

Aboria implements a STL container of particles in 3D space. The library is header-only.
The container supports random access of particles, as well as the normal STL algorithms.
Neighbourhood searches are possible, using a bucket search method (uniform bucket spacing).

The motivation behind Aboria is to provide a useful library for implementing particle-based numerical algorithms, for example Smoothed Particle Hydrodynamics or Molecular Dynamics. Each particle has a 3D position and user-defined data-package (for other variables such as velocity, density etc) and is optionally embedded within a cuboidal spatial domain (for neighbourhood searches) that can be periodic or not. Each particle also has its own random number generator that is seeded via its own unique id.

- [Creating Particles](#create)
- [Particle Objects](#particle)
- [Looping through a container](#looping)
- [Neighbourhood Searching](#neighbour)

Examples
--------

The *examples/* subdirectory contains a collection of examples for using Aboria. Currently these are:

- *examples/sph* - An Smoothed Particle Hydrodynamics example, simulating a 2D water column over a no-slip boundary. The *x* and *y* directions are periodic.
- *examples/dem* - An Discrete Element Model example, simulating 2 spherical particles falling onto an surface.
- *examples/dem_symbolic* - An Discrete Element Model example using the symbolic interface, simulating a polydisperse set of spherical particles falling onto an surface.
- *exampes/sphdem* - A coupled SPH and DEM example, simulating a single DEM particle falling down a water column
- *examples/bd* - Brownian dynamics of N particles within a reflecting sphere
- *examples/bd_symbolic* - Brownian dynamics of N point particles around a set of spheres, using the symbolic interface. The point particles reflect off the spheres as they diffuse.


A short sample from the DEM example, which shows what is possible with the library. This shows a `for_each`
loop over the DEM particles, calculating the contact forces between pairs of particles

```Cpp
std::for_each(dem->begin(),dem->end(),[&geometry,dem,dem_k,dem_gamma,dem_mass,dem_diameter](DemType::value_type& i) {
		const Vect3d& r = get<position>(i);
		Vect3d& f = get<force>(i);
		Vect3d& v = get<velocity>(i);

		f << 0,0,0;
		f = f + geometry(i);

		for (auto tpl: i.get_neighbours(dem)) {
			const Vect3d& dx = std::get<1>(tpl);
			const DemType::value_type& j = std::get<0>(tpl);
			const Vect3d& vj = get<velocity>(j);
			if (get<id>(i)==get<id>(j)) continue;

			const double r = dx.norm();
			const double overlap = dem_diameter-r;
			if (overlap>0) {
				const Vect3d normal = dx/r;
				const Vect3d dv = v-vj;
				const double overlap_dot = -dv.dot(normal);
				f += (dem_k*overlap + dem_gamma*overlap_dot)*normal/dem_mass;
			}
		}

	});
```

This can be further simplified to one line of code by using the symbolic interface, which provides a succinct way to specify accumulation loops over neighbours and expressions using particle variables.

```Cpp
dvdt = (// spring force between dem particles
        sum(b=dem, id_[a]!=id_[b] && norm_(dx)<r[a]+r[b], 
                   dem_k*((r[a]+r[b])/norm_(dx)-1)*dx  + dem_gamma*(v[b]-v[a]))
                
        // spring force between particles and bottom wall
        + if_else(r-p[2] > 0, dem_k*(r-p[2]), 0.0)*Vect3d(0,0,1) 

        // gravity force
        + Vect3d(0,0,-9.81)*m

       )/m;
```


<a name="create">Creating New Particles</a>
-------------------------------------------

The main particles data-structure, or container, is called `Particles`. It is templated on zero or more variable types. For example, the following creates a set of particles which each have (along with the standard variables such as position, id etc) a data package consisting of one `double` variable type named `scalar`.

```Cpp
using namespace Aboria;

ABORIA_VARIABLE(scalar,double,"my scalar")
typedef Particles<scalar> MyParticles;
MyParticles particles();
```

If you wanted each particle to have a `potential` variable held as a `double`, as well as a `velocity` variable held as a `Vect3d` vector class, then you would write the following

```Cpp
ABORIA_VARIABLE(potential,double,"potential energy")
ABORIA_VARIABLE(velocity,Vect3d,"velocity")
typedef Particles<potential,velocity> MyParticles;
MyParticles particles();
```

You can also give the `MyParticles` constructor a single `unsigned int` argument to set the random seed for the container:

```Cpp
MyParticles particles_with_seed(0);
```

To create new particles simply use the `value_type` of the container type. Each particle constructor takes a single `Vect3d` type for the particle position.

```Cpp
typedef MyParticles::value_type MyParticle;
particles.push_back(MyParticle(Vect3d(0,0,0)));
particles.push_back(MyParticle(Vect3d(1,0,0)));
```

<a name="particle">Particle Objects</a>
---------------------------------------

The `value_type` of the `Particles` container is a data-structure representing each particle. By default each particle has a position, a unique id and a boolean flag indicating if this particle is active or not. Use `get<position>()` to access the position, `get<id>()` for the id and `get<alive>()` for the alive flag.

```Cpp
MyParticle& particle = particles[0];
std::cout <<"Position = "<<get<position>(particle) << 
   ". Id = "<<get<id>(particle)<< ". Particle is ";
if (get<alive>(particle)) {
   std::cout << "alive. " << "\n";
} else {
   std::cout << "dead. " << "\n";
}
```

You can access the data by templating the `get` function with the variable type, for example

```Cpp
std::cout << "The scalar data is " << get<scalar>(particle) << "\n";
```

<a name="looping">Looping through the container</a>
---------------------------------------------------

You can use the indexing operator `Operator[]` to simply loop through the container

```Cpp
for (int i=0; i < particles.size(); i++) {
   std::cout << "Accessing particle with id = " << get<id>(particles[i]) << "\n";
}
```

Or you can use the normal STL `begin()` and `end()` functions that return random access iterators to the beginning and end of the container.

```Cpp
for (auto i = particles.begin(); i != particles.end(); i++) {
   std::cout << "Accessing particle with id = " << get<id>(*i) << "\n";
}
```

Or

```Cpp
for (auto i: particles) {
   std::cout << "Accessing particle with id = " << get<id>(i) << "\n";
}
```

Or you can use the STL algorithm `for_each`. If you are using a GCC compiler, you can turn on the parallel mode to enable this loop to be run in parallel

```Cpp
std::for_each(particles.begin(), particles.end(), [](MyParticle& i) {
   std::cout << "Accessing particle with id = " << get<id>(i) << "\n";
});
```

<a name="neighbour">Neighbourhood Searching</a>
-----------------------------------------------

The `Particles` container gives you neighbourhood searching functionality, using a simple Cell List or Linked-List approach. The domain is divided into a regular grid of cubes with side length equal to a constant lengthscale that is supplied by the user. Each particle in the container is assigned to the cell that contains its position. Neighbourhood search queries at a given point return all the particles within the cell that contains this point and the immediate cell neighbours.

Before you can use the neighbourhood searching, you need to initialise the domain using the `init_neighbour_search` function

```Cpp
Vect3d min(-1,-1,-1);
Vect3d max(1,1,1);
Vect3b periodic(true,true,true);
double diameter = 0.1;
particles.init_neighbour_search(min,max,diameter,periodic);
```

Here `diameter` is the lengthscale of the neighbourhood search. That is, any particles that are separated by more than `diameter` might not be classified as neighbours.

Once this is done you can begin using the neighbourhood search queries using the `get_neighbours` function. This returns a lightweight container with `begin()` and `end()` functions that return `const` forward only iterators to the particles that satisfy the neighbour search. For example, the following counts all the particles within a square domain of side length `diameter` of the point (0,0,0)

```Cpp
int count = 0;
for (auto i: particles.get_neighbours(Vect3d(0,0,0))) {
   count++;
}
std::cout << "There are "<< count << " particles.\n";
```

When dereferenced, the neighbourhood iterator returns a tuple of size 2 containing 

1. The found particle object
2. $dx$, a vector pointing to the found point from the query point. I.e. if $x_a$ is the query point and $x_b$ is the found point, then $dx = x_a - x_b$.

The latter is useful for periodic domains, the returned vector $dx$ takes periodic domains into account and returns the $dx$ with the smallest length. 

For example, 

```Cpp
for (auto i: particles.get_neighbours(Vect3d(0,0,0))) {
   const MyParticle& b = std::get<0>(tpl);
   const Vect3d& dx = std::get<1>(tpl);
   std::cout << "Found a particle with dx = " << dx << " and id = " << get<id>(b) << "\n";
}
```
