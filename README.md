Aboria
=====

Aboria implements a STL container of particles in 3D space. The library is header-only.
The container supports random access of particles, as well as the normal STL algorithms.
Neighbourhood searches are possible, using a bucket search method (uniform bucket spacing).

The motivation behind Aboria is to provide a useful library for implementing particle-based numerical algorithms, for example Smoothed Particle Hydrodynamics or Molecular Dynamics. Each particle has a 3D position and user-defined data-package (for other variables such as velocity, density etc) and is optionally embedded within a cuboidal spatial domain (for neighbourhood searches) that can be periodic or not. Each particle also has its own random number generator that is seeded via its own unique id.

Examples
--------

The *examples/* subdirectory contains a collection of examples for using Aboria. Currently these are:

- *examples/sph* - An Smoothed Particle Hydrodynamics example, simulating a 2D water column over a no-slip boundary. The *x* and *y* directions are periodic.
- *examples/dem* - An Discrete Element Model example, simulating 2 spherical particles falling onto an surface.
- *exampes/sphdem* - A coupled SPH and DEM example, simulating a single DEM particle falling down a water column
- *examples/bd* - Brownian dynamics of N particles within a reflecting sphere


A short sample from the DEM example, which shows what is possible with the library. This shows a `for_each`
loop over the DEM particles, calculating the contact forces between pairs of particles

```
std::for_each(dem->begin(),dem->end(),[&geometry,dem,dem_k,dem_gamma,dem_mass,dem_diameter](DemType::Value& i) {
		const Vect3d& r = i.get_position();
		Vect3d& f = std::get<DEM_FORCE>(i.get_data());
		Vect3d& v = std::get<DEM_VELOCITY>(i.get_data());

		f << 0,0,0;
		f = f + geometry(i);

		for (auto tpl: i.get_neighbours(dem)) {
			const Vect3d& dx = std::get<1>(tpl);
			const DemType::Value& j = std::get<0>(tpl);
			const Vect3d& vj = std::get<DEM_VELOCITY>(j.get_data());
			if (i.get_id()==j.get_id()) continue;

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


Creating New Particles
----------------------

The main particles data-structure is called `Particles`. It takes one template arguement, which is the type of the data package given to each particle. This type is restricted to being a tuple. So, for example, the following creates a set of particles which each have (along with the standard variables such as position, id etc) a data package consisting one one `double` variable.

```
using namespace Aboria;

typedef Particles<std::tuple<double> > MyParticles;
MyParticles particles1();
```

You can also give the `MyParticles` constructor a single `unsigned int` arguement to set the random seed for the container:

```
MyParticles particles2(0);
```

To create new particles simply use the `value_type` of the container type. Each particle constructor takes a single `Vect3d` type for the particle position.

```
typedef MyParticles::value_type MyParticle;
particles.push_back(MyParticle(Vect3d(0,0,0)));
particles.push_back(MyParticle(Vect3d(1,0,0)));
```
