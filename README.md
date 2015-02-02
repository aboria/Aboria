Aboria
=====

Aboria implements a STL container of particles in 3D space. The library is header-only.
The container supports random access of particles, as well as the normal STL algorithms.
Neighbourhood searches are possible, using a bucket search method (uniform bucket spacing).

The motivation behind Aboria is to provide a useful library for implementing particle-based numerical algorithms, for example Smoothed Particle Hydrodynamics or Molecular Dynamics. Each particle has a 3D position and user-defined data-package (for other variables such as velocity, density etc) and is optionally embedded within a cuboidal spatial domain (for neighbourhood searches) that can be periodic or not. Each particle also has its own random number generator that is seeded via its own unique id.

Examples
--------

The *examples/* subdirectory contains a collection of examples for using Aboria. Currently these are:

- *examples/sph* - An Smoothed Particle Hydrodynamics example, simulating a 3D water column over a no-slip boundary. The *x* and *y* directions are periodic.
- *examples/dem* - An Discrete Element Model example, simulating 2 spherical particles falling onto an surface.
- *exampes/sphdem* - A coupled SPH and DEM example, simulating a single DEM particle falling down a water column
- *examples/bd* - Brownian dynamics of N particles within a refelcting sphere


A short sample from the DEM example, which shows what is possible with the library. This shows a for_each 
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
