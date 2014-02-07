Aboria
=====

Aboria implements a STL container of particles in 3D space. The library is header-only.
The container supports random access of particles, as well as the normal STL algorithms.
Neighbourhood searches are possible, using a bucket search method (uniform bucket spacing)

This code is in alpha, and the API is currently in flux, use at your own risk. There is an example
of how to use the code in the example subdirectory, which implements a simple DEM algorithm, with 
linear spring contact forces. You will need cmake and VTK installed to compile the example.

A short sample from the DEM example, which shows what is possible with the library. This shows a for_each 
loop over the DEM particles, calculating the contact forces between pairs of particles

```
std::for_each(dem->begin(),dem->end(),[&geometry,dem,dem_k,dem_gamma,dem_mass,dem_diameter](DemType::Value& i) {
		const Vect3d& r = i.get_position();
		Vect3d& f = std::get<DEM_FORCE>(i.get_data());
		Vect3d& v = std::get<DEM_VELOCITY>(i.get_data());

		f << 0,0,0;
		f = f + geometry(i);

		for (auto tpl: i.get_in_radius(dem,dem_diameter)) {
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