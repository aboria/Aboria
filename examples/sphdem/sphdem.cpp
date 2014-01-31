/*
 * sphdem.cpp
 *
 *  Created on: 30 Jan 2014
 *      Author: mrobins
 */

#include "sphdem.h"

int main(int argc, char **argv) {

	auto sph = SphType::New();
	auto dem = DemType::New();
	auto params = ptr<ParamTuple>(new ParamTuple());

	params.get<PARAMS_DEM_DT>() = 0.1;
	params.get<PARAMS_DEM_DIAMETER>() = 1;
	params.get<PARAMS_DEM_GAMMA>() = 0;
	params.get<PARAMS_DEM_K>() = 1;
	params.get<PARAMS_DEM_MASS>() = 1;
	const int n = 1000;
	const int ndem = 100;

	auto geometry = [](Vect3d& r) {
		return r[0];
	};

	dem.create_particles(ndem,[](DemType::Value& i) {
		return Vect3d(i.get_id(),0,ndem);
	});

	for (int i = 0; i < n; ++i) {
		dem_start(dem,params,geometry);
		dem_end(dem,params,geometry);
	}


}
