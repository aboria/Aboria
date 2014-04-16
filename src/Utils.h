/*
 * Utils.h
 *
 *  Created on: 15 Apr 2014
 *      Author: mrobins
 */

#ifndef UTILS_H_
#define UTILS_H_

#include "Aboria.h"

namespace Aboria {

template<typename T>
ptr<std::vector<double> > radial_distribution_function(ptr<Particles<T> > particles,
							const double min, const double max,
							const int n) {
	auto out = ptr<std::vector<double> >(new std::vector<double>());

	out->resize(n,0);
	const double bsep = (max-min)/n;
	const double old_size = particles->get_lengthscale();
	particles->reset_neighbour_search(max);

	for(typename Particles<T>::Value& i: *particles) {
		for (auto tpl: i.get_neighbours(particles)) {
			const Vect3d& dx = std::get<1>(tpl);
			const typename Particles<T>::Value& j = std::get<0>(tpl);
			if (i.get_id()==j.get_id()) continue;
			const double r = dx.norm();
			const int index = (r-min)/bsep;
			if ((index>=0)&&(index<n)) (*out)[index] += 1.0;
		}
	}

	const double volume = (particles->get_high()-particles->get_low()).prod();
	const double rho = particles->size()/volume;

	for (int i = 0; i < n; ++i) {
		(*out)[i] /= particles->size()*(4.0/3.0)*PI*(pow((i+1)*bsep+min,3)-pow(i*bsep+min,3))*rho;
	}

	particles->reset_neighbour_search(old_size);

	return out;
}

}

#endif /* UTILS_H_ */
