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

template<typename PTYPE, typename GTYPE, typename RNDGTYPE>
void create_n_particles_with_rejection_sampling(const unsigned int n, PTYPE &particles, const GTYPE &generator_fun, const Vect3d &low, const Vect3d &high, RNDGTYPE &generator) {

    std::uniform_real_distribution<double> uni(0,1);
    const Vect3d h = (high-low)/1000;
    double max_generator_fun = 0;
    assert(h[0]>0);
    for (double x = low[0]+h[0]/2; x<high[0]; x += h[0]) {
        double yinc = h[1];
        if (yinc == 0) yinc = h[0];
        for (double y = low[1]+h[1]/2; y<=high[1]; y += yinc) {
            double zinc = h[2];
            if (zinc == 0) zinc = h[0];
            for (double z = low[2]+h[2]/2; z<=high[2]; z += zinc) {
                if (generator_fun(Vect3d(x,y,z)) > max_generator_fun) max_generator_fun = generator_fun(Vect3d(x,y,z));
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        Vect3d r = Vect3d(uni(generator),uni(generator),uni(generator))*(high-low) + low;
        while (uni(generator) > generator_fun(r)/max_generator_fun) {
            r = Vect3d(uni(generator),uni(generator),uni(generator))*(high-low) + low;
        }
        particles.push_back(r);
    }
}



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
			const int index = std::floor((r-min)/bsep);
			if ((index>=0)&&(index<n)) (*out)[index] += 1.0;
		}
	}

	const double volume = (particles->get_high()-particles->get_low()).prod();
	const double rho = particles->size()/volume;

	for (int i = 0; i < n; ++i) {
		(*out)[i] /= particles->size()*(4.0/3.0)*PI*(std::pow((i+1)*bsep+min,3)-std::pow(i*bsep+min,3))*rho;
	}

	particles->reset_neighbour_search(old_size);

	return out;
}

}

#endif /* UTILS_H_ */
