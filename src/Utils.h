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
void radial_distribution_function(Particles<T>& particles,
							const double min, const double max,
							std::vector<double>& out) {


}

template<typename T>
void copy (Particles<T>* from, Particles<T>* to) {
	to->clear();
	to->generate(from->size(),[from](Particles<T>::Value& i) {
		i.deep_copy(from[i.get_index()]);
	});
}


#endif /* UTILS_H_ */
