/*
 * particle_for_each.cpp
 *
 *  Created on: 29 Jan 2014
 *      Author: mrobins
 */
#include <Aboria.h>

int main(int argc, char **argv) {
	using namespace Aboria;

	Particles<1> p;

	p.add_particle(Vect3d(0,0,0));
	p.add_particle(Vect3d(1,0,0));
	p.add_particle(Vect3d(0,1,0));
	p.add_particle(Vect3d(0,0,1));
	p.add_particle(Vect3d(2,0,0));
	p.add_particle(Vect3d(0,2,0));
	p.add_particle(Vect3d(0,0,2));

	auto operation = for_each(p,[](Particles<1>::Value i) {std::cout << i.get_position() << std::endl;});

	operation->execute();
	int dt = 2;
	p.init_neighbour_search(Vect3d(0,0,0),Vect3d(2,2,2),1.1);
	auto operator2 = for_each(p,[&](Particles<1>::Value i) {
		std::cout << "looking for particles within radius 1.1 of "<<i.get_position() << std::endl;
		for (auto j: i.get_in_radius(p,1.1)) {
			std::cout << "found particles at "<<j<<std::endl;
		}

	});

	operator2->execute();

}

