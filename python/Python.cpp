/*
 * Python.cpp
 *
 *  Created on: 7 Sep 2014
 *      Author: mrobins
 */


#include "Python.h"

typedef std::tuple<
		{{for loop, i in looper(types)}}
			{{i}} {{if not loop.last}},{{endif}}
		{{endfor}}
		> SpeciesTuple;
typedef Particles<SpeciesTuple> SpeciesType;

namespace Aboria {
template <>
class DataNames<SpeciesTuple> {
public:
	std::string static get(unsigned int n) {
		{{py: n = len(names)}}
		std::string ret[{{n}}] = {
				{{for loop, i in looper(names)}}
					{{i}} {{if not loop.last}},{{endif}}
	{			{endfor}}};
		return ret[n];
	}
};
}


BOOST_PYTHON_MODULE({{name}}) {

	VTK_PYTHON_CONVERSION(vtkUnstructuredGrid);

	Vect3_from_python_list<double>();
	to_python_converter<
		Vect3d,
		Vect3_to_python<double> >();

	/*
	 * Particles
	 */
	class_<SpeciesType,typename Aboria::ptr<SpeciesType> >("Particles")
	        .def(boost::python::vector_indexing_suite<SpeciesType >())
	        .def("get_grid",&SpeciesType::get_grid)
	        .def("copy_from_vtk_grid",&SpeciesType::copy_from_vtk_grid)
	    ;


	/*
	 * Particle
	 */
#define ADD_PROPERTY(NUMBER) \
		.add_property(DataNames<SpeciesTuple>::get(NUMBER).c_str(), \
					make_function(&SpeciesType::value_type::get_data_elem<NUMBER>, \
									return_value_policy<copy_non_const_reference>()), \
					&SpeciesType::value_type::set_data_elem<NUMBER>)

	class_<SpeciesType::value_type, Aboria::ptr<SpeciesType::value_type> >("Particle",init<>())
		.add_property("id", &SpeciesType::value_type::get_id)
		.add_property("position",  make_function(&SpeciesType::value_type::get_position,
									return_value_policy<copy_const_reference>()),
							&SpeciesType::value_type::set_position)
		{{py:index_list = range(n)}}
		{{for a in index_list}}
			ADD_PROPERTY(a)
		{{endfor}}
		;

}
