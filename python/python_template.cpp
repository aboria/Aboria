/*
 * Python.cpp
 *
 *  Created on: 7 Sep 2014
 *      Author: mrobins
 */


#include "Python.h"

/*{% for variable in variables %}*/
ABORIA_VARIABLE(/*{{variable.type_name}}*/,/*{{variable.representation}}*/,"/*{{variable.description}}*/")
/*{% endfor %}*/

typedef /*{{variable_type_string}}*/ variable_tuple;
typedef Particles<variable_tuple> particles_type;

BOOST_PYTHON_MODULE(/*{{particles.name}}*/) {

	VTK_PYTHON_CONVERSION(vtkUnstructuredGrid);

	Vect3_from_python_list<double>();
	to_python_converter<
		Vect3d,
		Vect3_to_python<double> >();

	/*
	 * Particles
	 */
	class_<particles_type,typename std::shared_ptr<particles_type> >("/*{{particles.name}}*/")
	    .def(boost::python::vector_indexing_suite<particles_type>())
        /*{% for variable in variables %}*/
        .def("get_/*{{variable.type_name}}*/",&particles_type::get</*{{variable.type_name}}*/>)
        /*{% endfor %}*/
        .def("set",&particles_type::set)
	    ;


	/*
	 * Particle
	 */
	class_<particles_type::value_type, std::shared_ptr<particles_type::value_type> >("/*{{particles.particle_name}}*/",init<>())
        /*{% for variable in variables %}*/
		.add_property("/*{{variable.type_name}}*/", 
                make_function(&get</*{{variable.type_name}}*/,/*{{value_type_string}}*/ >,return_value_policy<copy_const_reference>()),
                &particles_type::set</*{{variable.type_name}}*/,/*{{value_type_string}}*/ >)
        /*{% endfor %}*/

}
