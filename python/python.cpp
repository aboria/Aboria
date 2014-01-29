/* 
 * python.cpp
 *
 * Copyright 2012 Martin Robinson
 *
 * This file is part of Tyche.
 *
 * Tyche is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Tyche is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with PDE_BD.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: Feb 9, 2013
 *      Author: mrobins
 */


#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <numpy/ndarrayobject.h>
#include "Tyche.h"
#include <numpy/arrayobject.h>

#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>


using namespace boost::python;

namespace Tyche {

struct ReactionEquation_from_python_list
{

	ReactionEquation_from_python_list() {
		boost::python::converter::registry::push_back(
				&convertible,
				&construct,
				boost::python::type_id<ReactionEquation>());
	}

	static void* convertible(PyObject* obj_ptr) {
		if (!PyList_Check(obj_ptr)) return 0;
		if (PyList_Size(obj_ptr) != 2) return 0;
		if (!PyList_Check(PyList_GetItem(obj_ptr,0)) || !PyList_Check(PyList_GetItem(obj_ptr,1))) return 0;
		return obj_ptr;
	}

	static void construct(
			PyObject* obj_ptr,
			boost::python::converter::rvalue_from_python_stage1_data* data) {

		// Extract the character data from the python string
		PyObject* reactants = PyList_GetItem(obj_ptr,0);
		PyObject* products = PyList_GetItem(obj_ptr,1);
		const int num_reactants = PyList_Size(reactants);
		const int num_products = PyList_Size(products);


		ReactionSide lhs;
		for (int i = 0; i < num_reactants; ++i) {
			Species* s = extract<Species*>(PyList_GetItem(reactants,i));
			lhs = lhs + *s;
		}
		ReactionSide rhs;
		for (int i = 0; i < num_products; ++i) {
			Species* s = extract<Species*>(PyList_GetItem(products,i));
			rhs = rhs + *s;
		}
		// Grab pointer to memory into which to construct the new QString
		void* storage = (
				(boost::python::converter::rvalue_from_python_storage<ReactionEquation>*)
				data)->storage.bytes;

		// in-place construct the new QString using the character data
		// extraced from the python object
		new (storage) ReactionEquation(lhs,rhs);

		// Stash the memory chunk pointer for later use by boost.python
		data->convertible = storage;
	}

};

template<typename T>
struct Vect3_from_python_list
{

	Vect3_from_python_list() {
		boost::python::converter::registry::push_back(
				&convertible,
				&construct,
				boost::python::type_id<Eigen::Matrix<T,3,1> >());
	}

	static void* convertible(PyObject* obj_ptr) {
		if (!PyList_Check(obj_ptr)) return 0;
		if (PyList_Size(obj_ptr) != 3) return 0;
		return obj_ptr;
	}

	static void construct(
			PyObject* obj_ptr,
			boost::python::converter::rvalue_from_python_stage1_data* data) {
		const int size = PyList_Size(obj_ptr);

		// Extract the character data from the python string
		const double x = extract<T>(PyList_GetItem(obj_ptr,0));
		const double y = extract<T>(PyList_GetItem(obj_ptr,1));
		const double z = extract<T>(PyList_GetItem(obj_ptr,2));


		// Grab pointer to memory into which to construct the new QString
		void* storage = (
				(boost::python::converter::rvalue_from_python_storage<Eigen::Matrix<T,3,1> >*)
				data)->storage.bytes;

		// in-place construct the new QString using the character data
		// extraced from the python object
		new (storage) Eigen::Matrix<T,3,1>(x,y,z);

		// Stash the memory chunk pointer for later use by boost.python
		data->convertible = storage;
	}

};

void python_init(boost::python::list& py_argv) {
	using boost::python::len;

	int argc = len(py_argv);
	char **argv = new char*[len(py_argv)];
	for (int i = 0; i < len(py_argv); ++i) {
		std::string this_argv = boost::python::extract<std::string>(py_argv[i]);
		argv[i] = new char[this_argv.size()];
		std::strcpy(argv[i],this_argv.c_str());
	}
	init(argc,argv);
	for (int i = 0; i < len(py_argv); ++i) {
		delete argv[i];
	}
	delete argv;

	// register the Vect3-to-python converter


}


BOOST_PYTHON_FUNCTION_OVERLOADS(new_uni_reaction_overloads, UniMolecularReaction::New, 2, 3);


std::auto_ptr<Operator> new_bi_reaction(const double rate, const ReactionEquation& eq,
					const double binding,
					const double unbinding,
					const double dt,
					const Vect3d& min, const Vect3d& max, const Vect3b& periodic,
					const bool reversible=false) {

	return BiMolecularReaction<BucketSort>::New(rate,eq,binding,unbinding,dt,min,max,periodic,reversible);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(new_bi_reaction_overloads, new_bi_reaction, 8, 9);

std::auto_ptr<Operator> new_bi_reaction2(const double rate, const ReactionEquation& eq,
					const double dt,
					const Vect3d& min, const Vect3d& max, const Vect3b& periodic,
					const bool reversible=false) {

	return BiMolecularReaction<BucketSort>::New(rate,eq,dt,min,max,periodic,reversible);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(new_bi_reaction_overloads2, new_bi_reaction2, 6, 7);


boost::python::object Species_get_compartments(Species& self) {
	Vect3i grid_size = self.grid->get_cells_along_axes();
	npy_intp size[3] = {grid_size[0],grid_size[1],grid_size[2]};
	PyObject *out = PyArray_SimpleNew(3, size, NPY_INT);
	for (int i = 0; i < grid_size[0]; ++i) {
		for (int j = 0; j < grid_size[1]; ++j) {
			for (int k = 0; k < grid_size[2]; ++k) {
				*((int *)PyArray_GETPTR3(out, i, j, k)) = self.copy_numbers[self.grid->vect_to_index(i,j,k)];
			}
		}
	}
	boost::python::handle<> handle(out);
	boost::python::numeric::array arr(handle);
	return arr;

}

boost::python::tuple Species_get_particles(Species& self) {
	const int n = self.mols.size();

	npy_intp size = {n};
	PyObject *out_x = PyArray_SimpleNew(1, &size, NPY_DOUBLE);
	PyObject *out_y = PyArray_SimpleNew(1, &size, NPY_DOUBLE);
	PyObject *out_z = PyArray_SimpleNew(1, &size, NPY_DOUBLE);

	for (int i = 0; i < n; ++i) {
		*((double *)PyArray_GETPTR1(out_x,i)) = self.mols.get_position(i)[0];
		*((double *)PyArray_GETPTR1(out_y,i)) = self.mols.get_position(i)[1];
		*((double *)PyArray_GETPTR1(out_z,i)) = self.mols.get_position(i)[2];
	}
	return boost::python::make_tuple(boost::python::handle<>(out_x),boost::python::handle<>(out_y),boost::python::handle<>(out_z));
}

boost::python::object Species_get_concentration1(Species& self, const Vect3d& min, const Vect3d& max, const Vect3i& n) {
	const Vect3d spacing = (max-min).cwiseQuotient(n.cast<double>());
	std::vector<double> concentration;
	StructuredGrid calc_grid(min,max,spacing);
	self.get_concentration(calc_grid,concentration);
	Vect3i grid_size = calc_grid.get_cells_along_axes();
	npy_intp size[3] = {grid_size[0],grid_size[1],grid_size[2]};

	PyObject *out = PyArray_SimpleNew(3, size, NPY_DOUBLE);
	for (int i = 0; i < grid_size[0]; ++i) {
		for (int j = 0; j < grid_size[1]; ++j) {
			for (int k = 0; k < grid_size[2]; ++k) {
				const int index = calc_grid.vect_to_index(i,j,k);
				*((double *)PyArray_GETPTR3(out, i, j, k)) = concentration[index];
			}
		}
	}
	boost::python::handle<> handle(out);
	boost::python::numeric::array arr(handle);
	return arr;
//	std::cout <<"c after set PyObject"<<std::endl;
//	object obj(handle<>(out));
//	return extract<numeric::array>(obj);
}

boost::python::object BindingReaction_get_state_sequence(BindingReaction& self) {
	boost::python::list retlist = boost::python::list();
	std::list<std::pair<int, double > > slist = self.get_state_sequence(true);
	std::list<std::pair<int, double > >::const_iterator iter;

	for (iter = slist.begin(); iter != slist.end(); ++iter) {
		int state = (*iter).first;
		double time = (*iter).second;
		retlist.append(boost::python::make_tuple(state, time));
	}

	return retlist;
}

struct BR_Python_Callback {
  boost::python::object py_callable;
  boost::python::tuple args;
  void operator()(int state) {
    py_callable(state, args);
  }
};

void BindingReaction_set_state_changed_cb(BindingReaction& self, boost::python::object& callable, boost::python::tuple& args)
{
  BR_Python_Callback cb {callable, args};
  self.set_state_changed_cb(boost::function<void(int)>(cb));
}

std::auto_ptr<Operator> group(const boost::python::list& ops) {
	OperatorList* result = new OperatorList();
	const int n = len(ops);
	for (int i = 0; i < n; ++i) {
		Operator *s = extract<Operator*>(ops[i]);
		result->push_back(s);
	}
	return std::auto_ptr<Operator>(result);
}


template<class T>
struct vtkSmartPointer_to_python {
	static PyObject *convert(const vtkSmartPointer<T> &p) {
		std::ostringstream oss;
		oss << (void*) p.GetPointer();
		std::string address_str = oss.str();

		using namespace boost::python;
		object obj = import("vtk").attr("vtkObject")(address_str);
		return incref(obj.ptr());
	}
};



BOOST_PYTHON_MODULE(pyTyche) {
	import_array();
	numeric::array::set_module_and_type("numpy", "ndarray");
	def("init", python_init);
	def("random_seed", random_seed);

	/*
	 * vector
	 */
	class_<std::vector<double> >("Vectd")
	        .def(boost::python::vector_indexing_suite<std::vector<double> >())
	    ;

	to_python_converter<vtkSmartPointer<vtkUnstructuredGrid>,vtkSmartPointer_to_python<vtkUnstructuredGrid> >();
	Vect3_from_python_list<double>();
	Vect3_from_python_list<int>();
	Vect3_from_python_list<bool>();
	ReactionEquation_from_python_list();

	/*
	 * Species
	 */

	class_<Species,typename std::auto_ptr<Species> >("Species",boost::python::init<double>())
			.def("fill_uniform",&Species::fill_uniform<Box>)
			.def("fill_uniform",&Species::fill_uniform<Sphere>)
			.def("fill_uniform",&Species::fill_uniform<Shell>)

			.def("count_mols_in",&Species::count_mols_in<Box>)
			.def("count_mols_in",&Species::count_mols_in<Sphere>)
			.def("count_mols_in",&Species::count_mols_in<Shell>)

			.def("get_concentration",Species_get_concentration1)
			.def("get_vtk",&Species::get_vtk)
			.def("get_compartments",Species_get_compartments)
			.def("get_particles",Species_get_particles)
			.def(self_ns::str(self_ns::self))
			;
	def("new_species",Species::New);


   /*
    * Operator
    */
	class_<Operator,typename std::auto_ptr<Operator> >("Operator",boost::python::no_init)
			.def("integrate",&Operator::operator())
			.def("reset",&Operator::reset)
			.def("add_species",&Operator::add_species)
			.def("get_species_index",&Operator::get_species_index)
			.def("integrate_for_time",&Operator::integrate_for_time)
			.def("get_active", &Operator::get_active)
	        .def("set_active", &Operator::set_active)
			.def(self_ns::str(self_ns::self))
			;

	def("group",group);



	/*
	 * Geometry
	 */
	def("new_xplane",xplane::New);
	def("new_yplane",yplane::New);
	def("new_zplane",zplane::New);
	def("new_xrect",xrect::New);
	def("new_yrect",xrect::New);
	def("new_zrect",xrect::New);


	class_<xplane,typename std::auto_ptr<xplane> >("Xplane",boost::python::no_init);
	class_<yplane,typename std::auto_ptr<yplane> >("Yplane",boost::python::no_init);
	class_<zplane,typename std::auto_ptr<zplane> >("Zplane",boost::python::no_init);

	class_<xrect,typename std::auto_ptr<xrect> >("Xrect",boost::python::no_init);
	class_<yrect,typename std::auto_ptr<yrect> >("Yrect",boost::python::no_init);
	class_<zrect,typename std::auto_ptr<zrect> >("Zrect",boost::python::no_init);

	def("new_box", Box::New);
	def("new_multiple_boxes", MultipleBoxes::New);

	class_<Box,typename std::auto_ptr<Box> >("Box",boost::python::no_init)
						;
	class_<MultipleBoxes,typename std::auto_ptr<MultipleBoxes> >("MultipleBoxes",boost::python::no_init)
				.def("add_hole",&MultipleBoxes::add_box)
				;

	def("new_sphere",Sphere::New);
	def("new_shell",Shell::New);

	class_<Sphere,typename std::auto_ptr<Sphere> >("Sphere",boost::python::no_init)
				.add_property("centre", boost::python::make_function( &Sphere::get_centre, return_internal_reference<>() ), &Sphere::set_centre)
				.add_property("radius", &Sphere::get_radius, &Sphere::set_radius)
				;
	class_<Shell,typename std::auto_ptr<Shell> >("Shell",boost::python::no_init)
				;

	/*
	 * Boundaries
	 */
	def("new_coupling_boundary",CouplingBoundary<xplane>::New);
	def("new_coupling_boundary",CouplingBoundary<yplane>::New);
	def("new_coupling_boundary",CouplingBoundary<zplane>::New);
	def("new_coupling_boundary",CouplingBoundary<xrect>::New);
	def("new_coupling_boundary",CouplingBoundary<yrect>::New);
	def("new_coupling_boundary",CouplingBoundary<Box>::New);


    def("new_reflective_boundary",ReflectiveBoundary<xplane>::New);
    def("new_reflective_boundary",ReflectiveBoundary<yplane>::New);
    def("new_reflective_boundary",ReflectiveBoundary<zplane>::New);
    def("new_reflective_boundary",ReflectiveBoundary<xrect>::New);
    def("new_reflective_boundary",ReflectiveBoundary<yrect>::New);
    def("new_reflective_boundary",ReflectiveBoundary<zrect>::New);


    def("new_jump_boundary",JumpBoundary<xplane>::New);
    def("new_jump_boundary",JumpBoundary<yplane>::New);
    def("new_jump_boundary",JumpBoundary<zplane>::New);
    def("new_jump_boundary",JumpBoundary<xrect>::New);
    def("new_jump_boundary",JumpBoundary<yrect>::New);
    def("new_jump_boundary",JumpBoundary<zrect>::New);


    /*
     * Diffusion
     */
    def("new_diffusion",Diffusion::New);

    /*
     * Reactions
     */
    def("new_zero_reaction",ZeroOrderMolecularReaction::New);


    def("new_uni_reaction",UniMolecularReaction::New, new_uni_reaction_overloads());

    def("new_bi_reaction",new_bi_reaction, new_bi_reaction_overloads());
    def("new_bi_reaction",new_bi_reaction2, new_bi_reaction_overloads2());
    def("new_tri_reaction",TriMolecularReaction::New);
    def("new_binding_reaction", BindingReaction::New);

	class_<BindingReaction, bases<Operator>, std::auto_ptr<BindingReaction> >("BindingReaction", boost::python::no_init)
		.def("get_site_state", &BindingReaction::get_site_state)
		.def("get_state_sequence", &BindingReaction_get_state_sequence)
	        .def("set_state_changed_cb", &BindingReaction_set_state_changed_cb);

    /*
     * Compartments
     */
    def("new_compartments",NextSubvolumeMethod::New);

    /*
     * NextSubvolume
     */
    class_<NextSubvolumeMethod, bases<Operator>, std::auto_ptr<NextSubvolumeMethod> >("NextSubvolumeMethod",boost::python::no_init)
    	.def("set_interface",&NextSubvolumeMethod::set_interface<xplane>)
    	.def("set_interface",&NextSubvolumeMethod::set_interface<yplane>)
    	.def("set_interface",&NextSubvolumeMethod::set_interface<zplane>)
    	.def("set_interface",&NextSubvolumeMethod::set_interface<xrect>)
    	.def("set_interface",&NextSubvolumeMethod::set_interface<yrect>)
    	.def("set_interface",&NextSubvolumeMethod::set_interface<zrect>)
    	.def("set_interface",&NextSubvolumeMethod::set_interface<Box>)
    	.def("set_interface",&NextSubvolumeMethod::set_interface<MultipleBoxes>)
    	.def("unset_interface",&NextSubvolumeMethod::unset_interface<xplane>)
    	.def("unset_interface",&NextSubvolumeMethod::unset_interface<yplane>)
    	.def("unset_interface",&NextSubvolumeMethod::unset_interface<zplane>)
    	.def("unset_interface",&NextSubvolumeMethod::unset_interface<xrect>)
    	.def("unset_interface",&NextSubvolumeMethod::unset_interface<yrect>)
    	.def("unset_interface",&NextSubvolumeMethod::unset_interface<zrect>)
    	.def("unset_interface",&NextSubvolumeMethod::unset_interface<Box>)
    	.def("unset_interface",&NextSubvolumeMethod::unset_interface<MultipleBoxes>)
    	.def("add_diffusion",&NextSubvolumeMethod::add_diffusion)
    	.def("add_diffusion_between",&NextSubvolumeMethod::add_diffusion_between<xplane>)
    	.def("add_diffusion_between",&NextSubvolumeMethod::add_diffusion_between<yplane>)
    	.def("add_diffusion_between",&NextSubvolumeMethod::add_diffusion_between<zplane>)
    	.def("add_diffusion_between",&NextSubvolumeMethod::add_diffusion_between<xrect>)
    	.def("add_diffusion_between",&NextSubvolumeMethod::add_diffusion_between<yrect>)
    	.def("add_diffusion_between",&NextSubvolumeMethod::add_diffusion_between<zrect>)
    	.def("add_reaction",&NextSubvolumeMethod::add_reaction)
    	.def("add_reaction_on",&NextSubvolumeMethod::add_reaction_on<xplane>)
    	.def("add_reaction_on",&NextSubvolumeMethod::add_reaction_on<yplane>)
    	.def("add_reaction_on",&NextSubvolumeMethod::add_reaction_on<zplane>)
    	.def("add_reaction_on",&NextSubvolumeMethod::add_reaction_on<xrect>)
    	.def("add_reaction_on",&NextSubvolumeMethod::add_reaction_on<yrect>)
    	.def("add_reaction_on",&NextSubvolumeMethod::add_reaction_on<zrect>)
    	.def("fill_uniform",&NextSubvolumeMethod::fill_uniform)
    	;

}

}
