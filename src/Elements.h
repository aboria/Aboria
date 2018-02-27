/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Aboria.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef ELEMENTS_H_
#define ELEMENTS_H_

#include <random>
#include <string>
#include <vector>
//#include <boost/array.hpp>
//#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/adaptors.hpp>

#include "CellList.h"
#include "Get.h"
#include "Traits.h"
#include "Variable.h"
#include "Vector.h"
#include "Zip.h"
//#include "OctTree.h"
#include "CudaInclude.h"

#ifdef HAVE_VTK
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>
#endif

#include "Particles.h"

namespace Aboria {

template <typename ParticlesType, typename VariableType,
          typename VAR = std::tuple<>, unsigned int SelfD = 2>
class Elements {
public:
  ///
  /// This type
  typedef Elements<ParticlesType, VariableType, VAR, SelfD> elements_type;

  typedef ParticlesType particles_type;
  typedef VariableType variable_type;

  ///
  /// The traits type used to build up the Elements container.
  /// Contains Level 0 vector class and dimension information
  typedef TraitsCommon<VAR, ParticlesType::dimension, SelfD,
                       typename ParticlesType::TRAITS_USER>
      traits_type;

  ///
  /// a tuple type containing value_types for each Variable
  typedef typename traits_type::value_type value_type;

  ///
  /// a tuple type containing references to value_types for each Variable
  typedef typename traits_type::reference reference;

  ///
  /// a tuple type containing const_references to value_types for each Variable
  typedef typename traits_type::const_reference const_reference;

  ///
  /// a tuple type containing raw references to value_types for each Variable
  typedef typename traits_type::raw_pointer raw_pointer;

  ///
  /// a tuple type containing raw references to value_types for each Variable
  typedef typename traits_type::raw_reference raw_reference;

  ///
  /// a tuple type containing raw references to value_types for each Variable
  typedef typename traits_type::raw_const_reference raw_const_reference;

  ///
  /// type used to hold data (a tuple of vectors)
  typedef typename traits_type::data_type data_type;

  ///
  /// type used to hold and return size information
  typedef typename traits_type::size_type size_type;

  ///
  /// type used to hold and return the difference of iterator, or distance
  /// between
  typedef typename traits_type::size_type difference_type;

  ///
  /// non-const iterator type
  /// \see zip_iterator
  typedef typename traits_type::iterator iterator;

  ///
  /// const iterator type
  /// \see zip_iterator
  typedef typename traits_type::const_iterator const_iterator;

  /// a boost mpl vector type containing a vector of Variable
  /// attached to the particles (includes position, id and
  /// alive flag as well as all user-supplied variables)
  typedef typename traits_type::mpl_type_vector mpl_type_vector;

  template <typename T>
  using elem_by_type = detail::get_elem_by_type<T, mpl_type_vector>;
  template <typename T>
  using return_type =
      typename detail::getter_helper<typename data_type::tuple_type>::
          template return_type<elem_by_type<T>::index>;

  ///
  /// the number of spatial dimensions
  static const unsigned int dimension = particles_type::dimension;
  static const unsigned int element_dimension = SelfD;

  ///
  /// a type to store a vector of doubles with given dimension
  typedef Vector<double, dimension> double_d;

  ///
  /// a type to store a vector of doubles with given dimension
  typedef Vector<double, dimension> int_d;

  ///
  /// a type to store a vector of bool with given dimension
  typedef Vector<bool, dimension> bool_d;

  ///
  /// the tag type for the default particles variable
  typedef typename traits_type::position particles;

  /// Contructs an empty container with no searching or id tracking enabled
  Elements(const particles_type &particles)
      : m_particles(&particles), m_next_id(0), m_seed(time(NULL)) {}

  /// push the particle \p val to the back of the container (if it is within
  /// the searchable domain)
  ///
  /// \param val the value_type to add
  void push_back(const value_type &val) {
    // add val to container
    traits_type::push_back(m_data, val);

    // overwrite id, alive and random generator
    reference i = *(end() - 1);
    Aboria::get<id>(i) = this->m_next_id++;
    Aboria::get<generator>(i) =
        generator_type((m_seed + uint32_t(Aboria::get<id>(i))));
    Aboria::get<alive>(i) = true;

    push_connections(end() - 1, end());
  }

  /// set the base seed of the container. Note that the random number generator
  /// for each particle is set to \p value plus the particle's id
  void set_seed(const uint32_t value) {
    m_seed = value;
    detail::for_each(begin(), end(),
                     detail::set_seed_lambda<raw_reference>(m_seed));
  }

  /// returns a reference to the element at index \p idx
  reference operator[](std::size_t idx) {
    return traits_type::index(m_data, idx);
  }

  /// returns a const_reference to the element at index \p idx
  const_reference operator[](std::size_t idx) const {
    return traits_type::index(m_data, idx);
  }

  /// returns an iterator to the beginning of the container
  iterator begin() { return traits_type::begin(m_data); }

  /// returns an iterator to the end of the container
  iterator end() { return traits_type::end(m_data); }

  /// returns an iterator to the beginning of the container
  const_iterator begin() const { return traits_type::cbegin(m_data); }

  /// returns an iterator to the end of the container
  const_iterator end() const { return traits_type::cend(m_data); }

  /// returns a const_iterator to the beginning of the container
  const_iterator cbegin() const { return traits_type::cbegin(m_data); }

  /// returns an iterator to the end of the container
  const_iterator cend() const { return traits_type::cend(m_data); }

  /// sets container to empty and deletes all particles
  void clear() {
    for (size_t i = 0; i < size(); ++i) {
      for (size_t d = 0; d < element_dimension; ++d) {
        size_t particle_id = get<particles>(m_data)[i][d];
        // TODO: this will not work with thrust
        typename ParticlesType::raw_pointer p =
            m_particles->get_query().find(particle_id);
        CHECK(p != iterator_to_raw_pointer(m_particles->end()),
              "particle " << particle_id << " does not exist");
        get<VariableType>(p).clear(i);
      }
    }

    return traits_type::clear(m_data);
  }

  /// return the total number of elements in the container
  size_t size() const { return Aboria::get<particles>(m_data).size(); }

  /// Update the particle connections for particles between
  /// \p update_begin and \p update_end (not including \p update_end).
  /// It is assumed that this range is all new elements
  ///
  void push_connections(iterator update_begin, iterator update_end) {
    for (size_t i = update_begin - begin(); i < update_end - end(); ++i) {
      for (size_t d = 0; d < element_dimension; ++d) {
        size_t particle_id = get<particles>(m_data)[i][d];
        // TODO: this will not work with thrust
        typename ParticlesType::raw_pointer p =
            m_particles->get_query().find(particle_id);
        CHECK(p != iterator_to_raw_pointer(m_particles->end()),
              "particle " << particle_id << " does not exist");
        get<VariableType>(p).push_back(i);
      }
    }
  }

  // Need to be mark as device to enable get functions being device/host
  CUDA_HOST_DEVICE
  const typename data_type::tuple_type &get_tuple() const {
    return m_data.get_tuple();
  }

  // Need to be mark as device to enable get functions being device/host
  CUDA_HOST_DEVICE
  typename data_type::tuple_type &get_tuple() { return m_data.get_tuple(); }

  /// stream output for element data
  friend std::ostream &operator<<(std::ostream &stream,
                                  const Elements &elements) {
    traits_type::header_to_stream(stream);
    stream << '\n';
    for (const_iterator i = elements.cbegin(); i != elements.cend(); ++i) {
      traits_type::to_stream(i, stream);
      stream << '\n';
    }
    return stream;
  }

  /// stream input for element data
  friend std::istream &operator>>(std::istream &stream, Elements &elements) {
    stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    for (iterator i = elements.begin(); i != elements.end(); ++i) {
      traits_type::from_stream(i);
    }
    return stream;
  }

  /// Boost serialization support
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    traits_type::serialize(m_data, ar, version);
  }

#ifdef HAVE_VTK

  /// get a vtk unstructured grid version of the particle container
  /// This grid is cached internally. The first time this function is
  /// called the grid is created and updated with the particle data.
  /// If the particle data is updated subsequently you need to set \p
  /// refresh=true for the grid to be updated
  vtkSmartPointer<vtkUnstructuredGrid> get_grid(bool refresh = false) {
    if (!m_cache_grid) {
      m_cache_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
      copy_to_vtk_grid(m_cache_grid);
    } else if (refresh == true) {
      copy_to_vtk_grid(m_cache_grid);
    }
    return m_cache_grid;
  }

  ///  copy the particle data to a VTK unstructured grid
  void copy_to_vtk_grid(vtkUnstructuredGrid *grid) {
    CHECK(element_dimension < 2,
          "cannot write elements with dimension greater than 1");
    size_t cell_type = VTK_VERTEX;
    if (element_dimension == 1) {
      cell_type = VTK_LINE;
    }

    vtkSmartPointer<vtkPoints> points = grid->GetPoints();
    if (!points) {
      points = vtkSmartPointer<vtkPoints>::New();
      grid->SetPoints(points);
    }
    vtkSmartPointer<vtkCellArray> cells = grid->GetCells();
    if (!cells) {
      cells = vtkSmartPointer<vtkCellArray>::New();
      grid->SetCells(cell_type, cells);
    }
    vtkSmartPointer<vtkUnsignedCharArray> cell_types =
        grid->GetCellTypesArray();

    const vtkIdType np = this->get_particles().size();
    const vtkIdType nc = this->size();

    constexpr size_t dn = mpl::size<mpl_type_vector>::type::value;
    vtkSmartPointer<vtkFloatArray> datas[dn];
    mpl::for_each<mpl::range_c<int, 1, dn>>(
        detail::setup_datas_for_writing<reference>(nc, datas, grid));
    for (int i = 0; i < dn; ++i) {
      grid->GetCellData()->AddArray(datas[i]);
    }
    points->SetNumberOfPoints(np);
    cells->Reset();
    cell_types->Reset();
    int j = 0;

    double write_point[3];
    const unsigned int max_d = std::min(3u, dimension);
    for (auto i : this->get_particles()) {
      const int index = j++;
      // std::cout <<"copying point at "<<i.get_position()<<" with id =
      // "<<i.get_id()<<std::endl;
      const double_d &r = get<typename particles_type::position>(i);
      for (int d = 0; d < max_d; ++d) {
        write_point[d] = r[d];
      }
      points->SetPoint(index, write_point);
    }
    for (auto i : *this) {
      cells->InsertNextCell(element_dimension);
      for (int j = 0; j < element_dimension; ++j) {
        auto p = this->get_particles().get_query().find(get<particles>(i)[j]);
        ASSERT(p != iterator_to_raw_pointer(this->get_particles().end()),
               "cannot find particle");
        const size_t index =
            get<typename particles_type::position>(p) -
            get<typename particles_type::position>(this->get_particles())
                .data();

        cells->InsertCellPoint(index);
      }
      cell_types->InsertNextTuple1(cell_type);

      mpl::for_each<mpl::range_c<int, 1, dn>>(
          detail::write_from_tuple<reference>(
              i.get_tuple(), index, datas,
              m_seed + uint32_t(Aboria::get<id>(i))));
    }
  }

  ///  update the particle data according to a supplied VTK unstructured grid
  ///  it is assumed that \p grid has been generated from copy_to_vtk_grid() or
  ///  get_grid()
  ///  \see get_grid()
  ///  \see copy_to_vtk_grid()
  void copy_from_vtk_grid(vtkUnstructuredGrid *grid) {
    vtkSmartPointer<vtkPoints> points = grid->GetPoints();
    CHECK(points, "No points in vtkUnstructuredGrid");
    vtkSmartPointer<vtkCellArray> cells = grid->GetCells();
    CHECK(points, "No cells in vtkUnstructuredGrid");
    vtkSmartPointer<vtkUnsignedCharArray> cell_types =
        grid->GetCellTypesArray();
    constexpr size_t dn = mpl::size<mpl_type_vector>::type::value;

    vtkSmartPointer<vtkFloatArray> datas[dn];

    const vtkIdType nc = cells->GetNumberOfCells();

    mpl::for_each<mpl::range_c<int, 1, dn>>(
        detail::setup_datas_for_reading<reference>(nc, datas, grid, false));

    this->clear();

    for (int j = 0; j < nc; ++j) {
      value_type element;
      mpl::for_each<mpl::range_c<int, 1, dn>>(
          detail::read_into_tuple<reference>(element.get_tuple(), j, datas));
      this->push_back(element);
    }

    push_connections(begin(), end());
  }
#endif

private:
  /// a pointer to the particles data structure that containes the element
  /// particles
  const particles_type *m_particles;

  /// Contains the element data, implemented as a std::tuple of Level 0 vectors
  data_type m_data;

  /// The next available id number
  int m_next_id;

  /// The base random seed for the container
  uint32_t m_seed;

#ifdef HAVE_VTK
  /// An vtkUnstructuredGrid to store particle data in (if neccessary)
  vtkSmartPointer<vtkUnstructuredGrid> m_cache_grid;
#endif
};

namespace detail {

template <unsigned int D, typename T> struct is_elements : std::false_type {};

template <unsigned int D, typename Particles, typename Variable, typename Var>
struct is_elements<D, Elements<Particles, Variable, Var, D>> : std::true_type {
};

} // namespace detail

} // namespace Aboria

#endif /* SPECIES_H_ */
