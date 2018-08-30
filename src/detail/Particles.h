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

#ifndef PARTICLES_DETAIL_H_
#define PARTICLES_DETAIL_H_

#include "CudaInclude.h"

namespace Aboria {

namespace detail {

template <typename T> struct is_particles : std::false_type {};

template <typename Variables, unsigned int DomainD,
          template <typename, typename> class Vector,
          template <typename> class SearchMethod, typename Traits>
struct is_particles<Particles<Variables, DomainD, Vector, SearchMethod, Traits>>
    : std::true_type {};

template <typename Reference> struct resize_lambda {
  uint32_t seed;
  int next_id;
  const size_t *start_id_pointer;

  resize_lambda(const uint32_t &seed, const int &next_id,
                const size_t *start_id_pointer)
      : seed(seed), next_id(next_id), start_id_pointer(start_id_pointer) {}

  CUDA_HOST_DEVICE
  void operator()(Reference i) const {
    Aboria::get<alive>(i) = true;

    const size_t index = &Aboria::get<id>(i) - start_id_pointer;
    Aboria::get<id>(i) = index + next_id;

    generator_type &gen = Aboria::get<generator>(i);
    gen.seed(seed + uint32_t(Aboria::get<id>(i)));
  }
};

template <typename Reference> struct set_seed_lambda {
  uint32_t seed;

  set_seed_lambda(const uint32_t &seed) : seed(seed) {}

  CUDA_HOST_DEVICE
  void operator()(Reference i) const {
    generator_type &gen = Aboria::get<generator>(i);
    gen.seed(seed + uint32_t(Aboria::get<id>(i)));
  }
};

template <typename ConstReference> struct is_alive {
  CUDA_HOST_DEVICE
  bool operator()(ConstReference i) const { return Aboria::get<alive>(i); }
};

#ifdef HAVE_VTK

template <typename reference> struct write_from_tuple {};

template <typename MplVector, typename... Types>
struct write_from_tuple<getter_type<std::tuple<Types...>, MplVector>> {
  typedef std::tuple<Types...> tuple_type;
  typedef MplVector mpl_type_vector;

  template <typename U>
  using non_ref_tuple_element = typename std::remove_reference<
      typename std::tuple_element<U::value, tuple_type>::type>::type;

  write_from_tuple(tuple_type write_from, int index,
                   vtkSmartPointer<vtkFloatArray> *datas, const uint32_t &seed)
      : write_from(write_from), index(index), datas(datas), seed(seed) {}

  template <typename U>
  typename boost::enable_if<
      boost::is_arithmetic<non_ref_tuple_element<U>>>::type
  operator()(U i) {
    if (datas[i] != nullptr) {
      datas[i]->SetValue(index, std::get<U::value>(write_from));
    }
  }

#ifdef HAVE_EIGEN
  template <typename U>
  typename boost::enable_if<is_eigen_vector<non_ref_tuple_element<U>>>::type
  operator()(U i) {
    if (datas[i] != nullptr) {
      Eigen::Matrix<float, non_ref_tuple_element<U>::RowsAtCompileTime, 1> tmp =
          std::get<U::value>(write_from).template cast<float>();
      datas[i]->SetTuple(index, tmp.data());
    }
  }
#endif

  template <typename U>
  typename boost::enable_if<is_vector<non_ref_tuple_element<U>>>::type
  operator()(U i) {
    if (datas[i] != nullptr) {
      Vector<float, non_ref_tuple_element<U>::size> tmp =
          std::get<U::value>(write_from);
      datas[i]->SetTuple(index, tmp.data());
    }
  }

  template <typename U>
  typename boost::enable_if<
      boost::is_same<non_ref_tuple_element<U>, generator_type>>::type
  operator()(U i) {
    // TODO: not sure what to write here, default to original seed
    // seed + uint32_t(Aboria::get<id>(i)
    if (datas[i] != nullptr) {
      datas[i]->SetValue(index, seed);
    }
  }

  template <typename U>
  typename boost::enable_if<
      boost::is_same<non_ref_tuple_element<U>, uint32_t>>::type
  operator()(U i) {
    if (datas[i] != nullptr) {
      datas[i]->SetValue(index, std::get<U::value>(write_from));
    }
  }

  template <typename U>
  typename std::enable_if<
      !std::is_same<non_ref_tuple_element<U>, uint32_t>::value &&
      !std::is_same<non_ref_tuple_element<U>, generator_type>::value &&
      !is_vector<non_ref_tuple_element<U>>::value &&
      !is_eigen_vector<non_ref_tuple_element<U>>::value &&
      !boost::is_arithmetic<non_ref_tuple_element<U>>::value>::type
  operator()(U i) {
    // do nothing
  }

  tuple_type write_from;
  int index;
  vtkSmartPointer<vtkFloatArray> *datas;
  uint32_t seed;
};

#ifdef HAVE_THRUST
template <typename MplVector, typename TT1, typename TT2, typename TT3,
          typename TT4, typename TT5, typename TT6, typename TT7, typename TT8,
          typename TT9>
struct write_from_tuple<getter_type<
    thrust::tuple<TT1, TT2, TT3, TT4, TT5, TT6, TT7, TT8, TT9>, MplVector>> {
  typedef thrust::tuple<TT1, TT2, TT3, TT4, TT5, TT6, TT7, TT8, TT9> tuple_type;
  typedef MplVector mpl_type_vector;

  template <typename U>
  using non_ref_tuple_element =
      typename std::remove_reference<typename thrust::detail::raw_reference<
          typename thrust::tuple_element<U::value, tuple_type>::type>::type>::
          type;

  write_from_tuple(tuple_type write_from, int index,
                   vtkSmartPointer<vtkFloatArray> *datas, const uint32_t &seed)
      : write_from(write_from), index(index), datas(datas), seed(seed) {}

  template <typename U>
  typename boost::enable_if<
      boost::is_arithmetic<non_ref_tuple_element<U>>>::type
  operator()(U i) {
    datas[i]->SetValue(index, thrust::get<U::value>(write_from));
  }

  template <typename U>
  typename boost::enable_if<is_vector<non_ref_tuple_element<U>>>::type
  operator()(U i) {
    Vector<float, non_ref_tuple_element<U>::size> tmp =
        static_cast<Vector<double, non_ref_tuple_element<U>::size>>(
            thrust::get<U::value>(write_from));
    datas[i]->SetTuple(index, tmp.data());
  }

#ifdef HAVE_EIGEN
  template <typename U>
  typename boost::enable_if<is_eigen_vector<non_ref_tuple_element<U>>>::type
  operator()(U i) {
    Eigen::Matrix<float, non_ref_tuple_element<U>::RowsAtCompileTime, 1> tmp =
        static_cast<Vector<double, non_ref_tuple_element<U>::size>>(
            thrust::get<U::value>(write_from))
            .template cast<float>();
    datas[i]->SetTuple(index, tmp.data());
  }
#endif

  template <typename U>
  typename boost::enable_if<
      boost::is_same<non_ref_tuple_element<U>, generator_type>>::type
  operator()(U i) {
    // TODO: not sure what to write here, default to original seed
    // seed + uint32_t(Aboria::get<id>(i)
    datas[i]->SetValue(index, seed);
  }

  template <typename U>
  typename boost::enable_if<
      boost::is_same<non_ref_tuple_element<U>, uint32_t>>::type
  operator()(U i) {
    datas[i]->SetValue(index, thrust::get<U::value>(write_from));
  }

  tuple_type write_from;
  int index;
  vtkSmartPointer<vtkFloatArray> *datas;
  uint32_t seed;
};
#endif

template <typename reference> struct read_into_tuple {
  typedef typename reference::tuple_type tuple_type;
  typedef typename reference::mpl_type_vector mpl_type_vector;
  template <typename U>
  using non_ref_tuple_element = typename std::remove_reference<
      typename std::tuple_element<U::value, tuple_type>::type>::type;

  read_into_tuple(tuple_type &read_into, int index,
                  vtkSmartPointer<vtkFloatArray> *datas)
      : read_into(read_into), index(index), datas(datas) {}

  template <typename U>
  typename boost::enable_if<
      boost::is_arithmetic<non_ref_tuple_element<U>>>::type
  operator()(U i) {
    std::get<U::value>(read_into) = datas[i]->GetValue(index);
  }

  template <typename U>
  typename boost::enable_if<is_vector<non_ref_tuple_element<U>>>::type
  operator()(U i) {
    datas[i]->GetTuple(index, std::get<U::value>(read_into).data());
  }

#ifdef HAVE_EIGEN
  template <typename U>
  typename boost::enable_if<is_eigen_vector<non_ref_tuple_element<U>>>::type
  operator()(U i) {
    datas[i]->GetTuple(index, std::get<U::value>(read_into).data());
  }
#endif

  template <typename U>
  typename boost::enable_if<
      boost::is_same<non_ref_tuple_element<U>, generator_type>>::type
  operator()(U i) {
    // TODO: not sure what to read here, abandon it for now
    // datas[i]->SetValue(index,seed);
  }

  tuple_type read_into;
  int index;
  vtkSmartPointer<vtkFloatArray> *datas;
};

template <typename reference> struct setup_datas_for_writing {};

template <typename MplVector, typename... Types>
struct setup_datas_for_writing<getter_type<std::tuple<Types...>, MplVector>> {
  typedef std::tuple<Types...> tuple_type;
  typedef MplVector mpl_type_vector;

  template <typename U>
  using non_ref_tuple_element = typename std::remove_reference<
      typename std::tuple_element<U::value, tuple_type>::type>::type;

  setup_datas_for_writing(size_t n, vtkSmartPointer<vtkFloatArray> *datas,
                          vtkUnstructuredGrid *grid)
      : n(n), datas(datas), grid(grid) {}

  template <typename U>
  typename std::enable_if<
      !is_vector<non_ref_tuple_element<U>>::value &&
      !is_eigen_vector<non_ref_tuple_element<U>>::value>::type
  operator()(U i) {
    typedef typename mpl::at<mpl_type_vector, U>::type variable_type;
    const char *name = variable_type().name;
    if (name[0] == '_') {
      return;
    }
    datas[i] =
        vtkFloatArray::SafeDownCast(grid->GetPointData()->GetArray(name));
    if (!datas[i]) {
      datas[i] = vtkSmartPointer<vtkFloatArray>::New();
      datas[i]->SetName(name);
    }
    datas[i]->SetNumberOfComponents(1);
    datas[i]->SetNumberOfTuples(n);
  }

  template <typename U>
  typename boost::enable_if<is_vector<non_ref_tuple_element<U>>>::type
  operator()(U i) {
    typedef typename mpl::at<mpl_type_vector, U>::type variable_type;
    const char *name = variable_type().name;
    if (name[0] == '_') {
      return;
    }
    datas[i] =
        vtkFloatArray::SafeDownCast(grid->GetPointData()->GetArray(name));
    if (!datas[i]) {
      datas[i] = vtkSmartPointer<vtkFloatArray>::New();
      datas[i]->SetName(name);
    }
    datas[i]->SetNumberOfComponents(non_ref_tuple_element<U>::size);
    datas[i]->SetNumberOfTuples(n);
  }

#ifdef HAVE_EIGEN
  template <typename U>
  typename boost::enable_if<is_eigen_vector<non_ref_tuple_element<U>>>::type
  operator()(U i) {
    typedef typename mpl::at<mpl_type_vector, U>::type variable_type;
    const char *name = variable_type().name;
    if (name[0] == '_') {
      return;
    }
    datas[i] =
        vtkFloatArray::SafeDownCast(grid->GetPointData()->GetArray(name));
    if (!datas[i]) {
      datas[i] = vtkSmartPointer<vtkFloatArray>::New();
      datas[i]->SetName(name);
    }
    datas[i]->SetNumberOfComponents(
        non_ref_tuple_element<U>::RowsAtCompileTime);
    datas[i]->SetNumberOfTuples(n);
  }
#endif

  size_t n;
  vtkSmartPointer<vtkFloatArray> *datas;
  vtkUnstructuredGrid *grid;
  bool particles;
};

#ifdef HAVE_THRUST
template <typename MplVector, typename TT1, typename TT2, typename TT3,
          typename TT4, typename TT5, typename TT6, typename TT7, typename TT8,
          typename TT9>
struct setup_datas_for_writing<getter_type<
    thrust::tuple<TT1, TT2, TT3, TT4, TT5, TT6, TT7, TT8, TT9>, MplVector>> {
  typedef thrust::tuple<TT1, TT2, TT3, TT4, TT5, TT6, TT7, TT8, TT9> tuple_type;
  typedef MplVector mpl_type_vector;

  template <typename U>
  using non_ref_tuple_element =
      typename std::remove_reference<typename thrust::detail::raw_reference<
          typename thrust::tuple_element<U::value, tuple_type>::type>::type>::
          type;

  setup_datas_for_writing(size_t n, vtkSmartPointer<vtkFloatArray> *datas,
                          vtkUnstructuredGrid *grid)
      : n(n), datas(datas), grid(grid) {}

  template <typename U>
  typename std::enable_if<
      !is_vector<non_ref_tuple_element<U>>::value &&
      !is_eigen_vector<non_ref_tuple_element<U>>::value>::type
  operator()(U i) {
    typedef typename mpl::at<mpl_type_vector, U>::type variable_type;
    const char *name = variable_type().name;
    if (name[0] == '_') {
      return;
    }
    datas[i] =
        vtkFloatArray::SafeDownCast(grid->GetPointData()->GetArray(name));
    if (!datas[i]) {
      datas[i] = vtkSmartPointer<vtkFloatArray>::New();
      datas[i]->SetName(name);
    }
    datas[i]->SetNumberOfComponents(1);
    datas[i]->SetNumberOfTuples(n);
  }

  template <typename U>
  typename boost::enable_if<is_vector<non_ref_tuple_element<U>>>::type
  operator()(U i) {
    typedef typename mpl::at<mpl_type_vector, U>::type variable_type;
    const char *name = variable_type().name;
    if (name[0] == '_') {
      return;
    }
    datas[i] =
        vtkFloatArray::SafeDownCast(grid->GetPointData()->GetArray(name));
    if (!datas[i]) {
      datas[i] = vtkSmartPointer<vtkFloatArray>::New();
      datas[i]->SetName(name);
    }
    datas[i]->SetNumberOfComponents(non_ref_tuple_element<U>::size);
    datas[i]->SetNumberOfTuples(n);
  }

#ifdef HAVE_EIGEN
  template <typename U>
  typename boost::enable_if<is_eigen_vector<non_ref_tuple_element<U>>>::type
  operator()(U i) {
    typedef typename mpl::at<mpl_type_vector, U>::type variable_type;
    const char *name = variable_type().name;
    if (name[0] == '_') {
      return;
    }
    datas[i] =
        vtkFloatArray::SafeDownCast(grid->GetPointData()->GetArray(name));
    if (!datas[i]) {
      datas[i] = vtkSmartPointer<vtkFloatArray>::New();
      datas[i]->SetName(name);
    }
    datas[i]->SetNumberOfComponents(
        non_ref_tuple_element<U>::RowsAtCompileTime);
    datas[i]->SetNumberOfTuples(n);
  }
#endif

  size_t n;
  vtkSmartPointer<vtkFloatArray> *datas;
  vtkUnstructuredGrid *grid;
};
#endif

template <typename reference> struct setup_datas_for_reading {
  typedef typename reference::mpl_type_vector mpl_type_vector;
  setup_datas_for_reading(size_t n, vtkSmartPointer<vtkFloatArray> *datas,
                          vtkUnstructuredGrid *grid, bool particles)
      : n(n), datas(datas), grid(grid), particles(particles) {}
  template <typename U> void operator()(U i) {
    typedef typename mpl::at<mpl_type_vector, U>::type variable_type;
    const char *name = variable_type().name;
    if (particles) {
      datas[i] =
          vtkFloatArray::SafeDownCast(grid->GetPointData()->GetArray(name));
    } else {
      datas[i] =
          vtkFloatArray::SafeDownCast(grid->GetCellData()->GetArray(name));
    }

    CHECK(datas[i], "No data array " << name << " in vtkUnstructuredGrid");
    CHECK(datas[i]->GetNumberOfTuples() == n,
          "Data array " << name << " has size != id array. data size = "
                        << datas[i]->GetNumberOfTuples()
                        << ". id size = " << n);
  }

  size_t n;
  bool particles;
  vtkSmartPointer<vtkFloatArray> *datas;
  vtkUnstructuredGrid *grid;
};

#endif

} // namespace detail
} // namespace Aboria

#endif
