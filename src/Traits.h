
#ifndef TRAITS_H_
#define TRAITS_H_

#include "CudaInclude.h"
#include "Get.h"
#include "Variable.h"
#include "Vector.h"
#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/permutation_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/vector.hpp>
#include <random>
#include <tuple>
#include <vector>

namespace Aboria {

namespace mpl = boost::mpl;

struct default_traits {
  template <typename T> struct vector_type { typedef std::vector<T> type; };

#ifdef ABORIA_THRUST_USE_THRUST_TUPLE
  template <typename T1 = thrust::null_type, typename T2 = thrust::null_type,
            typename T3 = thrust::null_type, typename T4 = thrust::null_type,
            typename T5 = thrust::null_type, typename T6 = thrust::null_type,
            typename T7 = thrust::null_type, typename T8 = thrust::null_type,
            typename T9 = thrust::null_type>
  struct tuple_type {
    typedef thrust::tuple<T1, T2, T3, T4, T5, T6, T7, T8, T9> type;
  };
  template <std::size_t I, class T> struct tuple_element {
    typedef thrust::tuple_element<I, T> type;
  };
  template <class T> struct tuple_size { typedef thrust::tuple_size<T> type; };
#else
 template <typename... T> struct tuple_type { typedef std::tuple<T...> type; };
 template <std::size_t I, class T> struct tuple_element {
    typedef std::tuple_element<I, T> type;
  };
  template <class T> struct tuple_size { typedef std::tuple_size<T> type; };

#endif

#ifdef HAVE_THRUST
  template <typename ElementIterator, typename IndexIterator>
  static auto make_permutation_iterator(ElementIterator e, IndexIterator i) {
    return thrust::make_permutation_iterator(e, i);
  }

  template <class AdaptableUnaryFunction, class Iterator>
  static auto make_transform_iterator(Iterator it, AdaptableUnaryFunction fun) {
    return thrust::make_transform_iterator(it, fun);
  }

  template <typename... T> static auto make_tuple(T... args) {
    return thrust::make_tuple(args...);
  }

  template <typename IteratorTuple>
  static auto make_zip_iterator(IteratorTuple arg) {
    return thrust::make_zip_iterator(arg);
  }

  template <typename Incrementable>
  static auto make_counting_iterator(Incrementable x) {
    return thrust::make_counting_iterator(x);
  }

#else

  template <typename ElementIterator, typename IndexIterator>
  static auto make_permutation_iterator(ElementIterator e, IndexIterator i) {
    return boost::make_permutation_iterator(e, i);
  }

  template <class AdaptableUnaryFunction, class Iterator>
  static auto make_transform_iterator(Iterator it, AdaptableUnaryFunction fun) {
    return boost::make_transform_iterator(it, fun);
  }

  template <typename... T> static auto make_tuple(T... args) {
    return boost::make_tuple(args...);
  }

  template <typename IteratorTuple>
  static auto make_zip_iterator(IteratorTuple arg) {
    return boost::make_zip_iterator(arg);
  }

  template <typename Incrementable>
  static auto make_counting_iterator(Incrementable x) {
    return boost::make_counting_iterator(x);
  }
#endif
};

template <template <typename, typename> class VECTOR> struct Traits {};

template <> struct Traits<std::vector> : public default_traits {};

#ifdef HAVE_THRUST
template <> struct Traits<thrust::device_vector> : public default_traits {
  template <typename T> struct vector_type {
    typedef thrust::device_vector<T> type;
  };
  template <typename T1 = thrust::null_type, typename T2 = thrust::null_type,
            typename T3 = thrust::null_type, typename T4 = thrust::null_type,
            typename T5 = thrust::null_type, typename T6 = thrust::null_type,
            typename T7 = thrust::null_type, typename T8 = thrust::null_type,
            typename T9 = thrust::null_type>
  struct tuple_type {
    typedef thrust::tuple<T1, T2, T3, T4, T5, T6, T7, T8, T9> type;
  };
  template <std::size_t I, class T> struct tuple_element {
    typedef thrust::tuple_element<I, T> type;
  };
  template <class T> struct tuple_size { typedef thrust::tuple_size<T> type; };

  template <typename ElementIterator, typename IndexIterator>
  static auto make_permutation_iterator(ElementIterator e, IndexIterator i) {
    return thrust::make_permutation_iterator(e, i);
  }

  template <class AdaptableUnaryFunction, class Iterator>
  static auto make_transform_iterator(Iterator it, AdaptableUnaryFunction fun) {
    return thrust::make_transform_iterator(it, fun);
  }

  template <typename... T> static auto make_tuple(T... args) {
    return thrust::make_tuple(args...);
  }

  template <typename IteratorTuple>
  static auto make_zip_iterator(IteratorTuple arg) {
    return thrust::make_zip_iterator(arg);
  }

  template <typename Incrementable>
  static auto make_counting_iterator(Incrementable x) {
    return thrust::make_counting_iterator(x);
  }
};

template <> struct Traits<thrust::host_vector> : public default_traits {
  template <typename T> struct vector_type {
    typedef thrust::host_vector<T> type;
  };
  template <typename T1 = thrust::null_type, typename T2 = thrust::null_type,
            typename T3 = thrust::null_type, typename T4 = thrust::null_type,
            typename T5 = thrust::null_type, typename T6 = thrust::null_type,
            typename T7 = thrust::null_type, typename T8 = thrust::null_type,
            typename T9 = thrust::null_type>
  struct tuple_type {
    typedef thrust::tuple<T1, T2, T3, T4, T5, T6, T7, T8, T9> type;
  };
  template <std::size_t I, class T> struct tuple_element {
    typedef thrust::tuple_element<I, T> type;
  };
  template <class T> struct tuple_size { typedef thrust::tuple_size<T> type; };

  template <typename ElementIterator, typename IndexIterator>
  static auto make_permutation_iterator(ElementIterator e, IndexIterator i) {
    return thrust::make_permutation_iterator(e, i);
  }

  template <class AdaptableUnaryFunction, class Iterator>
  static auto make_transform_iterator(Iterator it, AdaptableUnaryFunction fun) {
    return thrust::make_transform_iterator(it, fun);
  }

  template <typename... T> static auto make_tuple(T... args) {
    return thrust::make_tuple(args...);
  }

  template <typename IteratorTuple>
  static auto make_zip_iterator(IteratorTuple arg) {
    return thrust::make_zip_iterator(arg);
  }

  template <typename Incrementable>
  static auto make_counting_iterator(Incrementable x) {
    return thrust::make_counting_iterator(x);
  }
};
#endif

template <typename ARG, unsigned int DomainD, unsigned int SelfD,
          typename TRAITS>
struct TraitsCommon {
  typedef typename ARG::
      ERROR_FIRST_TEMPLATE_ARGUMENT_TO_PARTICLES_MUST_BE_A_STD_TUPLE_TYPE error;
};

template <typename traits, unsigned int DomainD, unsigned int SelfD,
          typename... TYPES>
struct TraitsCommon<std::tuple<TYPES...>, DomainD, SelfD, traits>
    : public traits {

  template <typename... T>
  using tuple = typename traits::template tuple_type<T...>::type;

  template <size_t I, typename T>
  using tuple_element = typename traits::template tuple_element<I, T>::type;

  template <typename T>
  using vector = typename traits::template vector_type<T>::type;

  // TODO: use vector below

  const static unsigned int dimension = DomainD;
  typedef vector<Vector<double, DomainD>> vector_double_d;
  typedef typename vector_double_d::iterator vector_double_d_iterator;
  typedef
      typename vector_double_d::const_iterator vector_double_d_const_iterator;
  typedef typename traits::template vector_type<Vector<int, DomainD>>::type
      vector_int_d;
  typedef
      typename traits::template vector_type<Vector<unsigned int, DomainD>>::type
          vector_unsigned_int_d;
  typedef
      typename vector_unsigned_int_d::iterator vector_unsigned_int_d_iterator;
  typedef typename vector_unsigned_int_d::const_iterator
      vector_unsigned_int_d_const_iterator;
  typedef typename traits::template vector_type<Vector<bool, DomainD>>::type
      vector_bool_d;
  typedef typename traits::template vector_type<int>::type vector_int;
  typedef typename traits::template vector_type<double>::type vector_double;
  typedef typename traits::template vector_type<size_t>::type vector_size_t;
  typedef typename traits::template vector_type<unsigned int>::type
      vector_unsigned_int;
  typedef typename vector_unsigned_int::iterator vector_unsigned_int_iterator;
  typedef typename vector_unsigned_int::const_iterator
      vector_unsigned_int_const_iterator;
  typedef
      typename traits::template vector_type<Vector<int, 2>>::type vector_int2;

  typedef Vector<double, dimension> double_d;
  typedef Vector<int, dimension> int_d;
  typedef Vector<unsigned int, dimension> unsigned_int_d;
  typedef Vector<bool, dimension> bool_d;

  typedef typename std::conditional<(SelfD > 1), particles_d<SelfD>,
                                    position_d<dimension>>::type position;
  typedef typename position::value_type position_value_type;
  typedef alive::value_type alive_value_type;
  typedef id::value_type id_value_type;
  typedef generator::value_type random_value_type;

  typedef typename traits::template vector_type<position_value_type>::type
      position_vector_type;
  typedef typename traits::template vector_type<alive_value_type>::type
      alive_vector_type;
  typedef
      typename traits::template vector_type<id_value_type>::type id_vector_type;
  typedef typename traits::template vector_type<random_value_type>::type
      random_vector_type;

  typedef traits traits_type;
  typedef mpl::vector<position, id, alive, generator, TYPES...> mpl_type_vector;

  typedef tuple<typename position_vector_type::iterator,
                typename id_vector_type::iterator,
                typename alive_vector_type::iterator,
                typename random_vector_type::iterator,
                typename traits::template vector_type<
                    typename TYPES::value_type>::type::iterator...>
      tuple_of_iterators_type;

  typedef tuple<typename position_vector_type::const_iterator,
                typename id_vector_type::const_iterator,
                typename alive_vector_type::const_iterator,
                typename random_vector_type::const_iterator,
                typename traits::template vector_type<
                    typename TYPES::value_type>::type::const_iterator...>
      tuple_of_const_iterators_type;

  // need a std::tuple here, rather than a thrust one...
  typedef std::tuple<position_vector_type, id_vector_type, alive_vector_type,
                     random_vector_type,
                     typename traits::template vector_type<
                         typename TYPES::value_type>::type...>
      vectors_data_type;

  typedef
      typename Aboria::zip_iterator<tuple_of_iterators_type, mpl_type_vector>
          iterator;
  typedef typename Aboria::zip_iterator<tuple_of_const_iterators_type,
                                        mpl_type_vector>
      const_iterator;

  typedef typename iterator::reference reference;
  typedef typename iterator::value_type value_type;
  typedef typename iterator::pointer pointer;
  typedef typename iterator::getter_raw_pointer raw_pointer;
  typedef typename iterator::getter_raw_reference raw_reference;
  typedef typename const_iterator::getter_raw_reference raw_const_reference;
  typedef typename const_iterator::reference const_reference;
  typedef Aboria::getter_type<vectors_data_type, mpl_type_vector> data_type;

  // need a std::tuple here, rather than a thrust one...
  const static size_t N = std::tuple_size<vectors_data_type>::value;

  template <std::size_t... I>
  static iterator begin_impl(data_type &data, detail::index_sequence<I...>) {
    return iterator(get_by_index<I>(data).begin()...);
  }

  template <std::size_t... I>
  static iterator end_impl(data_type &data, detail::index_sequence<I...>) {
    return iterator(get_by_index<I>(data).end()...);
  }

  template <std::size_t... I>
  static const_iterator cbegin_impl(const data_type &data,
                                    detail::index_sequence<I...>) {
    return const_iterator(get_by_index<I>(data).cbegin()...);
  }

  template <std::size_t... I>
  static const_iterator cend_impl(const data_type &data,
                                  detail::index_sequence<I...>) {
    return const_iterator(get_by_index<I>(data).cend()...);
  }

  template <std::size_t... I>
  static reference index_impl(data_type &data, const size_t i,
                              detail::index_sequence<I...>, std::true_type) {
    return reference(get_by_index<I>(data)[i]...);
  }

  template <std::size_t... I>
  static reference index_impl(data_type &data, const size_t i,
                              detail::index_sequence<I...>, std::false_type) {
    return reference(get_by_index<I>(data)[i]...);
  }

  template <std::size_t... I>
  static const_reference index_const_impl(const data_type &data, const size_t i,
                                          detail::index_sequence<I...>,
                                          std::true_type) {
    return const_reference(get_by_index<I>(data)[i]...);
  }

  template <std::size_t... I>
  static const_reference index_const_impl(const data_type &data, const size_t i,
                                          detail::index_sequence<I...>,
                                          std::false_type) {
    return const_reference(get_by_index<I>(data)[i]...);
  }

  template <std::size_t... I>
  static void clear_impl(data_type &data, detail::index_sequence<I...>) {
    int dummy[] = {0, (get_by_index<I>(data).clear(), void(), 0)...};
    static_cast<void>(dummy);
  }

  template <std::size_t... I>
  static void resize_impl(data_type &data, const size_t new_size,
                          detail::index_sequence<I...>) {
    int dummy[] = {0, (get_by_index<I>(data).resize(new_size), void(), 0)...};
    static_cast<void>(dummy);
  }

  template <std::size_t... I>
  static void push_back_impl(data_type &data, const value_type &val,
                             detail::index_sequence<I...>) {
    int dummy[] = {0, (get_by_index<I>(data).push_back(get_by_index<I>(val)),
                       void(), 0)...};
    // int dummy[] = { 0,
    // (detail::get_impl<I>(data.get_tuple()).push_back(detail::get_impl<I>(val.get_tuple())),void(),0)...
    // };
    static_cast<void>(dummy); // Avoid warning for unused variable.
  }

  template <std::size_t... I>
  static void pop_back_impl(data_type &data, detail::index_sequence<I...>) {
    int dummy[] = {0, (get_by_index<I>(data).pop_back(), void(), 0)...};
    static_cast<void>(dummy); // Avoid warning for unused variable.
  }

  template <std::size_t... I>
  static iterator erase_impl(data_type &data, iterator position,
                             detail::index_sequence<I...>) {
    return iterator(get_by_index<I>(data).erase(get_by_index<I>(position))...);
  }

  template <std::size_t... I>
  static iterator erase_impl(data_type &data, iterator first, iterator last,
                             detail::index_sequence<I...>) {
    return iterator(get_by_index<I>(data).erase(get_by_index<I>(first),
                                                get_by_index<I>(last))...);
  }

  template <std::size_t... I>
  static iterator insert_impl(data_type &data, iterator position,
                              const value_type &val,
                              detail::index_sequence<I...>) {
    return iterator(get_by_index<I>(data).insert(get_by_index<I>(position),
                                                 get_by_index<I>(val))...);
  }

  template <std::size_t... I>
  static void insert_impl(data_type &data, iterator position, size_t n,
                          const value_type &val, detail::index_sequence<I...>) {
    int dummy[] = {0, (get_by_index<I>(data).insert(position, n,
                                                    get_by_index<I>(val)))...};
    static_cast<void>(dummy);
  }

  template <class InputIterator, std::size_t... I>
  static iterator insert_impl(data_type &data, iterator position,
                              InputIterator first, InputIterator last,
                              detail::index_sequence<I...>) {
    return iterator(get_by_index<I>(data).insert(get_by_index<I>(position),
                                                 get_by_index<I>(first),
                                                 get_by_index<I>(last))...);
  }

  template <class InputIterator, std::size_t... I>
  static data_type construct_impl(InputIterator first, InputIterator last,
                                  detail::index_sequence<I...>) {
    return data_type(typename tuple_element<I, vectors_data_type>::type(
        get_by_index<I>(first), get_by_index<I>(last))...);
  }

  template <std::size_t... I>
  static void header_to_stream_impl(std::ostream &os,
                                    detail::index_sequence<I...>) {
    int dummy[] = {
        0, (void(os << (I == 0 ? "" : ", ") <<
                 typename mpl::at<mpl_type_vector, mpl::int_<I>>::type().name),
            0)...};
    static_cast<void>(dummy); // Avoid warning for unused variable.
  }

  template <class InputIterator, std::size_t... I>
  static void to_stream_impl(InputIterator i, std::ostream &os,
                             detail::index_sequence<I...>) {
    int dummy[] = {
        0, (void(os << (I == 0 ? "" : ", ") << *get_by_index<I>(i)), 0)...};
    static_cast<void>(dummy); // Avoid warning for unused variable.
  }

  template <class InputIterator, std::size_t... I>
  static void from_stream_impl(InputIterator i, std::istream &is,
                               detail::index_sequence<I...>) {
    int dummy[] = {0, (is >> get_by_index<I>(i))...};
    static_cast<void>(dummy); // Avoid warning for unused variable.
  }

  template <class Archive, std::size_t... I>
  static void serialize_impl(const data_type &data, Archive &ar,
                             const unsigned int version,
                             detail::index_sequence<I...>) {
    int dummy[] = {
        0,
        (void(ar &
              boost::serialization::make_nvp(
                  typename mpl::at<mpl_type_vector, mpl::int_<I>>::type().name,
                  get_by_index<I>(data))),
         0)...};
    static_cast<void>(dummy); // Avoid warning for unused variable.
  }

  template <typename Indices = detail::make_index_sequence<N>>
  static iterator begin(data_type &data) {
    return begin_impl(data, Indices());
  }

  template <typename Indices = detail::make_index_sequence<N>>
  static iterator end(data_type &data) {
    return end_impl(data, Indices());
  }

  template <typename Indices = detail::make_index_sequence<N>>
  static const_iterator cbegin(const data_type &data) {
    return cbegin_impl(data, Indices());
  }

  template <typename Indices = detail::make_index_sequence<N>>
  static const_iterator cend(const data_type &data) {
    return cend_impl(data, Indices());
  }

  template <typename Indices = detail::make_index_sequence<N>>
  static reference index(data_type &data, const size_t i) {
    return index_impl(data, i, Indices(),
                      std::is_reference<decltype(get<id>(data)[0])>());
  }

  template <typename Indices = detail::make_index_sequence<N>>
  static const_reference index(const data_type &data, const size_t i) {
    return index_const_impl(data, i, Indices(),
                            std::is_reference<decltype(get<id>(data)[0])>());
  }

  template <typename Indices = detail::make_index_sequence<N>>
  static void clear(data_type &data) {
    clear_impl(data, Indices());
  }

  template <typename Indices = detail::make_index_sequence<N>>
  static void resize(data_type &data, const size_t new_size) {
    resize_impl(data, new_size, Indices());
  }

  template <typename Indices = detail::make_index_sequence<N>>
  static void push_back(data_type &data, const value_type &val) {
    push_back_impl(data, val, Indices());
  }

  template <typename Indices = detail::make_index_sequence<N>>
  static void pop_back(data_type &data) {
    pop_back_impl(data, Indices());
  }

  template <typename Indices = detail::make_index_sequence<N>>
  static iterator erase(data_type &data, iterator pos) {
    return erase_impl(data, pos, Indices());
  }

  template <typename Indices = detail::make_index_sequence<N>>
  static iterator erase(data_type &data, iterator first, iterator last) {
    return erase_impl(data, first, last, Indices());
  }

  template <typename Indices = detail::make_index_sequence<N>>
  static iterator insert(data_type &data, iterator pos, const value_type &val) {
    return insert_impl(data, pos, val, Indices());
  }

  template <typename Indices = detail::make_index_sequence<N>>
  static void insert(data_type &data, iterator position, size_t n,
                     const value_type &val) {
    insert_impl(data, position, val, n, Indices());
  }

  template <class InputIterator,
            typename Indices = detail::make_index_sequence<N>>
  static iterator insert(data_type &data, iterator pos, InputIterator first,
                         InputIterator last) {
    return insert_impl(data, pos, first, last, Indices());
  }

  template <class InputIterator,
            typename Indices = detail::make_index_sequence<N>>
  static data_type construct(InputIterator first, InputIterator last) {
    return construct_impl(first, last, Indices());
  }

  template <typename Indices = detail::make_index_sequence<N>>
  static void header_to_stream(std::ostream &os) {
    header_to_stream_impl(os, Indices());
  }

  template <class InputIterator,
            typename Indices = detail::make_index_sequence<N>>
  static void to_stream(InputIterator i, std::ostream &os) {
    to_stream_impl(i, os, Indices());
  }

  template <class InputIterator,
            typename Indices = detail::make_index_sequence<N>>
  static void from_stream(InputIterator i, std::istream &is) {
    from_stream_impl(i, is, Indices());
  }

  template <class Archive, typename Indices = detail::make_index_sequence<N>>
  static void serialize(const data_type &data, Archive &ar,
                        const unsigned int version) {
    serialize_impl(data, ar, version, Indices());
  }

  typedef typename position_vector_type::size_type size_type;
  typedef typename position_vector_type::difference_type difference_type;
};

} // namespace Aboria

#endif // TRAITS_H_
