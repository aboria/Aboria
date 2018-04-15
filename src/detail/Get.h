
#ifndef GET_DETAIL_H_
#define GET_DETAIL_H_

#include <algorithm>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/vector.hpp>
#include <tuple>
#include <type_traits>

#include "CudaInclude.h"
#include "Get.h"
#include "Log.h"

#if HAVE_THRUST
namespace thrust {

template <typename mpl_vector_type, typename tuple_type>
struct iterator_system<Aboria::zip_iterator<tuple_type, mpl_vector_type>> {};

template <typename mpl_vector_type, typename T0, typename... T>
struct iterator_system<
    Aboria::zip_iterator<std::tuple<T0, T...>, mpl_vector_type>> {
  typedef typename iterator_system<T0>::type type;
};

template <typename mpl_vector_type, typename T0, typename T1, typename T2,
          typename T3, typename T4, typename T5, typename T6, typename T7,
          typename T8, typename T9>
struct iterator_system<Aboria::zip_iterator<
    thrust::tuple<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9>, mpl_vector_type>> {
  typedef typename iterator_system<T0>::type type;
};

/*
template <typename TUPLE, typename mpl_vector_type>
struct getter_type;
*/

// what follows is a copy of thrust's detail/raw_reference_cast.h for Aboria's
// getter type
namespace detail {

/*
template <>
struct pointer_traits<thrust::null_type> {
typedef null_type raw_pointer;
};
*/

// specialize is_tuple_of_iterator_references to for getter_type to
// device_reference
template <typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6, typename T7, typename T8, typename T9, typename T10,
          typename MplVector>
struct is_tuple_of_iterator_references<
    Aboria::getter_type<thrust::tuple<
                            // this seems dangerous, matches values as well,
                            // swap to something like below?
                            T1, T2, T3, T4, T5, T6, T7, T8, T9, T10>,
                        MplVector>
    /*
    thrust::device_reference<T1>,
    thrust::device_reference<T2>,
    thrust::device_reference<T3>,
    thrust::device_reference<T4>,
    thrust::device_reference<T5>,
    thrust::device_reference<T6>,
    thrust::device_reference<T7>,
    thrust::device_reference<T8>,
    thrust::device_reference<T9>,
    thrust::device_reference<T10>>, MplVector>
    */
    > : thrust::detail::true_type {};

template <typename MplVector, typename... T>
struct is_tuple_of_iterator_references<
    Aboria::getter_type<std::tuple<
                            // this seems dangerous, matches values as well,
                            // swap to something like below?
                            T...>,
                        MplVector>> : thrust::detail::true_type {};

// specialize is_unwrappable
template <typename TUPLE, typename mpl_vector_type>
struct is_unwrappable<Aboria::getter_type<TUPLE, mpl_vector_type>>
    : is_unwrappable<TUPLE> {};

namespace raw_reference_detail {

// recurse on tuples
template <typename mpl_vector_type, typename T0, typename T1, typename T2,
          typename T3, typename T4, typename T5, typename T6, typename T7,
          typename T8, typename T9>
struct raw_reference_tuple_helper<Aboria::getter_type<
    thrust::tuple<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9>, mpl_vector_type>> {
  typedef Aboria::getter_type<
      thrust::tuple<typename raw_reference_tuple_helper<T0>::type,
                    typename raw_reference_tuple_helper<T1>::type,
                    typename raw_reference_tuple_helper<T2>::type,
                    typename raw_reference_tuple_helper<T3>::type,
                    typename raw_reference_tuple_helper<T4>::type,
                    typename raw_reference_tuple_helper<T5>::type,
                    typename raw_reference_tuple_helper<T6>::type,
                    typename raw_reference_tuple_helper<T7>::type,
                    typename raw_reference_tuple_helper<T8>::type,
                    typename raw_reference_tuple_helper<T9>::type>,
      mpl_vector_type>
      type;
};

template <typename mpl_vector_type, typename... T>
struct raw_reference_tuple_helper<
    Aboria::getter_type<std::tuple<T...>, mpl_vector_type>> {
  typedef Aboria::getter_type<
      std::tuple<typename raw_reference_tuple_helper<T>::type...>,
      mpl_vector_type>
      type;
};

} // namespace raw_reference_detail

template <typename TUPLE, typename mpl_vector_type>
struct raw_reference<Aboria::getter_type<TUPLE, mpl_vector_type>> {
private:
  typedef TUPLE tuple_type;

public:
  typedef typename eval_if<
      is_unwrappable<tuple_type>::value,
      identity_<
          Aboria::getter_type<typename raw_reference_detail::
                                  raw_reference_tuple_helper<tuple_type>::type,
                              mpl_vector_type>>,
      add_reference<Aboria::getter_type<TUPLE, mpl_vector_type>>>::type type;
};

} // namespace detail

template <typename TUPLE, typename mpl_vector_type>
__host__ __device__ typename detail::enable_if_unwrappable<
    Aboria::getter_type<TUPLE, mpl_vector_type>,
    typename detail::raw_reference<
        Aboria::getter_type<TUPLE, mpl_vector_type>>::type>::type
raw_reference_cast(Aboria::getter_type<TUPLE, mpl_vector_type> t);

namespace detail {
namespace aboria_addition {

struct raw_reference_caster {
  template <typename T>
  __host__ __device__ typename detail::raw_reference<T>::type
  operator()(T &ref) {
    return thrust::raw_reference_cast(ref);
  }

  template <typename T>
  __host__ __device__ typename detail::raw_reference<const T>::type
  operator()(const T &ref) {
    return thrust::raw_reference_cast(ref);
  }

  template <typename TUPLE, typename mpl_vector_type>
  __host__ __device__ typename detail::raw_reference<
      Aboria::getter_type<TUPLE, mpl_vector_type>>::type
  operator()(
      Aboria::getter_type<TUPLE, mpl_vector_type> t,
      typename enable_if<is_unwrappable<
          Aboria::getter_type<TUPLE, mpl_vector_type>>::value>::type * = 0) {
    return thrust::raw_reference_cast(t);
  }
}; // end raw_reference_caster

} // namespace aboria_addition
} // namespace detail

template <typename TUPLE, typename mpl_vector_type>
__host__ __device__ typename detail::enable_if_unwrappable<
    Aboria::getter_type<TUPLE, mpl_vector_type>,
    typename detail::raw_reference<
        Aboria::getter_type<TUPLE, mpl_vector_type>>::type>::type
raw_reference_cast(Aboria::getter_type<TUPLE, mpl_vector_type> t) {
  typedef typename detail::raw_reference<
      Aboria::getter_type<TUPLE, mpl_vector_type>>::type return_type;

  thrust::detail::raw_reference_caster f;

  // note that we pass raw_reference_tuple_helper, not raw_reference as the
  // unary metafunction the different way that raw_reference_tuple_helper
  // unwraps tuples is important
  return return_type(thrust::detail::tuple_host_device_transform<
                     detail::raw_reference_detail::raw_reference_tuple_helper>(
      t.get_tuple(), f));
} // end raw_reference_cast

} // namespace thrust
#endif

namespace Aboria {

namespace detail {

template <typename T> struct is_std_getter_type {
  const static bool value = false;
  typedef std::false_type type;
};

template <typename M, typename... T>
struct is_std_getter_type<getter_type<std::tuple<T...>, M>> {
  const static bool value = true;
  typedef std::true_type type;
};

template <typename T> struct is_thrust_getter_type {
  const static bool value = false;
  typedef std::false_type type;
};

#ifdef HAVE_THRUST
template <typename M, typename... T>
struct is_thrust_getter_type<getter_type<thrust::tuple<T...>, M>> {
  const static bool value = true;
  typedef std::true_type type;
};
#endif

template <typename T> struct is_std_zip_iterator {
  const static bool value = false;
  typedef std::false_type type;
};

template <typename M, typename... T>
struct is_std_zip_iterator<zip_iterator<std::tuple<T...>, M>> {
  const static bool value = true;
  typedef std::true_type type;
};

template <typename T> struct is_thrust_zip_iterator {
  const static bool value = false;
  typedef std::false_type type;
};

#ifdef HAVE_THRUST
template <typename M, typename... T>
struct is_thrust_zip_iterator<zip_iterator<thrust::tuple<T...>, M>> {
  const static bool value = true;
  typedef std::true_type type;
};
#endif

ABORIA_HOST_DEVICE_IGNORE_WARN
template <size_t I, typename... T>
CUDA_HOST_DEVICE const typename std::tuple_element<I, std::tuple<T...>>::type &
get_impl(const std::tuple<T...> &arg) {
#if defined(__CUDA_ARCH__)
  ERROR_CUDA("Cannot use get_impl on `std::tuple` in device code");
  return std::get<I>(arg);
#else
  return std::get<I>(arg);
#endif
}

ABORIA_HOST_DEVICE_IGNORE_WARN
template <size_t I, typename... T>
CUDA_HOST_DEVICE typename std::tuple_element<I, std::tuple<T...>>::type &
get_impl(std::tuple<T...> &arg) {
#if defined(__CUDA_ARCH__)
  ERROR_CUDA("Cannot use get_impl on `std::tuple` in device code");
  return std::get<I>(arg);
#else
  return std::get<I>(arg);
#endif
}

ABORIA_HOST_DEVICE_IGNORE_WARN
template <size_t I, typename... T>
CUDA_HOST_DEVICE typename std::tuple_element<I, std::tuple<T...>>::type &&
get_impl(std::tuple<T...> &&arg) {
#if defined(__CUDA_ARCH__)
  ERROR_CUDA("Cannot use get_impl on `std::tuple` in device code");
  return std::get<I>(arg);
#else
  return std::get<I>(arg);
#endif
}

#ifdef HAVE_THRUST
/*
ABORIA_HOST_DEVICE_IGNORE_WARN
template<size_t I, typename ... T>
typename thrust::tuple_element<I,thrust::tuple<T...>>::type const &
get_impl(const thrust::tuple<T...>& arg) {
    return thrust::get<I>(arg);
}

ABORIA_HOST_DEVICE_IGNORE_WARN
template<size_t I, typename ... T>
typename thrust::tuple_element<I,thrust::tuple<T...>>::type &
get_impl(thrust::tuple<T...>& arg) {
    return thrust::get<I>(arg);
}

ABORIA_HOST_DEVICE_IGNORE_WARN
template<size_t I, typename ... T>
typename thrust::tuple_element<I,thrust::tuple<T...>>::type &
get_impl(thrust::tuple<T...>&& arg) {
    return thrust::get<I>(arg);
}
*/

// why do I need this? why doesnt the above work... works on a simpler test
// case...
ABORIA_HOST_DEVICE_IGNORE_WARN
template <size_t I, typename T0, typename T1, typename T2, typename T3,
          typename T4, typename T5, typename T6, typename T7, typename T8,
          typename T9>
CUDA_HOST_DEVICE typename thrust::tuple_element<
    I, thrust::tuple<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9>>::type const &
get_impl(const thrust::tuple<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9> &arg) {
  return thrust::get<I>(arg);
}

ABORIA_HOST_DEVICE_IGNORE_WARN
template <size_t I, typename T0, typename T1, typename T2, typename T3,
          typename T4, typename T5, typename T6, typename T7, typename T8,
          typename T9>
CUDA_HOST_DEVICE typename thrust::tuple_element<
    I, thrust::tuple<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9>>::type &
get_impl(thrust::tuple<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9> &arg) {
  return thrust::get<I>(arg);
}

ABORIA_HOST_DEVICE_IGNORE_WARN
template <size_t I, typename T0, typename T1, typename T2, typename T3,
          typename T4, typename T5, typename T6, typename T7, typename T8,
          typename T9>
CUDA_HOST_DEVICE typename thrust::tuple_element<
    I, thrust::tuple<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9>>::type &
get_impl(thrust::tuple<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9> &&arg) {
  return thrust::get<I>(arg);
}

#endif

template <typename IteratorTuple, typename MplVector>
struct zip_iterator_base {};

template <typename mpl_vector_type, typename... Types>
struct zip_iterator_base<std::tuple<Types...>, mpl_vector_type> {
  typedef std::tuple<Types...> iterator_tuple_type;

  typedef getter_type<
      typename zip_helper<iterator_tuple_type>::tuple_value_type,
      mpl_vector_type>
      value_type;
  typedef getter_type<typename zip_helper<iterator_tuple_type>::tuple_reference,
                      mpl_vector_type>
      reference;

public:
  typedef boost::iterator_facade<
      zip_iterator<iterator_tuple_type, mpl_vector_type>, value_type,
      typename zip_helper<iterator_tuple_type>::iterator_category, reference,
      typename zip_helper<iterator_tuple_type>::difference_type>
      type;

}; // end zip_iterator_base

#ifdef HAVE_THRUST
template <typename mpl_vector_type, typename... Types>
struct zip_iterator_base<thrust::tuple<Types...>, mpl_vector_type> {
  typedef thrust::tuple<Types...> iterator_tuple_type;

  typedef getter_type<
      typename zip_helper<iterator_tuple_type>::tuple_value_type,
      mpl_vector_type>
      value_type;
  typedef getter_type<typename zip_helper<iterator_tuple_type>::tuple_reference,
                      mpl_vector_type>
      reference;

public:
  typedef thrust::iterator_facade<
      zip_iterator<iterator_tuple_type, mpl_vector_type>, value_type,
      typename zip_helper<iterator_tuple_type>::system,
      typename zip_helper<iterator_tuple_type>::iterator_category, reference,
      typename zip_helper<iterator_tuple_type>::difference_type>
      type;

}; // end zip_iterator_base
#endif

/*
__aboria_hd_warning_disable__
template<typename pointer, typename tuple_type, std::size_t... I>
CUDA_HOST_DEVICE
static pointer make_pointer(tuple_type&& tuple, index_sequence<I...>) {
    return pointer(&(get_impl<I>(std::forward(tuple)))...);
}
*/

template <typename iterator_tuple_type, typename mpl_vector_type>
typename zip_iterator<iterator_tuple_type, mpl_vector_type>::getter_raw_pointer
iterator_to_raw_pointer(
    const zip_iterator<iterator_tuple_type, mpl_vector_type> &arg,
    std::true_type) {
  typedef typename zip_helper<iterator_tuple_type>::index_type index_type;
  typedef typename zip_iterator<iterator_tuple_type,
                                mpl_vector_type>::getter_raw_pointer
      getter_raw_pointer;
  return getter_raw_pointer(zip_helper<iterator_tuple_type>::make_raw_pointer(
      arg.get_tuple(), index_type()));
}

template <typename Iterator>
typename std::iterator_traits<Iterator>::value_type *
single_iterator_to_raw_pointer(const Iterator &arg,
                               std::random_access_iterator_tag) {
  return &*arg;
}

#ifdef HAVE_THRUST
template <typename Iterator>
typename std::iterator_traits<Iterator>::value_type *
single_iterator_to_raw_pointer(const Iterator &arg,
                               thrust::random_access_device_iterator_tag) {
  return thrust::raw_pointer_cast(&*arg);
}
#endif

template <typename Iterator>
typename std::iterator_traits<Iterator>::value_type *
iterator_to_raw_pointer(const Iterator &arg, std::false_type) {
#ifdef HAVE_THRUST
  return thrust::raw_pointer_cast(&*arg);
#else
  return &*arg;
#endif
  // return single_iterator_to_raw_pointer(arg,
  //        typename std::iterator_traits<Iterator>::iterator_category());
}

template <typename T> struct is_zip_iterator {
  typedef std::false_type type;
  static const bool value = false;
};

template <typename tuple_type, typename mpl_vector_type>
struct is_zip_iterator<zip_iterator<tuple_type, mpl_vector_type>> {
  typedef std::true_type type;
  static const bool value = true;
};

} // namespace detail
} // namespace Aboria

#endif // GET_DETAIL_H_
