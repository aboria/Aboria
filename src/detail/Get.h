
#ifndef GET_DETAIL_H_ 
#define GET_DETAIL_H_ 

#include <boost/mpl/vector.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <algorithm>
#include <tuple>
#include <type_traits>

#include "Log.h"


#if defined(__CUDACC__)
namespace thrust {
    template <>
    struct iterator_traits<thrust::null_type> {
        typedef thrust::null_type value_type;
        typedef thrust::null_type reference;
        typedef thrust::null_type pointer;
    };

    template <typename mpl_vector_type, typename T0, typename ... T>
    struct iterator_system<Aboria::zip_iterator<tuple_ns::tuple<T0,T...>,mpl_vector_type>> {
        typedef typename iterator_system<T0>::type type;
    };

template <typename TUPLE, typename mpl_vector_type> 
struct getter_type;

// what follows is a copy of thrust's detail/raw_reference_cast.h for Aboria's getter type
/*
namespace detail {

// specialize is_unwrappable
template <typename TUPLE, typename mpl_vector_type> 
  struct is_unwrappable<Aboria::getter_type<TUPLE,mpl_vector_type>>
    : is_unwrappable<TUPLE>
    {};

namespace raw_reference_detail
{


// recurse on tuples
template <typename mpl_vector_type, typename ... T> 
struct raw_reference_tuple_helper<
    Aboria::getter_type<tuple_ns::tuple<T ...>,mpl_vector_type>
    > {
  typedef Aboria::getter_type<
        tuple_ns::tuple<typename raw_reference_tuple_helper<T>::type ...>
        ,mpl_vector_type> type;
};

} //namespace raw_reference_detail

template <typename TUPLE, typename mpl_vector_type> 
struct raw_reference<Aboria::getter_type<TUPLE,mpl_vector_type>> {
  private:
    typedef TUPLE tuple_type;

  public:
    typedef typename eval_if<
      is_unwrappable<tuple_type>::value,
      raw_reference_detail::raw_reference_tuple_helper<tuple_type>,
      add_reference<Aboria::getter_type<TUPLE,mpl_vector_type>>
    >::type type;
};


} //namespace detail

template <typename TUPLE, typename mpl_vector_type> 
__host__ __device__
typename detail::enable_if_unwrappable<
  Aboria::getter_type<TUPLE,mpl_vector_type>,
  typename detail::raw_reference<
    Aboria::getter_type<TUPLE,mpl_vector_type>
  >::type
>::type
raw_reference_cast(Aboria::getter_type<TUPLE,mpl_vector_type> t);

namespace detail  {
namespace aboria_addition {

struct raw_reference_caster
{
  template<typename T>
  __host__ __device__
  typename detail::raw_reference<T>::type operator()(T &ref)
  {
    return thrust::raw_reference_cast(ref);
  }

  template<typename T>
  __host__ __device__
  typename detail::raw_reference<const T>::type operator()(const T &ref)
  {
    return thrust::raw_reference_cast(ref);
  }


  template <typename TUPLE, typename mpl_vector_type> 
  __host__ __device__
  typename detail::raw_reference<
    Aboria::getter_type<TUPLE,mpl_vector_type>
  >::type
  operator()(Aboria::getter_type<TUPLE,mpl_vector_type> t,
             typename enable_if<
               is_unwrappable<Aboria::getter_type<TUPLE,mpl_vector_type>>::value
             >::type * = 0)
  {
    return thrust::raw_reference_cast(t);
  }
}; // end raw_reference_caster


} //namespace aboria_addition
} //namespace detail

template <typename TUPLE, typename mpl_vector_type> 
__host__ __device__
typename detail::enable_if_unwrappable<
  Aboria::getter_type<TUPLE,mpl_vector_type>,
  typename detail::raw_reference<
    Aboria::getter_type<TUPLE,mpl_vector_type>
  >::type
>::type
raw_reference_cast(Aboria::getter_type<TUPLE,mpl_vector_type> t)
{
  thrust::detail:aboria_addition::raw_reference_caster f;

  // note that we pass raw_reference_tuple_helper, not raw_reference as the unary metafunction
  // the different way that raw_reference_tuple_helper unwraps tuples is important
  return thrust::detail::tuple_host_device_transform<detail::raw_reference_detail::raw_reference_tuple_helper>(t, f);
} // end raw_reference_cast

*/

} //namespace thrust
#endif

namespace Aboria {


namespace detail {

namespace mpl = boost::mpl;

// implementation of c++11 make_index_sequence 
// copied from http://stackoverflow.com/questions/17424477/implementation-c14-make-integer-sequence

template<size_t...> struct index_sequence { using type = index_sequence; };
template<typename T1, typename T2> struct concat;
template<size_t... I1, size_t... I2> struct concat<index_sequence<I1...>, index_sequence<I2...>>: index_sequence<I1..., (sizeof...(I1) + I2)...> {};

template<size_t N> struct make_index_sequence;
template<size_t N> struct make_index_sequence: concat<typename make_index_sequence<N/2>::type, typename make_index_sequence<N-N/2>::type>::type {};
template <> struct make_index_sequence<0>: index_sequence<>{};
template <> struct make_index_sequence<1>: index_sequence<0>{};

/// 
/// helper class to find an element of mpl_type_vector from a Variable type T
template<typename T, typename mpl_type_vector>
struct get_elem_by_type {
    typedef T type;
    typedef typename T::value_type value_type;

    /// 
    /// iter is a boost mpl iterator to the found element in Variable T
    typedef typename mpl::find<mpl_type_vector,T>::type iter;

    /// 
    /// index contains the index of the found element
    static const size_t index = iter::pos::value;
    static_assert(index < mpl::size<mpl_type_vector>::type::value,"variable not found in particle");
};


/// helper class to find an element of mpl_type_vector from an
/// unsigned int index I
template<unsigned int I, typename mpl_type_vector>
struct get_elem_by_index {
    BOOST_MPL_ASSERT_RELATION( (mpl::size<mpl_type_vector>::type::value), >, I );
    typedef typename mpl::at<mpl_type_vector,mpl::int_<I> > type;

    /// 
    /// value_type is the variable's value_type at index I 
    typedef typename type::value_type value_type;
    static const size_t index = I;
};


template<typename T>
struct remove_pointer_or_reference_for_null_type {
    typedef T type;
};

#if defined(__CUDACC__)
template<>
struct remove_pointer_or_reference_for_null_type<thrust::null_type*> {
    typedef thrust::null_type type;
};

template<>
struct remove_pointer_or_reference_for_null_type<thrust::null_type&> {
    typedef thrust::null_type type;
};
#endif

template <typename T>
struct is_std_getter_type {
    const static bool value = false;
    typedef std::false_type type;
};

template <typename M, typename ... T>
struct is_std_getter_type<getter_type<std::tuple<T...>,M>> {
    const static bool value = true;
    typedef std::true_type type;
};

template <typename T>
struct is_thrust_getter_type {
    const static bool value = false;
    typedef std::false_type type;
};

template <typename M, typename ... T>
struct is_thrust_getter_type<getter_type<thrust::tuple<T...>,M>> {
    const static bool value = true;
    typedef std::true_type type;
};

template<size_t I, typename ... T>
__aboria_hd_warning_disable__ 
CUDA_HOST_DEVICE
typename std::tuple_element<I,std::tuple<T...>>::type const &
get_impl(const std::tuple<T...>& arg)
{
#if defined(__CUDA_ARCH__)
    ERROR_CUDA("Cannot use get_impl on `std::tuple` in device code");
    return std::get<I>(arg);
#else
    return std::get<I>(arg);
#endif
}

template<size_t I, typename ... T>
__aboria_hd_warning_disable__ 
CUDA_HOST_DEVICE
typename std::tuple_element<I,std::tuple<T...>>::type &
get_impl(std::tuple<T...>& arg)
{
#if defined(__CUDA_ARCH__)
    ERROR_CUDA("Cannot use get_impl on `std::tuple` in device code");
    return std::get<I>(arg);
#else
    return std::get<I>(arg);
#endif
}

template<size_t I, typename ... T>
__aboria_hd_warning_disable__ 
CUDA_HOST_DEVICE
typename std::tuple_element<I,std::tuple<T...>>::type &
get_impl(std::tuple<T...>&& arg)
{
#if defined(__CUDA_ARCH__)
    ERROR_CUDA("Cannot use get_impl on `std::tuple` in device code");
    return std::get<I>(arg);
#else
    return std::get<I>(arg);
#endif
}

template<size_t I, typename ... T>
CUDA_HOST_DEVICE
typename thrust::tuple_element<I,thrust::tuple<T...>>::type const &
get_impl(const thrust::tuple<T...>& arg) {
    return thrust::get<I>(arg);
}

template<size_t I, typename ... T>
CUDA_HOST_DEVICE
typename thrust::tuple_element<I,thrust::tuple<T...>>::type &
get_impl(thrust::tuple<T...>& arg) {
    return thrust::get<I>(arg);
}

template<size_t I, typename ... T>
CUDA_HOST_DEVICE
typename thrust::tuple_element<I,thrust::tuple<T...>>::type &
get_impl(thrust::tuple<T...>&& arg) {
    return thrust::get<I>(arg);
}


template<typename tuple_of_iterators>
struct zip_helper {};


template <typename ... T>
struct zip_helper<std::tuple<T ...>> {
    //typedef std::false_type is_thrust;
    typedef std::tuple<T...> tuple_iterator_type; 
    typedef std::tuple<typename std::iterator_traits<T>::value_type ...> tuple_value_type; 
    typedef std::tuple<typename std::iterator_traits<T>::reference ...> tuple_reference; 
    typedef std::tuple<typename std::iterator_traits<T>::pointer...> tuple_pointer; 

#if defined(__CUDACC__)
    typedef thrust::tuple<
        typename detail::remove_pointer_or_reference_for_null_type<
            typename std::iterator_traits<T>::value_type*>::type...
        > tuple_raw_pointer; 
#else
typedef std::tuple<
        typename detail::remove_pointer_or_reference_for_null_type<
            typename std::iterator_traits<T>::value_type*>::type...
        > tuple_raw_pointer; 
#endif


    typedef std::tuple<
        typename detail::remove_pointer_or_reference_for_null_type<
            typename std::iterator_traits<T>::value_type&>::type...
        > tuple_raw_reference; 
    typedef typename std::tuple<T...> iterator_tuple_type;

    template <unsigned int N>
    using tuple_element = std::tuple_element<N,iterator_tuple_type>;

    typedef typename std::iterator_traits<typename tuple_element<0>::type>::difference_type difference_type;
    typedef typename std::iterator_traits<typename tuple_element<0>::type>::iterator_category iterator_category;
    typedef make_index_sequence<std::tuple_size<iterator_tuple_type>::value> index_type;

};

template <typename ... T>
struct zip_helper<thrust::tuple<T ...>> {
    //typedef std::false_type is_thrust;
    typedef thrust::tuple<T...> tuple_iterator_type; 
    typedef thrust::tuple<typename thrust::iterator_traits<T>::value_type ...> tuple_value_type; 
    typedef thrust::tuple<typename thrust::iterator_traits<T>::reference ...> tuple_reference; 
    typedef thrust::tuple<typename thrust::iterator_traits<T>::pointer...> tuple_pointer; 

    typedef thrust::tuple<
        typename detail::remove_pointer_or_reference_for_null_type<
            typename thrust::iterator_traits<T>::value_type*>::type...
        > tuple_raw_pointer; 
    typedef thrust::tuple<
        typename detail::remove_pointer_or_reference_for_null_type<
            typename thrust::iterator_traits<T>::value_type&>::type...
        > tuple_raw_reference; 
    typedef typename thrust::tuple<T...> iterator_tuple_type;

    template <unsigned int N>
    using tuple_element = thrust::tuple_element<N,iterator_tuple_type>;

    typedef typename thrust::iterator_traits<typename tuple_element<0>::type>::difference_type difference_type;
    typedef typename thrust::iterator_traits<typename tuple_element<0>::type>::iterator_category iterator_category;
    typedef typename thrust::iterator_system<typename tuple_element<0>::type> system;
    typedef make_index_sequence<thrust::tuple_size<iterator_tuple_type>::value> index_type;

};

template<typename IteratorTuple, typename MplVector>
struct zip_iterator_base {};

template <typename mpl_vector_type, typename ... Types> 
struct zip_iterator_base<std::tuple<Types...>, mpl_vector_type>{
    typedef std::tuple<Types...> iterator_tuple_type;
 
    typedef getter_type<typename zip_helper<iterator_tuple_type>::tuple_value_type,mpl_vector_type> value_type;
    typedef getter_type<typename zip_helper<iterator_tuple_type>::tuple_reference,mpl_vector_type> reference;
  
 public:
  
typedef boost::iterator_facade<
    zip_iterator<iterator_tuple_type,mpl_vector_type>,
    value_type,  
    typename zip_helper<iterator_tuple_type>::iterator_category,
    reference,
    typename zip_helper<iterator_tuple_type>::difference_type
> type;

}; // end zip_iterator_base

#ifdef HAVE_THRUST
template <typename mpl_vector_type, typename ... Types> 
struct zip_iterator_base<thrust::tuple<Types...>, mpl_vector_type>{
    typedef thrust::tuple<Types...> iterator_tuple_type;
 
    typedef getter_type<typename zip_helper<iterator_tuple_type>::tuple_value_type,mpl_vector_type> value_type;
    typedef getter_type<typename zip_helper<iterator_tuple_type>::tuple_reference,mpl_vector_type> reference;
  
 public:
  
typedef iterator_facade_ns::iterator_facade<
        zip_iterator<iterator_tuple_type,mpl_vector_type>,
        value_type,  
        typename zip_helper<iterator_tuple_type>::system,
        typename zip_helper<iterator_tuple_type>::iterator_category,
        reference,
        typename zip_helper<iterator_tuple_type>::difference_type
    > type;


}; // end zip_iterator_base
#endif



template<typename tuple_type>
struct getter_helper {};

#ifdef HAVE_THRUST
template <typename ... T>
struct getter_helper<thrust::tuple<T ...>> {
    typedef typename thrust::tuple<T...> tuple_type; 
    typedef typename thrust::tuple<T& ...> tuple_reference; 
    template <unsigned int N>
    using return_type = thrust::tuple_element<N,tuple_type>;
    typedef typename thrust::tuple_element<0,tuple_type>::type first_type; 
    typedef typename std::is_reference<first_type> is_reference;

    template<std::size_t... I>
    static tuple_reference make_reference(tuple_type& tuple, detail::index_sequence<I...>) {
        return thrust::tie(thrust::get<I>(tuple)...);
    }

};
#endif

template <typename ... T>
struct getter_helper<std::tuple<T ...>> {
    typedef typename std::tuple<T...> tuple_type; 
    typedef typename std::tuple<T& ...> tuple_reference; 
    template <unsigned int N>
    using return_type = std::tuple_element<N,tuple_type>;
    typedef typename std::tuple_element<0,tuple_type>::type first_type; 
    typedef typename std::is_reference<first_type> is_reference;

    template<std::size_t... I>
    static tuple_reference make_reference(tuple_type& tuple, detail::index_sequence<I...>) {
        return std::tie(std::get<I>(tuple)...);
    }

};

__aboria_hd_warning_disable__
template<typename reference, typename iterator_tuple_type, std::size_t... I>
CUDA_HOST_DEVICE
static reference make_reference(const iterator_tuple_type& tuple, index_sequence<I...>) {
    return reference(*(get_impl<I>(tuple))...);
}

__aboria_hd_warning_disable__
template<typename pointer, typename tuple_type, std::size_t... I>
CUDA_HOST_DEVICE
static pointer make_pointer(tuple_type&& tuple, index_sequence<I...>) {
    return pointer(&(get_impl<I>(std::forward(tuple)))...);
}

__aboria_hd_warning_disable__
template<typename iterator_tuple_type, std::size_t... I>
CUDA_HOST_DEVICE
static void increment_impl(iterator_tuple_type& tuple, index_sequence<I...>) {
    //using expander = int[];
    //(void)expander { 0, (++std::get<I>(tuple),0)...};
    int dummy[] = { 0, (++get_impl<I>(tuple),0)...};
    static_cast<void>(dummy);
}

__aboria_hd_warning_disable__
template<typename iterator_tuple_type, std::size_t... I>
CUDA_HOST_DEVICE
static void decrement_impl(iterator_tuple_type& tuple, index_sequence<I...>) {
    int dummy[] = { 0, (--get_impl<I>(tuple),0)...};
    static_cast<void>(dummy);
}

__aboria_hd_warning_disable__
template<typename iterator_tuple_type, std::size_t... I>
CUDA_HOST_DEVICE
static void advance_impl(iterator_tuple_type& tuple, 
        const typename zip_helper<iterator_tuple_type>::difference_type n,  index_sequence<I...>) {
    int dummy[] = { 0, (get_impl<I>(tuple)+=n,0)...};
    static_cast<void>(dummy);
}

template <typename ZipIterator, std::size_t... I>
typename ZipIterator::tuple_raw_pointer 
iterator_to_raw_pointer_impl(const ZipIterator& arg, index_sequence<I...>) {
#ifdef HAVE_THRUST
    return typename ZipIterator::tuple_raw_pointer(thrust::raw_pointer_cast(&*get_impl<I>(arg.get_tuple()))...);
#else
    return typename ZipIterator::tuple_raw_pointer(&*get_impl<I>(arg.get_tuple())...);
#endif
}
    
template <typename iterator_tuple_type, typename mpl_vector_type>
typename zip_iterator<iterator_tuple_type,mpl_vector_type>::tuple_raw_pointer
iterator_to_raw_pointer(const zip_iterator<iterator_tuple_type,mpl_vector_type>& arg, std::true_type) {
    typedef typename zip_helper<iterator_tuple_type>::index_type index_type;
    return iterator_to_raw_pointer_impl(arg,index_type());
}

template <typename Iterator>
typename std::iterator_traits<Iterator>::value_type*
iterator_to_raw_pointer(const Iterator& arg, std::false_type) {
#ifdef HAVE_THRUST
    return thrust::raw_pointer_cast(&*arg);
#else
    return &*arg;
#endif
}


template <typename T>
struct is_zip_iterator {
    typedef std::false_type type;
    static const bool value = false; 
};

template <typename tuple_type, typename mpl_vector_type>
struct is_zip_iterator<zip_iterator<tuple_type,mpl_vector_type>> {
    typedef std::true_type type;
    static const bool value = true; 
};



}
}



#endif //GET_DETAIL_H_
