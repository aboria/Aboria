
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


template<typename tuple_of_iterators>
struct zip_helper {};


template <typename ... T>
struct zip_helper<tuple_ns::tuple<T ...>> {
    //typedef std::false_type is_thrust;
    typedef tuple_ns::tuple<T...> tuple_iterator_type; 
    typedef tuple_ns::tuple<typename tuple_ns::iterator_traits<T>::value_type ...> tuple_value_type; 
    typedef tuple_ns::tuple<typename tuple_ns::iterator_traits<T>::reference ...> tuple_reference; 
    typedef tuple_ns::tuple<typename tuple_ns::iterator_traits<T>::pointer...> tuple_pointer; 

    typedef tuple_ns::tuple<
        typename detail::remove_pointer_or_reference_for_null_type<
            typename tuple_ns::iterator_traits<T>::value_type*>::type...
        > tuple_raw_pointer; 
    typedef tuple_ns::tuple<
        typename detail::remove_pointer_or_reference_for_null_type<
            typename tuple_ns::iterator_traits<T>::value_type&>::type...
        > tuple_raw_reference; 
    typedef typename tuple_ns::tuple<T...> iterator_tuple_type;
    template <unsigned int N>
    using tuple_element = tuple_ns::tuple_element<N,iterator_tuple_type>;
    typedef typename tuple_ns::iterator_traits<typename tuple_element<0>::type>::difference_type difference_type;
    typedef typename tuple_ns::iterator_traits<typename tuple_element<0>::type>::iterator_category iterator_category;
#if defined(__CUDACC__)
    typedef typename thrust::iterator_system<typename tuple_element<0>::type> system;
#endif
    typedef make_index_sequence<tuple_ns::tuple_size<iterator_tuple_type>::value> index_type;
};

template<typename IteratorTuple, typename MplVector>
struct zip_iterator_base {

    typedef getter_type<typename zip_helper<IteratorTuple>::tuple_value_type,MplVector> value_type;
    typedef getter_type<typename zip_helper<IteratorTuple>::tuple_reference,MplVector> reference;
  
 public:
  
#if defined(__CUDACC__)
    typedef iterator_facade_ns::iterator_facade<
        zip_iterator<IteratorTuple,MplVector>,
        value_type,  
        typename zip_helper<IteratorTuple>::system,
        typename zip_helper<IteratorTuple>::iterator_category,
        reference,
        typename zip_helper<IteratorTuple>::difference_type
    > type;
#else
    typedef iterator_facade_ns::iterator_facade<
        zip_iterator<IteratorTuple,MplVector>,
        value_type,  
        typename zip_helper<IteratorTuple>::iterator_category,
        reference,
        typename zip_helper<IteratorTuple>::difference_type
    > type;
#endif

}; // end zip_iterator_base


template<typename tuple_type>
struct getter_helper {};

template <typename ... T>
struct getter_helper<tuple_ns::tuple<T ...>> {
    typedef typename tuple_ns::tuple<T...> tuple_type; 
    typedef typename tuple_ns::tuple<T& ...> tuple_reference; 
    template <unsigned int N>
    using return_type = tuple_ns::tuple_element<N,tuple_type>;
    typedef typename tuple_ns::tuple_element<0,tuple_type>::type first_type; 
    typedef typename std::is_reference<first_type> is_reference;

    template<std::size_t... I>
    static tuple_reference make_reference(tuple_type& tuple, detail::index_sequence<I...>) {
        return tuple_ns::tie(tuple_ns::get<I>(tuple)...);
    }

};

__aboria_hd_warning_disable__
template<typename reference, typename iterator_tuple_type, std::size_t... I>
CUDA_HOST_DEVICE
static reference make_reference(const iterator_tuple_type& tuple, index_sequence<I...>) {
    return reference(*(tuple_ns::get<I>(tuple))...);
}

__aboria_hd_warning_disable__
template<typename pointer, typename tuple_type, std::size_t... I>
CUDA_HOST_DEVICE
static pointer make_pointer(tuple_type&& tuple, index_sequence<I...>) {
    return pointer(&(tuple_ns::get<I>(std::forward(tuple)))...);
}

__aboria_hd_warning_disable__
template<typename iterator_tuple_type, std::size_t... I>
CUDA_HOST_DEVICE
static void increment_impl(iterator_tuple_type& tuple, index_sequence<I...>) {
    //using expander = int[];
    //(void)expander { 0, (++std::get<I>(tuple),0)...};
    int dummy[] = { 0, (++tuple_ns::get<I>(tuple),0)...};
    static_cast<void>(dummy);
}

__aboria_hd_warning_disable__
template<typename iterator_tuple_type, std::size_t... I>
CUDA_HOST_DEVICE
static void decrement_impl(iterator_tuple_type& tuple, index_sequence<I...>) {
    int dummy[] = { 0, (--tuple_ns::get<I>(tuple),0)...};
    static_cast<void>(dummy);
}

__aboria_hd_warning_disable__
template<typename iterator_tuple_type, std::size_t... I>
CUDA_HOST_DEVICE
static void advance_impl(iterator_tuple_type& tuple, 
        const typename zip_helper<iterator_tuple_type>::difference_type n,  index_sequence<I...>) {
    int dummy[] = { 0, (tuple_ns::get<I>(tuple)+=n,0)...};
    static_cast<void>(dummy);
}

template <typename ZipIterator, std::size_t... I>
typename ZipIterator::tuple_raw_pointer 
iterator_to_raw_pointer_impl(const ZipIterator& arg, index_sequence<I...>) {
#if defined(__CUDACC__)
    return typename ZipIterator::tuple_raw_pointer(thrust::raw_pointer_cast(&*thrust::get<I>(arg.get_tuple()))...);
#else
    return typename ZipIterator::tuple_raw_pointer(&*std::get<I>(arg.get_tuple())...);
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
#if defined(__CUDACC__)
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
