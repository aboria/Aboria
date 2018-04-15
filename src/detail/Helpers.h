#ifndef HELPERS_DETAIL_H_ 
#define HELPERS_DETAIL_H_ 


#include <boost/iterator/iterator_facade.hpp>
#include "boost/mpl/contains.hpp"

#include "../CudaInclude.h"

#ifdef HAVE_THRUST
#include <thrust/iterator/iterator_traits.h>

namespace thrust {
    template <>
    struct iterator_traits<thrust::null_type> {
        typedef thrust::null_type value_type;
        typedef thrust::null_type reference;
        typedef thrust::null_type pointer;
    };
   

namespace detail {

    template <>
    struct pointer_traits<thrust::null_type> {
        typedef null_type raw_pointer;
    };

}
}
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




template<typename tuple_type>
struct getter_helper {};

#ifdef HAVE_THRUST

template<typename T>
struct is_thrust_device_reference : std::false_type { };

template<typename T>
struct is_thrust_device_reference<thrust::device_reference<T>> : std::true_type { };

template <typename ... T>
struct getter_helper<thrust::tuple<T ...>> {
    typedef typename thrust::tuple<T...> tuple_type; 

    typedef typename thrust::tuple<
                    typename std::conditional<
                        std::is_same<T,thrust::null_type>::value,
                        T,
                        T&>::type ...> tuple_reference; 

    typedef typename thrust::tuple<
                    typename std::conditional<
                        std::is_same<T,thrust::null_type>::value,
                        T,
                        thrust::device_reference<   
                            typename std::remove_reference<T>::type>
                            >::type ...> tuple_device_reference; 

    /* not needed anymore
    typedef typename thrust::tuple<
                    typename std::conditional<
                        std::is_same<T,thrust::null_type>::value,
                        T,
                        typename std::conditional<
                            is_thrust_device_reference<T>::value,
                            typename T::value_type,
                            typename std::decay<T>::type
                        >::type
                    >::type ...> tuple_value_type;
                    */
                        


    template <unsigned int N>
    using return_type = thrust::tuple_element<N,tuple_type>;
    typedef typename thrust::tuple_element<0,tuple_type>::type first_type; 
    typedef typename std::is_reference<first_type> is_reference;

    typedef make_index_sequence<thrust::tuple_size<tuple_type>::value> index_type;

    template<std::size_t... I>
    CUDA_HOST_DEVICE
    static tuple_reference make_reference(tuple_type& tuple, detail::index_sequence<I...>) {
        return thrust::tie(thrust::get<I>(tuple)...);
    }

    template<std::size_t... I>
    CUDA_HOST_DEVICE
    static tuple_reference raw_reference_cast(const tuple_device_reference& tuple, 
                                              detail::index_sequence<I...>) {
        return tuple_reference(thrust::raw_reference_cast(thrust::get<I>(tuple))...);
    }

};
#endif

template <typename ... T>
struct getter_helper<std::tuple<T ...>> {
    typedef typename std::tuple<T...> tuple_type; 
    typedef typename std::tuple<T& ...> tuple_reference; 
    typedef typename std::tuple<
                        typename std::decay<T>::type ...> tuple_value_type; 

    template <unsigned int N>
    using return_type = std::tuple_element<N,tuple_type>;
    typedef typename std::tuple_element<0,tuple_type>::type first_type; 
    typedef typename std::is_reference<first_type> is_reference;
    typedef make_index_sequence<std::tuple_size<tuple_type>::value> index_type;

    template<std::size_t... I>
    static tuple_reference make_reference(tuple_type& tuple, detail::index_sequence<I...>) {
        return std::tie(std::get<I>(tuple)...);
    }

};

template<typename T>
struct remove_pointer_or_reference_for_null_type {
    typedef T type;
};

#if HAVE_THRUST
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
struct zip_helper<std::tuple<T ...>> {
    //typedef std::false_type is_thrust;
    typedef std::tuple<T...> tuple_iterator_type; 
    typedef std::tuple<typename std::iterator_traits<T>::value_type ...> tuple_value_type; 
    typedef std::tuple<typename std::iterator_traits<T>::reference ...> tuple_reference; 
    typedef std::tuple<typename std::iterator_traits<T>::pointer...> tuple_pointer; 
    typedef tuple_pointer tuple_raw_pointer;
    typedef tuple_reference tuple_raw_reference;


    /*
typedef std::tuple<
        typename detail::remove_pointer_or_reference_for_null_type<
            typename std::iterator_traits<T>::value_type*>::type...
        > tuple_raw_pointer; 


    typedef std::tuple<
        typename detail::remove_pointer_or_reference_for_null_type<
            typename std::iterator_traits<T>::value_type&>::type...
        > tuple_raw_reference; 
    */
    typedef typename std::tuple<T...> iterator_tuple_type;

    template <unsigned int N>
    using tuple_element = std::tuple_element<N,iterator_tuple_type>;

    typedef typename std::iterator_traits<typename tuple_element<0>::type>::difference_type difference_type;
    typedef typename std::iterator_traits<typename tuple_element<0>::type>::iterator_category iterator_category;
    typedef make_index_sequence<std::tuple_size<iterator_tuple_type>::value> index_type;
    typedef boost::iterator_core_access iterator_core_access;

    template<std::size_t... I>
    static void increment_impl(tuple_iterator_type& tuple, index_sequence<I...>) {
        int dummy[] = { 0, (++std::get<I>(tuple),0)...};
        static_cast<void>(dummy);
    }

    template<std::size_t... I>
    static void decrement_impl(tuple_iterator_type& tuple, index_sequence<I...>) {
        int dummy[] = { 0, (--std::get<I>(tuple),0)...};
        static_cast<void>(dummy);
    }

    template<std::size_t... I>
    static void advance_impl(tuple_iterator_type& tuple, 
                difference_type n,  index_sequence<I...>) {
        int dummy[] = { 0, (std::get<I>(tuple)+=n,0)...};
        static_cast<void>(dummy);
    }

    template<std::size_t... I>
    static tuple_reference make_reference(const tuple_iterator_type& tuple, 
                                        index_sequence<I...>) {
        return tuple_reference(*(std::get<I>(tuple))...);
    }

    template <std::size_t... I>
    static tuple_raw_pointer make_raw_pointer(const tuple_iterator_type& arg, 
            index_sequence<I...>) {
        return tuple_raw_pointer(
                &*std::get<I>(arg)...
                );
    }

};


#ifdef HAVE_THRUST
template <typename T0,
          typename T1,
          typename T2,
          typename T3,
          typename T4,
          typename T5,
          typename T6,
          typename T7,
          typename T8,
          typename T9
          >
struct zip_helper<thrust::tuple<T0,T1,T2,T3,T4,T5,T6,T7,T8,T9>> {
    typedef thrust::tuple<T0,T1,T2,T3,T4,T5,T6,T7,T8,T9> tuple_iterator_type; 
    typedef thrust::tuple<
        typename thrust::iterator_traits<T0>::value_type, 
        typename thrust::iterator_traits<T1>::value_type, 
        typename thrust::iterator_traits<T2>::value_type, 
        typename thrust::iterator_traits<T3>::value_type, 
        typename thrust::iterator_traits<T4>::value_type, 
        typename thrust::iterator_traits<T5>::value_type, 
        typename thrust::iterator_traits<T6>::value_type, 
        typename thrust::iterator_traits<T7>::value_type, 
        typename thrust::iterator_traits<T8>::value_type, 
        typename thrust::iterator_traits<T9>::value_type 
        > tuple_value_type; 
    typedef thrust::tuple<
        typename thrust::iterator_traits<T0>::reference, 
        typename thrust::iterator_traits<T1>::reference, 
        typename thrust::iterator_traits<T2>::reference, 
        typename thrust::iterator_traits<T3>::reference, 
        typename thrust::iterator_traits<T4>::reference, 
        typename thrust::iterator_traits<T5>::reference, 
        typename thrust::iterator_traits<T6>::reference, 
        typename thrust::iterator_traits<T7>::reference, 
        typename thrust::iterator_traits<T8>::reference, 
        typename thrust::iterator_traits<T9>::reference 
        > tuple_reference; 
    typedef thrust::tuple<
        typename thrust::iterator_traits<T0>::pointer, 
        typename thrust::iterator_traits<T1>::pointer, 
        typename thrust::iterator_traits<T2>::pointer, 
        typename thrust::iterator_traits<T3>::pointer, 
        typename thrust::iterator_traits<T4>::pointer, 
        typename thrust::iterator_traits<T5>::pointer, 
        typename thrust::iterator_traits<T6>::pointer, 
        typename thrust::iterator_traits<T7>::pointer, 
        typename thrust::iterator_traits<T8>::pointer, 
        typename thrust::iterator_traits<T9>::pointer 
        > tuple_pointer; 

    typedef thrust::tuple<
        typename remove_pointer_or_reference_for_null_type< 
            typename thrust::iterator_traits<T0>::value_type*
            >::type,
        typename remove_pointer_or_reference_for_null_type< 
            typename thrust::iterator_traits<T1>::value_type*
            >::type,
        typename remove_pointer_or_reference_for_null_type< 
            typename thrust::iterator_traits<T2>::value_type*
            >::type,
        typename remove_pointer_or_reference_for_null_type< 
            typename thrust::iterator_traits<T3>::value_type*
            >::type,
        typename remove_pointer_or_reference_for_null_type< 
            typename thrust::iterator_traits<T4>::value_type*
            >::type,
        typename remove_pointer_or_reference_for_null_type< 
            typename thrust::iterator_traits<T5>::value_type*
            >::type,
        typename remove_pointer_or_reference_for_null_type< 
            typename thrust::iterator_traits<T6>::value_type*
            >::type,
        typename remove_pointer_or_reference_for_null_type< 
            typename thrust::iterator_traits<T7>::value_type*
            >::type,
        typename remove_pointer_or_reference_for_null_type< 
            typename thrust::iterator_traits<T8>::value_type*
            >::type,
        typename remove_pointer_or_reference_for_null_type< 
            typename thrust::iterator_traits<T9>::value_type*
            >::type
                > tuple_raw_pointer; 

    typedef thrust::tuple<
        typename remove_pointer_or_reference_for_null_type< 
            typename thrust::detail::raw_reference<
                typename thrust::iterator_traits<T0>::reference
                >::type
            >::type,
        typename remove_pointer_or_reference_for_null_type< 
            typename thrust::detail::raw_reference<
                typename thrust::iterator_traits<T1>::reference
                >::type
            >::type,
        typename remove_pointer_or_reference_for_null_type< 
            typename thrust::detail::raw_reference<
                typename thrust::iterator_traits<T2>::reference
                >::type
            >::type,
        typename remove_pointer_or_reference_for_null_type< 
            typename thrust::detail::raw_reference<
                typename thrust::iterator_traits<T3>::reference
                >::type
            >::type,
        typename remove_pointer_or_reference_for_null_type< 
            typename thrust::detail::raw_reference<
                typename thrust::iterator_traits<T4>::reference
                >::type
            >::type,
        typename remove_pointer_or_reference_for_null_type< 
            typename thrust::detail::raw_reference<
                typename thrust::iterator_traits<T5>::reference
                >::type
            >::type,
        typename remove_pointer_or_reference_for_null_type< 
            typename thrust::detail::raw_reference<
                typename thrust::iterator_traits<T6>::reference
                >::type
            >::type,
        typename remove_pointer_or_reference_for_null_type< 
            typename thrust::detail::raw_reference<
                typename thrust::iterator_traits<T7>::reference
                >::type
            >::type,
        typename remove_pointer_or_reference_for_null_type< 
            typename thrust::detail::raw_reference<
                typename thrust::iterator_traits<T8>::reference
                >::type
            >::type,
        typename remove_pointer_or_reference_for_null_type< 
            typename thrust::detail::raw_reference<
                typename thrust::iterator_traits<T9>::reference
                >::type
            >::type
        > tuple_raw_reference; 

    typedef tuple_iterator_type iterator_tuple_type;

    template <unsigned int N>
    using tuple_element = thrust::tuple_element<N,iterator_tuple_type>;

    typedef typename thrust::iterator_traits<typename tuple_element<0>::type>::difference_type difference_type;
    typedef typename thrust::iterator_traits<typename tuple_element<0>::type>::iterator_category iterator_category;
    typedef typename thrust::iterator_system<typename tuple_element<0>::type> system;
    typedef make_index_sequence<thrust::tuple_size<iterator_tuple_type>::value> index_type;
    typedef thrust::iterator_core_access iterator_core_access;

    template<std::size_t... I>
    CUDA_HOST_DEVICE
    static void increment_impl(tuple_iterator_type& tuple, index_sequence<I...>) {
        int dummy[] = { 0, (++thrust::get<I>(tuple),0)...};
        static_cast<void>(dummy);
    }

    template<std::size_t... I>
    CUDA_HOST_DEVICE
    static void decrement_impl(tuple_iterator_type& tuple, index_sequence<I...>) {
        int dummy[] = { 0, (--thrust::get<I>(tuple),0)...};
        static_cast<void>(dummy);
    }

    template<std::size_t... I>
    CUDA_HOST_DEVICE
    static void advance_impl(tuple_iterator_type& tuple, 
                difference_type n,  index_sequence<I...>) {
        int dummy[] = { 0, (thrust::get<I>(tuple)+=n,0)...};
        static_cast<void>(dummy);
    }

    template<std::size_t... I>
    CUDA_HOST_DEVICE
    static tuple_reference make_reference(const tuple_iterator_type& tuple, 
                                        index_sequence<I...>) {
        return tuple_reference(*(thrust::get<I>(tuple))...);
    }

    template <std::size_t... I>
    CUDA_HOST_DEVICE
    static tuple_raw_pointer make_raw_pointer(const tuple_iterator_type& arg, 
            index_sequence<I...>) {
        return tuple_raw_pointer(
                thrust::raw_pointer_cast(&*thrust::get<I>(arg))...
                );
    }

};
#endif


}
}
#endif
