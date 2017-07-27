
#ifndef GET_H_ 
#define GET_H_ 

#include <boost/mpl/vector.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/range_c.hpp>
#include <algorithm>
#include <tuple>
#include <type_traits>

#include "detail/GetterType.h"

namespace Aboria {

template <typename Tuple, typename MplVector> 
using getter_type = detail::getter_type_base<
                detail::pointer_check<Tuple>::value, MplVector, Tuple>;

template <typename iterator_tuple_type, typename mpl_vector_type>
class zip_iterator;

}
 
#include "detail/Get.h"
#include "Log.h"


namespace Aboria {

/*
template <typename tuple_type, typename mpl_vector_type> 
void
swap(getter_type<tuple_type,mpl_vector_type> x,
     getter_type<tuple_type,mpl_vector_type> y) {
    x.swap(y);
}
*/

template <typename tuple_type, typename tuple_type2, typename mpl_vector_type> 
typename std::enable_if<
    detail::getter_helper<tuple_type>::is_reference::value
    && !detail::getter_helper<tuple_type2>::is_reference::value
    >::type
swap(getter_type<tuple_type,mpl_vector_type> x,
     getter_type<tuple_type2,mpl_vector_type>& y) {
    x.swap_via_tie(y.data,
                detail::make_index_sequence<tuple_ns::tuple_size<tuple_type2>::value>());
}

template <typename tuple_type, typename tuple_type2, typename mpl_vector_type> 
typename std::enable_if<
    !detail::getter_helper<tuple_type>::is_reference::value
    && detail::getter_helper<tuple_type2>::is_reference::value
    >::type
swap(getter_type<tuple_type,mpl_vector_type>& x,
     getter_type<tuple_type2,mpl_vector_type> y) {
    y.swap_via_tie(x.data,
                detail::make_index_sequence<tuple_ns::tuple_size<tuple_type2>::value>());

}

template <typename tuple_type, typename mpl_vector_type> 
void swap(getter_type<tuple_type,mpl_vector_type> x,
          getter_type<tuple_type,mpl_vector_type> y) {
    y.swap(x);
}

/*
template <typename tuple_type, typename tuple_type2, typename mpl_vector_type> 
typename std::enable_if<
    detail::getter_helper<tuple_type>::is_reference::value
    && detail::getter_helper<tuple_type2>::is_reference::value
    >::type
swap(getter_type<tuple_type,mpl_vector_type> x,
     getter_type<tuple_type2,mpl_vector_type> y) {
    y.swap(x);
}

template <typename tuple_type, typename tuple_type2, typename mpl_vector_type> 
typename std::enable_if<
    !detail::getter_helper<tuple_type>::is_reference::value
    && !detail::getter_helper<tuple_type2>::is_reference::value
    >::type
swap(getter_type<tuple_type,mpl_vector_type>& x,
     getter_type<tuple_type2,mpl_vector_type>& y) {
    y.swap(x);
}
*/


template <typename iterator_tuple_type, typename mpl_vector_type>
class zip_iterator: 
    public detail::zip_iterator_base<iterator_tuple_type,mpl_vector_type>::type {
public:
    typedef iterator_tuple_type tuple_type;
    typedef typename detail::zip_iterator_base<tuple_type,mpl_vector_type>::value_type value_type;
    typedef typename detail::zip_iterator_base<tuple_type,mpl_vector_type>::reference reference;

    typedef typename detail::zip_helper<iterator_tuple_type>::difference_type difference_type;
    typedef typename detail::zip_helper<iterator_tuple_type>::iterator_category iterator_category;

    typedef getter_type<typename detail::zip_helper<iterator_tuple_type>::tuple_pointer,mpl_vector_type> pointer;

    // Note: Can't call it raw_pointer or thrust thinks it is a trivial iterator!
    typedef getter_type<typename detail::zip_helper<iterator_tuple_type>::tuple_raw_pointer,mpl_vector_type> tuple_raw_pointer;
    typedef getter_type<typename detail::zip_helper<iterator_tuple_type>::tuple_raw_reference,mpl_vector_type> tuple_raw_reference;

    template <typename T>
    using elem_by_type = detail::get_elem_by_type<T,mpl_vector_type>;

    template<typename T>
    struct return_type {
        static const size_t N = elem_by_type<T>::index;
        typedef const typename detail::zip_helper<iterator_tuple_type>::template tuple_element<N>::type type;
    };

    zip_iterator() {}

    explicit zip_iterator(iterator_tuple_type iter) : iter(iter) {}

    template <typename ...T>
    explicit zip_iterator(T... args) : iter(args...) {}

    CUDA_HOST_DEVICE
    const iterator_tuple_type & get_tuple() const { 
        #if defined(__CUDA_ARCH__)
        ERROR_CUDA("Cannot use `zip_iterator` in device code");
        #endif
        return iter; 
    }

    CUDA_HOST_DEVICE
    iterator_tuple_type & get_tuple() { 
        #if defined(__CUDA_ARCH__)
        ERROR_CUDA("Cannot use `zip_iterator` in device code");
        #endif
        return iter; 
    }

private:

    typedef typename detail::zip_helper<iterator_tuple_type>::index_type index_type;

    void increment() { 
        detail::zip_helper<iterator_tuple_type>::increment_impl(iter,index_type()); 
    }
    
    void decrement() { 
        detail::zip_helper<iterator_tuple_type>::decrement_impl(iter,index_type()); 
    }

    bool equal(zip_iterator const& other) const { 
        return detail::get_impl<0>(other.iter) 
            == detail::get_impl<0>(iter);
    }

    reference dereference() const { 
        return reference(
                detail::zip_helper<iterator_tuple_type>::make_reference(
                    iter,index_type())
                ); 
    }

    difference_type distance_to(zip_iterator const& other) const { 
        return detail::get_impl<0>(other.iter) 
             - detail::get_impl<0>(iter);
    }

    void advance(difference_type n) { 
        detail::zip_helper<iterator_tuple_type>::advance_impl(iter,n,index_type()); 
    }

    iterator_tuple_type iter;
    friend class boost::iterator_core_access;
};



template <typename Iterator>
auto iterator_to_raw_pointer(const Iterator& arg) ->
    decltype(detail::iterator_to_raw_pointer(arg,
                typename detail::is_zip_iterator<Iterator>::type()
                           ))

{
    return detail::iterator_to_raw_pointer(arg,
            typename detail::is_zip_iterator<Iterator>::type()
                           );
}


//
// Particle getters/setters
//
//


/// get a variable from a particle \p arg
/// \param T Variable type
/// \param arg the particle
/// \return a const reference to a T::value_type holding the variable data
///

/*
template<typename T, typename ValueType>
CUDA_HOST_DEVICE
auto get(const ValueType& arg) ->
    decltype(detail::get_impl<ValueType::template elem_by_type<T>::index>(arg.get_tuple()))
{
    //std::cout << "get const reference" << std::endl;
    return detail::get_impl<ValueType::template elem_by_type<T>::index>(arg.get_tuple());
    //return arg.template get<T>();
}
*/

template<typename T, typename ValueType>
CUDA_HOST_DEVICE
typename ValueType::template return_type<T>::type const & 
get(const ValueType& arg) {
    //std::cout << "get reference" << std::endl;
    return detail::get_impl<ValueType::template elem_by_type<T>::index>(arg.get_tuple());
    //return arg.template get<T>();
}

template<typename T, typename ValueType>
CUDA_HOST_DEVICE
typename ValueType::template return_type<T>::type & 
get(ValueType& arg) {
    //std::cout << "get reference" << std::endl;
    return detail::get_impl<ValueType::template elem_by_type<T>::index>(arg.get_tuple());
    //return arg.template get<T>();
}

template<typename T, typename ValueType>
CUDA_HOST_DEVICE
typename ValueType::template return_type<T>::type & 
get(ValueType&& arg) {
    //std::cout << "get reference" << std::endl;
    return detail::get_impl<ValueType::template elem_by_type<T>::index>(arg.get_tuple());
    //return arg.template get<T>();
}

template<unsigned int N, typename ValueType>
CUDA_HOST_DEVICE
const typename detail::getter_helper<typename ValueType::tuple_type>::template return_type<N>::type &
get_by_index(const ValueType& arg) {
    return detail::get_impl<N>(arg.get_tuple());
}

template<unsigned int N, typename ValueType>
CUDA_HOST_DEVICE
typename detail::getter_helper<typename ValueType::tuple_type>::template return_type<N>::type &
get_by_index(ValueType& arg) {
    return detail::get_impl<N>(arg.get_tuple());
}

template<unsigned int N, typename ValueType>
CUDA_HOST_DEVICE
typename detail::getter_helper<typename ValueType::tuple_type>::template return_type<N>::type &
get_by_index(ValueType&& arg) {
    return tuple_ns::get<N>(arg.get_tuple());        
}

}




#endif //GET_H_
