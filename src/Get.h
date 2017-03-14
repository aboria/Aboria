
#ifndef GET_H_ 
#define GET_H_ 

#include <boost/mpl/vector.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/range_c.hpp>
#include <algorithm>
#include <tuple>
#include <type_traits>


namespace Aboria {

template <typename TUPLE, typename mpl_vector_type> 
struct getter_type;
template <typename iterator_tuple_type, typename mpl_vector_type>
class zip_iterator;

}
 
#include "detail/Get.h"


namespace Aboria {

template <typename TUPLE, typename mpl_vector_type> 
struct getter_type{
    typedef mpl_vector_type mpl_type_vector;
    typedef TUPLE tuple_type;
    typedef typename detail::getter_helper<tuple_type>::tuple_reference tuple_reference;
    template <typename T>
    using elem_by_type = detail::get_elem_by_type<T,mpl_vector_type>;
    template <typename T>
    using return_type = typename detail::getter_helper<tuple_type>::template return_type<elem_by_type<T>::index>;

    //typedef typename detail::zip_helper<tuple_type>::pointer pointer;

    CUDA_HOST_DEVICE
    getter_type() {}
    CUDA_HOST_DEVICE
    explicit getter_type(const tuple_type& data):data(data) {}
    CUDA_HOST_DEVICE
    getter_type(const getter_type& other):data(other.data) {}
    CUDA_HOST_DEVICE
    getter_type(getter_type&& other):data(std::move(other.data)) {}

    template <typename T=tuple_reference, typename = typename
    std::enable_if<
        !std::is_same<T,tuple_type>::value
        >::type>
    CUDA_HOST_DEVICE
    getter_type(const getter_type<tuple_reference,mpl_vector_type>& other):
        data(other.data) {}

    template <typename T=tuple_reference, typename = typename
    std::enable_if<
        !std::is_same<T,tuple_type>::value
        >::type>
    CUDA_HOST_DEVICE
    getter_type(getter_type<tuple_reference,mpl_vector_type>&& other):
        data(std::move(other.data)) {}

    template <typename T, typename = typename
    std::enable_if<
        (!std::is_same<T,tuple_type>::value)
        && (!std::is_same<T,tuple_reference>::value)
        >::type>
    CUDA_HOST_DEVICE
    getter_type(const getter_type<T,mpl_vector_type>& other):
        data(other.data) {}

    template <typename T, typename = typename
    std::enable_if<
        (!std::is_same<T,tuple_type>::value)
        && (!std::is_same<T,tuple_reference>::value)
        >::type>
    CUDA_HOST_DEVICE
    getter_type(getter_type<T,mpl_vector_type>&& other):
        data(std::move(other.data)) {}

    
    template <typename T1, typename T2, typename... T3>
    CUDA_HOST_DEVICE
    getter_type(T1&& arg1, T2&& arg2, T3&&... args):data(std::forward<T1>(arg1),std::forward<T2>(arg2),std::forward<T3>(args)...) {}

    template <typename T>
    int throw_away(const T& in) {
        return 0;
    }

    template<typename Tuple, std::size_t... I>
    void copy_impl(const Tuple& other_data, detail::index_sequence<I...>) {
        int dummy[] = { 0, throw_away(tuple_ns::get<I>(data) = tuple_ns::get<I>(other_data))... };
        static_cast<void>(dummy);
    }

    CUDA_HOST_DEVICE
    getter_type& operator=( const getter_type& other ) {
        //copy_impl(other.data,detail::make_index_sequence<tuple_ns::tuple_size<tuple_type>::value>());
        data = other.data;
        return *this;
    }
    CUDA_HOST_DEVICE
    getter_type& operator=( getter_type&& other ) {
        //copy_impl(other.data,detail::make_index_sequence<tuple_ns::tuple_size<tuple_type>::value>());
        data = other.data;
        return *this;
    }

    template <typename T=tuple_reference, typename = typename
    std::enable_if<
        !std::is_same<T,tuple_type>::value
        >::type>
    CUDA_HOST_DEVICE
    getter_type& operator=( const getter_type<tuple_reference,mpl_vector_type>& other) {
        //copy_impl(other.data,detail::make_index_sequence<tuple_ns::tuple_size<tuple_type>::value>());
        data = other.data;
        return *this;
    }

    template <typename T=tuple_reference, typename = typename
    std::enable_if<
        !std::is_same<T,tuple_type>::value
        >::type>
    CUDA_HOST_DEVICE
    getter_type& operator=(getter_type<tuple_reference,mpl_vector_type>&& other) {
        //copy_impl(other.data,detail::make_index_sequence<tuple_ns::tuple_size<tuple_type>::value>());
        data = other.data;
        return *this;
    }

    template <typename T, typename = typename
    std::enable_if<
        (!std::is_same<T,tuple_reference>::value)
        && (!std::is_same<T,tuple_type>::value)
        >::type>
    CUDA_HOST_DEVICE
    getter_type& operator=( const getter_type<T,mpl_vector_type>& other) {
        //copy_impl(other.data,detail::make_index_sequence<tuple_ns::tuple_size<tuple_type>::value>());
        data = other.data;
        return *this;
    }

    template <typename T, typename = typename
    std::enable_if<
        (!std::is_same<T,tuple_reference>::value)
        && (!std::is_same<T,tuple_type>::value)
        >::type>
    CUDA_HOST_DEVICE
    getter_type& operator=( getter_type<T,mpl_vector_type>&& other) {
        //copy_impl(other.data,detail::make_index_sequence<tuple_ns::tuple_size<tuple_type>::value>());
        data = other.data;
        return *this;
    }

#if defined(__aboria_use_thrust_algorithms__) || defined(__CUDACC__)
    template <typename T1, typename T2, typename PointerType, typename DerivedType> 
    CUDA_HOST_DEVICE
    getter_type& operator=( const thrust::reference<getter_type<T1,T2>,PointerType,DerivedType>& other) {
        data = static_cast<getter_type<T1,T2>>(other).data;
        return *this;
    }
#endif
    
    CUDA_HOST_DEVICE
    void swap(getter_type &other) {
        data.swap(other.data);
    }

    template <typename tuple_type2,std::size_t... I>
    CUDA_HOST_DEVICE
    void swap_via_tie(tuple_type2 &tuple, detail::index_sequence<I...>) {
        tuple_type tmp = tuple_ns::tie(tuple_ns::get<I>(tuple)...);
        data.swap(tmp);
    }

    CUDA_HOST_DEVICE
    const tuple_type & get_tuple() const { return data; }
    CUDA_HOST_DEVICE
    tuple_type & get_tuple() { return data; }


    tuple_type data;
};

/*
 * specialisation for tuple of pointers. Note hack to detect tuple of pointers.. could be better?
 * TODO: api is not good, consider doing an iterator_facade type thing
 */
template <typename mpl_vector_type, typename FirstType, typename ... OtherTypes> 
struct getter_type<tuple_ns::tuple<FirstType*, OtherTypes...>, mpl_vector_type>{
    typedef tuple_ns::tuple<FirstType*, OtherTypes...> tuple_type;
    typedef mpl_vector_type mpl_type_vector;
    template <typename T>
    using elem_by_type = detail::get_elem_by_type<T,mpl_vector_type>;
    template <typename T>
    using return_type = typename detail::getter_helper<tuple_type>::template return_type<elem_by_type<T>::index>;

    typedef getter_type<typename detail::zip_helper<tuple_type>::tuple_reference,mpl_vector_type> reference;
    typedef typename detail::zip_helper<tuple_type>::difference_type difference_type;
    typedef typename detail::zip_helper<tuple_type>::index_type index_type;


    CUDA_HOST_DEVICE
    getter_type() {}
    CUDA_HOST_DEVICE
    explicit getter_type(const tuple_type& data):data(data) {}
    CUDA_HOST_DEVICE
    getter_type(const getter_type& other):data(other.data) {}
    CUDA_HOST_DEVICE
    getter_type(getter_type&& other):data(std::move(other.data)) {}

    template <typename tuple_type2, typename = typename
    std::enable_if<std::is_convertible<tuple_type2,tuple_type>::value>::type>
    CUDA_HOST_DEVICE
    getter_type(const getter_type<tuple_type2,mpl_vector_type>& other):data(other.data) {}

    template <typename tuple_type2, typename = typename
    std::enable_if<std::is_convertible<tuple_type2,tuple_type>::value>::type>
    CUDA_HOST_DEVICE
    getter_type(const getter_type<tuple_type2,mpl_vector_type>&& other):data(std::move(other.data)) {}

    template <typename T1, typename T2, typename... T3>
    CUDA_HOST_DEVICE
    getter_type(T1&& arg1, T2&& arg2, T3&&... args):data(std::forward<T1>(arg1),std::forward<T2>(arg2),std::forward<T3>(args)...) {}

    CUDA_HOST_DEVICE
    getter_type& operator=( const getter_type& other ) {
        data = other.data;
        return *this;
    }
    CUDA_HOST_DEVICE
    getter_type& operator=( getter_type&& other ) {
        data = std::move(other.data);
        return *this;
    }
    template <typename T1, typename T2> 
    CUDA_HOST_DEVICE
    getter_type& operator=( const getter_type<T1,T2>& other) {
        data = other.data;
        return *this;
    }
    template <typename T1, typename T2> 
    CUDA_HOST_DEVICE
    getter_type& operator=( getter_type<T1,T2>&& other) {
        data = std::move(other.data);
        return *this;
    }
    
    CUDA_HOST_DEVICE
    void swap(getter_type &other) {
        data.swap(other.data);
    }

    CUDA_HOST_DEVICE
    const tuple_type & get_tuple() const { return data; }
    CUDA_HOST_DEVICE
    tuple_type & get_tuple() { return data; }

    
    CUDA_HOST_DEVICE
    bool operator==(const getter_type& other) const {
        return equal(other);
    }

    CUDA_HOST_DEVICE
    getter_type& operator++() {
        advance(1);
        return *this;
    }
    
    CUDA_HOST_DEVICE
    getter_type operator++(int) {
        getter_type temp(*this);
        advance(1);
        return temp; // return saved state
    }


    CUDA_HOST_DEVICE
    getter_type operator+(const difference_type& n) const {
        getter_type ret(*this);
        ret.advance(n);
        return ret;
    }

    CUDA_HOST_DEVICE
    difference_type operator-(const getter_type& other) const {
        return other.distance_to(*this);
    }

    CUDA_HOST_DEVICE
    getter_type& operator--() {
        decrement();
        return *this;
    }
    
    CUDA_HOST_DEVICE
    getter_type operator--(int) {
        getter_type temp(*this);
        decrement();
        return temp; // return saved state
    }


    CUDA_HOST_DEVICE
    reference operator*() const {
        return dereference();
    }

    CUDA_HOST_DEVICE
    void increment() { 
        detail::increment_impl(data,index_type()); 
    }
    CUDA_HOST_DEVICE
    void decrement() { detail::decrement_impl(data,index_type()); }

    CUDA_HOST_DEVICE
    bool equal(getter_type const& other) const { 
        return tuple_ns::get<0>(other.data) == tuple_ns::get<0>(data);
    }

    CUDA_HOST_DEVICE
    reference dereference() const { return detail::make_reference<reference>(data,index_type()); }

    CUDA_HOST_DEVICE
    difference_type distance_to(getter_type const& other) const { 
        return tuple_ns::get<0>(other.data) - tuple_ns::get<0>(data);
    }

    CUDA_HOST_DEVICE
    void advance(difference_type n) { detail::advance_impl(data,n,index_type()); }

    tuple_type data;
};

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

    CUDA_HOST_DEVICE
    zip_iterator() {}

    CUDA_HOST_DEVICE
    explicit zip_iterator(iterator_tuple_type iter) : iter(iter) {}

    template <typename ...T>
    CUDA_HOST_DEVICE
    explicit zip_iterator(T... args) : iter(args...) {}

    CUDA_HOST_DEVICE
    const iterator_tuple_type & get_tuple() const { return iter; }

    CUDA_HOST_DEVICE
    iterator_tuple_type & get_tuple() { return iter; }

private:

    typedef typename detail::zip_helper<iterator_tuple_type>::index_type index_type;

    CUDA_HOST_DEVICE
    void increment() { 
        detail::increment_impl(iter,index_type()); 
    }
    
    CUDA_HOST_DEVICE
    void decrement() { detail::decrement_impl(iter,index_type()); }

    __aboria_hd_warning_disable__
    CUDA_HOST_DEVICE
    bool equal(zip_iterator const& other) const { 
        return tuple_ns::get<0>(other.iter) == tuple_ns::get<0>(iter);
    }

    CUDA_HOST_DEVICE
    reference dereference() const { return detail::make_reference<reference>(iter,index_type()); }

    __aboria_hd_warning_disable__
    CUDA_HOST_DEVICE
    difference_type distance_to(zip_iterator const& other) const { 
        return tuple_ns::get<0>(other.iter)-tuple_ns::get<0>(iter);
    }

    CUDA_HOST_DEVICE
    void advance(difference_type n) { detail::advance_impl(iter,n,index_type()); }

    iterator_tuple_type iter;
    friend class iterator_facade_ns::iterator_core_access;
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

/// get a variable from a particle \p arg
/// \param T Variable type
/// \param arg the particle
/// \return a const reference to a T::value_type holding the variable data
///
template<typename T, typename value_type>
CUDA_HOST_DEVICE
typename value_type::template return_type<T>::type const & 
get(const value_type& arg) {
    //std::cout << "get const reference" << std::endl;
    return tuple_ns::get<value_type::template elem_by_type<T>::index>(arg.get_tuple());
    //return arg.template get<T>();
}

template<typename T, typename value_type>
CUDA_HOST_DEVICE
typename value_type::template return_type<T>::type & 
get(value_type& arg) {
    //std::cout << "get reference" << std::endl;
    return tuple_ns::get<value_type::template elem_by_type<T>::index>(arg.get_tuple());
    //return arg.template get<T>();
}

template<typename T, typename value_type>
CUDA_HOST_DEVICE
typename value_type::template return_type<T>::type & 
get(value_type&& arg) {
    //std::cout << "get reference" << std::endl;
    return tuple_ns::get<value_type::template elem_by_type<T>::index>(arg.get_tuple());
    //return arg.template get<T>();
}

template<unsigned int N, typename value_type>
CUDA_HOST_DEVICE
const typename detail::getter_helper<typename value_type::tuple_type>::template return_type<N>::type &
get_by_index(const value_type& arg) {
    return tuple_ns::get<N>(arg.get_tuple());
}

template<unsigned int N, typename value_type>
CUDA_HOST_DEVICE
typename detail::getter_helper<typename value_type::tuple_type>::template return_type<N>::type &
get_by_index(value_type& arg) {
    return tuple_ns::get<N>(arg.get_tuple());        
}

template<unsigned int N, typename value_type>
CUDA_HOST_DEVICE
typename detail::getter_helper<typename value_type::tuple_type>::template return_type<N>::type &
get_by_index(value_type&& arg) {
    return tuple_ns::get<N>(arg.get_tuple());
}

}




#endif //GET_H_
