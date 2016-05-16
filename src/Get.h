
#ifndef GET_H_ 
#define GET_H_ 

#include <boost/mpl/vector.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <tuple>


namespace Aboria {

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
struct elem_by_type {
    typedef T type;
    typedef typename T::value_type value_type;

    /// 
    /// iter is a boost mpl iterator to the found element in Variable T
    typedef typename mpl::find<mpl_type_vector,T>::type iter;
    BOOST_MPL_ASSERT_NOT(( boost::is_same< typename mpl::end<mpl_type_vector>::type, typename iter::type > ));

    /// 
    /// index contains the index of the found element
    static const size_t index = iter::pos::value;
};


/// helper class to find an element of mpl_type_vector from an
/// unsigned int index I
template<unsigned int I, typename mpl_type_vector>
struct elem_by_index {
    BOOST_MPL_ASSERT_RELATION( (mpl::size<mpl_type_vector>::type::value), >, I );
    typedef typename mpl::at<mpl_type_vector,mpl::int_<I> > type;

    /// 
    /// value_type is the variable's value_type at index I 
    typedef typename type::value_type value_type;
    static const size_t index = I;
};

template<size_t N, typename tuple_type>
struct element {
    typedef typename boost::tuples::element<N,tuple_type>::type type;
};

template<size_t N, typename... tuple_args>
struct element<N,std::tuple<tuple_args...>> {
    typedef typename std::tuple_element<N,std::tuple<tuple_args...>>::type type;
};

/*
template<typename ARG>
struct unpack_tuple_types {};

template <typename ... TYPES>
struct unpack_tuple_types< std::tuple<TYPES...> > {
    typedef mpl::vector<typename std::remove_reference<TYPES>::type...> mpl_type_vector;
};

template <typename ... TYPES>
struct unpack_tuple_types< boost::tuple<TYPES...> > {
    typedef mpl::vector<typename std::remove_reference<TYPES>::type...> mpl_type_vector;
};
*/

template <typename tuple_type, typename mpl_vector_type> 
struct getter_type {
    typedef mpl_vector_type mpl_type_vector;
    template <typename T>
    using elem_by_type = elem_by_type<T,mpl_type_vector>;
    template <typename T>
    using return_type = element<elem_by_type<T>::index,tuple_type>;

    getter_type() {}
    getter_type(const tuple_type& data):data(data) {}

    template<typename T>
    typename return_type<T>::type get() const {
        return std::get<elem_by_type<T>::index>(data);        
    }

    const tuple_type & get_tuple() const { return data; }
    tuple_type & get_tuple() { return data; }
            
    tuple_type data;
};


template <typename iterator_tuple_type, typename value_type, typename reference_type, typename mpl_vector_type>
class zip_iterator: public boost::iterator_facade<zip_iterator<iterator_tuple_type,value_type,reference_type,mpl_vector_type>,
    value_type,
    boost::random_access_traversal_tag,
    reference_type> {

    typedef mpl_vector_type mpl_type_vector;
    template <typename T>
    using elem_by_type = elem_by_type<T,mpl_type_vector>;
    template <typename T>
    using return_type = element<elem_by_type<T>::index,iterator_tuple_type>;

    typedef typename element<0,iterator_tuple_type>::type first_type;
public:
    explicit zip_iterator(iterator_tuple_type iter) : iter(iter) {}
    typedef typename std::iterator_traits<first_type>::difference_type difference_type;
private:
    typedef make_index_sequence<std::tuple_size<iterator_tuple_type>::value> index_type;

    template<std::size_t... I>
    static reference_type make_reference(iterator_tuple_type& tuple, index_sequence<I...>) {
        return reference_type(std::get<I>(tuple)...);
    }

    template<std::size_t... I>
    static void increament_impl(iterator_tuple_type& tuple, index_sequence<I...>) {
        using expander = int[];
        (void)expander { 0, (++std::get<I>(tuple))...};
    }

    template<std::size_t... I>
    static void decrement_impl(iterator_tuple_type& tuple, index_sequence<I...>) {
        using expander = int[];
        (void)expander { 0, (--std::get<I>(tuple))...};
    }

    template<std::size_t... I>
    static void advance_impl(iterator_tuple_type& tuple, const difference_type n,  index_sequence<I...>) {
        using expander = int[];
        (void)expander { 0, (std::get<I>(tuple) += n)...};
    }

    void increment() { increament_impl(iter,index_type()); }
    void decrement() { decrement_impl(iter,index_type()); }

    bool equal(zip_iterator const& other) const { return std::get<0>(iter) == std::get<0>(other.iter); }
    reference_type dereference() const { return make_reference(iter,index_type()); }
    difference_type distance_to(zip_iterator const& other) const { return std::get<0>(other.iter) - std::get<0>(iter); }
    void advance(difference_type n) { advance_impl(iter,index_type()); }

    template<typename T>
    typename return_type<T>::type & get() const {
        return std::get<elem_by_type<T>::index>(iter);        
    }

    iterator_tuple_type iter;
    friend class boost::iterator_core_access;
};

template<template<typename> class iterator_traits,typename tuple_of_iterators>
struct zip_helper {};

template <template<typename> class iterator_traits, typename ... T>
struct zip_helper<iterator_traits, std::tuple<T ...>> {
    typedef std::tuple<T...> tuple_iterator_type; 
    typedef std::tuple<typename iterator_traits<T>::value_type ...> value_type; 
    typedef std::tuple<typename iterator_traits<T>::reference ...> reference; 
};

template <template<typename> class iterator_traits, typename tuple_type, typename mpl_vector_type>
zip_iterator<typename zip_helper<iterator_traits,tuple_type>::tuple_iterator_type, 
             typename zip_helper<iterator_traits,tuple_type>::value_type, 
             typename zip_helper<iterator_traits,tuple_type>::reference,
             mpl_vector_type> 
make_zip_iterator(tuple_type arg) {
    return zip_iterator<typename zip_helper<iterator_traits,tuple_type>::tuple_iterator_type, 
             typename zip_helper<iterator_traits,tuple_type>::value_type, 
             typename zip_helper<iterator_traits,tuple_type>::reference, 
             mpl_vector_type>(arg);
}

//
// Particle getters/setters
//

/// get a variable from a particle \p arg
/// \param T Variable type
/// \param arg the particle
/// \return a const reference to a T::value_type holding the variable data
template<typename T, typename value_type>
typename value_type::template return_type<T>::type get(const value_type& arg) {
    return arg.template get<T>();
}

/// set a variable attached to a particle \p arg
/// \param T Variable type
/// \param arg the particle
/// \param data a reference to a T::value_type. This will be copied to
/// the variable attached to particle \p arg
template<typename T, typename value_type>
void set(const value_type& arg, const typename T::value_type & data) {
    arg.template get<T>() = data;
}



}


#endif //GET_H_
