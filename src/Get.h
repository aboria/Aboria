
#ifndef GET_H_ 
#define GET_H_ 

#include <boost/mpl/vector.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <algorithm>
#include <tuple>
#include <type_traits>



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
    typedef typename boost::tuples::element<N,tuple_type>::type raw_type;
    typedef typename std::remove_reference<raw_type>::type type;
    typedef const type const_type;
};

template<size_t N, typename... tuple_args>
struct element<N,std::tuple<tuple_args...>> {
    typedef typename std::tuple_element<N,std::tuple<tuple_args...>>::type raw_type;
    typedef typename std::remove_reference<raw_type>::type type;
    typedef const type const_type;
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


/*
template<typename T>
struct getter_reference {
    typedef T& type;
};

template <typename tuple_type, typename mpl_vector_type> 
struct getter_reference<getter_type<tuple_type,mpl_vector_type>> {
    typedef typename getter_type<tuple_type,mpl_vector_type>::reference type;
};

template<typename tuple_type>
struct tuple_helper {};

template <typename ... T>
struct tuple_helper<std::tuple<T ...>> {
    typedef std::tuple<T...> value_type; 
    typedef std::tuple<typename getter_reference<T>::type...> reference; 
};
*/


template <typename tuple_type, typename mpl_vector_type> 
struct getter_type{};

template <typename mpl_vector_type, typename ... tuple_elements> 
struct getter_type<std::tuple<tuple_elements...>,mpl_vector_type> {
    typedef std::tuple<tuple_elements...> tuple_type;

    typedef mpl_vector_type mpl_type_vector;
    template <typename T>
    using elem_by_type = elem_by_type<T,mpl_type_vector>;
    template <typename T>
    using return_type = element<elem_by_type<T>::index,tuple_type>;

    //typedef getter_type<typename tuple_helper<tuple_type>::reference,mpl_type_vector> reference;

    getter_type() {}
    explicit getter_type(const tuple_type& data):data(data) {}
    getter_type(const getter_type& other):data(other.data) {}
    getter_type(getter_type&& other):data(std::move(other.data)) {}
    //getter_type(const reference& other):data(other.data) {}

    /*
    template<typename... _UElements, typename = typename
        enable_if<__and_<is_convertible<const _UElements&,
					_Elements>...>::value>::type>
        constexpr tuple(const tuple<_UElements...>& __in)
        : _Inherited(static_cast<const _Tuple_impl<0, _UElements...>&>(__in))
        { }

      template<typename... _UElements, typename = typename
        enable_if<__and_<is_convertible<_UElements,
					_Elements>...>::value>::type>
        constexpr tuple(tuple<_UElements...>&& __in)
        : _Inherited(static_cast<_Tuple_impl<0, _UElements...>&&>(__in)) { }
        */



    template <typename T1, typename... T2,typename = typename
	std::enable_if<std::__and_<std::is_convertible<const T2&, tuple_elements>...>::value>::type>
    getter_type(const getter_type<std::tuple<T2...>,T1>& other):data(other.data) {}

    template <typename T1, typename... T2,typename = typename
	std::enable_if<std::__and_<std::is_convertible<T2, tuple_elements>...>::value>::type>
    getter_type(const getter_type<std::tuple<T2...>,T1>&& other):data(std::move(other.data)) {}

    
    getter_type& operator=( const getter_type& other ) {
        data = other.data;
        return *this;
    }
    getter_type& operator=( getter_type&& other ) {
        data = std::move(other.data);
        return *this;
    }
    template <typename T1, typename T2> 
    getter_type& operator=( const getter_type<T1,T2>& other) {
        data = other.data;
        return *this;
    }
    template <typename T1, typename T2> 
    getter_type& operator=( getter_type<T1,T2>&& other) {
        data = std::move(other.data);
        return *this;
    }
    
    void swap(getter_type &other) {
        data.swap(other.data);
    }

    template<typename T>
    typename return_type<T>::const_type & get() const {
        return std::get<elem_by_type<T>::index>(data);        
    }

    template<typename T>
    typename return_type<T>::type & get() {
        return std::get<elem_by_type<T>::index>(data);        
    }

    const tuple_type & get_tuple() const { return data; }
    tuple_type & get_tuple() { return data; }
            
    tuple_type data;
};


template <typename tuple_type, typename mpl_vector_type> 
void swap(getter_type<tuple_type,mpl_vector_type> x,
          getter_type<tuple_type,mpl_vector_type> y) {
    x.swap(y);
}




template<typename tuple_of_iterators>
struct zip_helper {};

template <typename ... T>
struct zip_helper<std::tuple<T ...>> {
    typedef std::tuple<T...> tuple_iterator_type; 
    typedef std::tuple<typename std::iterator_traits<T>::value_type ...> tuple_value_type; 
    typedef std::tuple<typename std::iterator_traits<T>::reference ...> tuple_reference; 
};


template <typename iterator_tuple_type, typename mpl_vector_type=mpl::vector<int>>
class zip_iterator: public boost::iterator_facade<zip_iterator<iterator_tuple_type,mpl_vector_type>,
    getter_type<typename zip_helper<iterator_tuple_type>::tuple_value_type,mpl_vector_type>,
    boost::random_access_traversal_tag,
    getter_type<typename zip_helper<iterator_tuple_type>::tuple_reference,mpl_vector_type>
        > {

    typedef mpl_vector_type mpl_type_vector;
    template <typename T>
    using elem_by_type = elem_by_type<T,mpl_type_vector>;

public:
    typedef getter_type<typename zip_helper<iterator_tuple_type>::tuple_value_type,mpl_vector_type> value_type;
    typedef getter_type<typename zip_helper<iterator_tuple_type>::tuple_reference,mpl_vector_type> reference;
        
    template<typename T>
    struct return_type {
        static const size_t N = elem_by_type<T>::index;
        typedef const typename std::tuple_element<N,iterator_tuple_type>::type type;
        typedef type const_type;
    };


    typedef typename element<0,iterator_tuple_type>::type first_type;
    zip_iterator() {}
    explicit zip_iterator(iterator_tuple_type iter) : iter(iter) {}
    typedef typename std::iterator_traits<first_type>::difference_type difference_type;


    template<typename T>
    const typename return_type<T>::type & get() const {
        return std::get<elem_by_type<T>::index>(iter);        
    }

private:
    typedef make_index_sequence<std::tuple_size<iterator_tuple_type>::value> index_type;

    template<std::size_t... I>
    static reference make_reference(const iterator_tuple_type& tuple, index_sequence<I...>) {
        typedef typename zip_helper<iterator_tuple_type>::tuple_reference tuple_reference;
        return reference(tuple_reference(*(std::get<I>(tuple))...));
    }
    template<std::size_t... I>
    static void increment_impl(iterator_tuple_type& tuple, index_sequence<I...>) {
        using expander = int[];
        (void)expander { 0, (++std::get<I>(tuple),0)...};
    }

    template<std::size_t... I>
    static void decrement_impl(iterator_tuple_type& tuple, index_sequence<I...>) {
        using expander = int[];
        (void)expander { 0, (--std::get<I>(tuple),void(),0)...};
    }

    template<std::size_t... I>
    static void advance_impl(iterator_tuple_type& tuple, const difference_type n,  index_sequence<I...>) {
        using expander = int[];
        (void)expander { 0, (std::get<I>(tuple) += n,void(),0)...};
    }

    void increment() { increment_impl(iter,index_type()); }
    void decrement() { decrement_impl(iter,index_type()); }

    bool equal(zip_iterator const& other) const { return std::get<0>(iter) == std::get<0>(other.iter); }
    reference dereference() const { return make_reference(iter,index_type()); }
    difference_type distance_to(zip_iterator const& other) const { return std::get<0>(other.iter) - std::get<0>(iter); }
    void advance(difference_type n) { advance_impl(iter,n,index_type()); }

    iterator_tuple_type iter;
    friend class boost::iterator_core_access;
};


//
// Particle getters/setters
//

/// get a variable from a particle \p arg
/// \param T Variable type
/// \param arg the particle
/// \return a const reference to a T::value_type holding the variable data
///
template<typename T, typename value_type>
typename value_type::template return_type<T>::const_type & get(const value_type& arg) {
    return arg.template get<T>();
}

template<typename T, typename value_type>
typename value_type::template return_type<T>::type & get(value_type& arg) {
    return arg.template get<T>();
}

/// set a variable attached to a particle \p arg
/// \param T Variable type
/// \param arg the particle
/// \param data a reference to a T::value_type. This will be copied to
/// the variable attached to particle \p arg
template<typename T, typename value_type>
void set(value_type arg, const typename T::value_type & data) {
    arg.template get<T>() = data;
}



}




#endif //GET_H_
