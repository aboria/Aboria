
#ifndef GET_H_ 
#define GET_H_ 

#include <boost/mpl/vector.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/range_c.hpp>
#include "boost/tuple/tuple.hpp"

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
struct value_type {
        typedef mpl_vector_type mpl_type_vector;
        template <typename T>
        using elem_by_type = elem_by_type<T,mpl_type_vector>;
        template <typename T>
        using return_type = element<elem_by_type<T>::index,tuple_type>;

        value_type() {}
        value_type(const tuple_type& data):data(data) {}

        template<typename T>
        typename return_type<T>::type get() const {
            return data.template get<elem_by_type<T>::index>();
        }

        const tuple_type & get_tuple() const { return data; }
        tuple_type & get_tuple() { return data; }
                
        tuple_type data;
};

template <typename tuple_type, typename mpl_vector_type> 
struct reference_type {
        typedef mpl_vector_type mpl_type_vector;
        template <typename T>
        using elem_by_type = elem_by_type<T,mpl_type_vector>;
        template <typename T>
        using return_type = element<elem_by_type<T>::index,tuple_type>;

        reference_type() {}
        reference_type(tuple_type data):data(data) {}

        template<typename T>
        typename return_type<T>::type get() const {
            return data.template get<elem_by_type<T>::index>();
        }

        const tuple_type & get_tuple() const { return data; }
        tuple_type & get_tuple() { return data; }
                
        tuple_type data;
};

template <typename zip_type, typename iterator_tuple_type, typename reference_type, typename mpl_vector_type>
struct zip_iterator {
    typedef zip_iterator<zip_type,iterator_tuple_type,reference_type,mpl_vector_type> my_type;
    typedef mpl_vector_type mpl_type_vector;
    template <typename T>
    using elem_by_type = elem_by_type<T,mpl_type_vector>;
    template <typename T>
    using return_type = element<elem_by_type<T>::index,iterator_tuple_type>;


    typedef reference_type reference;
    zip_iterator() {}
    zip_iterator(iterator_tuple_type iterators): zip(iterators) {}
    zip_iterator(const zip_type& other): zip(other) {}

    reference operator*() const {
        return reference(*zip);
    }
    zip_iterator& operator++() {
        ++zip;
        return *this;
    }
    zip_iterator operator++(int) {
        return my_type(zip++);
    }
    zip_iterator& operator--() {
        zip--;
        return *this;
    }
    zip_iterator operator-(size_t n) {
        return my_type(zip-n);
    }
    size_t operator-(const my_type& other) const {
        return zip-other.zip;
    }
    zip_iterator operator+(size_t n) {
        return my_type(zip+n);
    }
    bool operator!=(const my_type& rhs) { 
        return (zip != rhs.zip); 
    }
    bool operator==(const my_type& rhs) { 
        return (zip == rhs.zip); 
    }

    template<typename T>
    typename return_type<T>::type get() const {
        return zip.get_iterator_tuple().template get<elem_by_type<T>::index>();
    }

    zip_type zip;
};


template <typename tuple_type, typename mpl_vector_type> 
struct tuple_of_vectors_type {
        typedef mpl_vector_type mpl_type_vector;
        template <typename T>
        using elem_by_type = elem_by_type<T,mpl_type_vector>;
        template <typename T>
        using return_type = element<elem_by_type<T>::index,tuple_type>;

        tuple_of_vectors_type() {}
        tuple_of_vectors_type(const tuple_type& data):data(data) {}

        template<typename T>
        typename return_type<T>::type & get() const {
            return std::get<elem_by_type<T>::index>(data);        
        }

        const tuple_type & get_tuple() const { return data; }
        tuple_type & get_tuple() { return data; }

                
        tuple_type data;
};



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
