
#ifndef GET_H_ 
#define GET_H_ 

#include <boost/mpl/vector.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/range_c.hpp>

namespace Aboria {

namespace mpl = boost::mpl;
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

template<typename ARG>
struct unpack_tuple_types {};

template <typename ... TYPES>
struct unpack_tuple_types< std::tuple<TYPES...> > {
    typedef mpl::vector<typename std::remove_reference<TYPES>::type...> mpl_type_vector;
};


//
// Particle getters/setters
//

/// get a variable from a particle \p arg
/// \param T Variable type
/// \param arg the particle
/// \return a const reference to a T::value_type holding the variable data
template<typename T, typename value_type>
const typename T::value_type & get(const value_type& arg) {
    typedef typename unpack_tuple_types<value_type>::mpl_type_vector mpl_type_vector;
    return std::get<elem_by_type<T,mpl_type_vector>::index>(arg);
}

/// get a variable from a particle \p arg
/// \param T Variable type
/// \param arg the particle
/// \return a reference to a T::value_type holding the variable data
template<typename T, typename value_type>
typename T::value_type & get(value_type& arg) {
    typedef typename unpack_tuple_types<value_type>::mpl_type_vector mpl_type_vector;
    return std::get<elem_by_type<T,mpl_type_vector>::index>(arg);
}

/// set a variable attached to a particle \p arg
/// \param T Variable type
/// \param arg the particle
/// \param data a reference to a T::value_type. This will be copied to
/// the variable attached to particle \p arg
template<typename T, typename value_type>
void set(value_type& arg, const typename T::value_type & data) {
    typedef typename unpack_tuple_types<value_type>::mpl_type_vector mpl_type_vector;
    std::get<elem_by_type<T,mpl_type_vector>::index>(arg) = data;
}



}


#endif //GET_H_
