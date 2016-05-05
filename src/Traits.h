
#ifndef TRAITS_H_ 
#define TRAITS_H_ 

#include "Vector.h"

namespace Aboria {

template<template<typename,typename> class VECTOR>
struct Traits {};

template <>
struct Traits<std::vector> {
    template <typename T>
    struct vector_type {
        typedef std::vector<T> type;
    };
    template <typename ... TupleTypes >
    struct tuple_type {
        typedef std::tuple<TupleTypes...> type;
    };
    template <typename T>
    struct zip_type {
        typedef boost::zip_iterator<T> type;
    };
};

#ifdef HAVE_THRUST
template <>
struct Traits<thrust::device_vector> {
    template <typename T>
    struct vector_type {
        typedef thrust::device_vector<T> type;
    };
    template <typename ... TupleTypes >
    struct tuple_type {
        typedef trust::tuple<TupleTypes...> type;
    };
    template <typename T>
    struct zip_type {
        typedef thrust::zip_iterator<T> type;
    };

    };
#endif

template<typename ARG,unsigned int D, typename TRAITS>
struct TraitsCommon {};

template <typename traits, unsigned int D, typename ... TYPES>
struct TraitsCommon<std::tuple<TYPES...>,D,traits> {

    const static unsigned int dimension = D;
    typedef typename traits::template vector_type<Vector<double,D> >::type vector_double_d;
    typedef typename traits::template vector_type<Vector<int,D> >::type vector_int_d;
    typedef typename traits::template vector_type<Vector<bool,D> >::type vector_bool_d;
    typedef typename traits::template vector_type<int>::type vector_int;
    typedef typename traits::template vector_type<int>::type vector_unsigned_int;
    typedef typename traits::template vector_type<Vector<int,2> >::type vector_int2;

    typedef Vector<double,dimension> double_d;
    typedef Vector<int,dimension> int_d;
    typedef Vector<unsigned int,dimension> unsigned_int_d;
    typedef Vector<bool,dimension> bool_d;

    ABORIA_VARIABLE(position,double_d,"position")
    ABORIA_VARIABLE(alive,bool,"is alive")
    ABORIA_VARIABLE(id,size_t,"id")

    typedef std::tuple<TYPES...> variable_type;
    typedef traits traits_type;
    typedef mpl::vector<position,id,alive,TYPES...> mpl_type_vector;
    typedef typename traits::template tuple_type<
        traits::template vector_type<position::value_type>::type::iterator,
        traits::template vector_type<id::value_type>::type::iterator,
        traits::template vector_type<alive::value_type>::type::iterator,
        typename traits::template vector_type<typename TYPES::value_type>::type::iterator...
            >::type tuple_of_iterators_type;

    typedef typename traits::template tuple_type<
        traits::template vector_type<position::value_type>::type::const_iterator,
        traits::template vector_type<id::value_type>::type::const_iterator,
        traits::template vector_type<alive::value_type>::type::const_iterator,
        typename traits::template vector_type<typename TYPES::value_type>::type::const_iterator...
            >::type tuple_of_const_iterators_type;

    typedef typename traits::template zip_type<tuple_of_iterators_type>::type iterator;
    typedef typename traits::template zip_type<tuple_of_const_iterators_type>::type const_iterator;

    typedef typename traits::template tuple_type< 
        typename position::value_type&,
        typename id::value_type&,
        typename alive::value_type&,
        typename TYPES::value_type&...
            >::type value_type;

    typedef typename traits::template tuple_type<
        traits::template vector_type<position::value_type>::type::iterator,
        traits::template vector_type<id::value_type>::type::iterator,
        traits::template vector_type<alive::value_type>::type::iterator,
        typename traits::template vector_type<typename TYPES::value_type>::type...
            >::type data_type;
};

#define UNPACK_TRAITS(traits)      \
    const static unsigned int dimension = traits::dimension;                          \
    typedef typename traits::vector_double_d vector_double_d;                           \
    typedef typename vector_double_d::const_iterator vector_double_d_const_iterator;    \
    typedef typename traits::vector_int_d vector_int_d;                                 \
    typedef typename traits::vector_unsigned_int_d vector_unsigned_int_d;               \
    typedef typename traits::vector_bool_d vector_bool_d;                               \
    typedef typename traits::vector_int2 vector_int2;                               \
    typedef Vector<double,dimension> double_d;                                        \
    typedef Vector<int,dimension> int_d;                                              \
    typedef Vector<unsigned int,dimension> unsigned_int_d;                            \
    typedef Vector<bool,dimension> bool_d;                                            \
    typedef Vector<int,2> int2;                                            \
    typedef typename traits::vector_int vector_int;                                     \
    typedef typename traits::vector_unsigned_int vector_unsigned_int;                   \
    typedef typename vector_unsigned_int::iterator vector_unsigned_int_iterator;        \
    typedef typename traits::iterator particles_iterator;                               \
    typedef typename traits::const_iterator const_particles_iterator;                   \
    typedef typename traits::value_type particles_value_type;                           \
    typedef typename traits::position position;                                         \
    typedef typename traits::id id;                                                     \
    typedef typename traits::alive alive;                                               \


}


#endif //TRAITS_H_
