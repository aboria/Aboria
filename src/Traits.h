
#ifndef TRAITS_H_ 
#define TRAITS_H_ 

#include "Vector.h"
#include <tuple>


namespace Aboria {

template<template<typename,typename> class VECTOR>
struct Traits {};

template <>
struct Traits<std::vector> {
    template <typename T>
    struct vector_type {
        typedef std::vector<T> type;
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
struct TraitsCommon {
    typedef typename ARG::ERROR_FIRST_TEMPLATE_ARGUMENT_TO_PARTICLES_MUST_BE_A_STD_TUPLE_TYPE error; 
    };

template <typename traits, unsigned int D, typename ... TYPES>
struct TraitsCommon<std::tuple<TYPES...>,D,traits> {

    template <typename ... TupleTypes >
    struct tuple_type {
        typedef boost::tuple<TupleTypes...> type;
    };
    template <typename T>
    struct zip_type {
        typedef boost::zip_iterator<T> type;
    };

    const static unsigned int dimension = D;
    typedef typename traits::template vector_type<Vector<double,D> >::type vector_double_d;
    typedef typename traits::template vector_type<Vector<int,D> >::type vector_int_d;
    typedef typename traits::template vector_type<Vector<unsigned int,D> >::type vector_unsigned_int_d;
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
    typedef typename position::value_type position_value_type;
    typedef typename alive::value_type alive_value_type;
    typedef typename id::value_type id_value_type;
    typedef typename traits::template vector_type<position_value_type>::type position_vector_type;
    typedef typename traits::template vector_type<alive_value_type>::type alive_vector_type;
    typedef typename traits::template vector_type<id_value_type>::type id_vector_type;

    typedef std::tuple<TYPES...> variable_type;
    typedef traits traits_type;
    typedef mpl::vector<position,id,alive,TYPES...> mpl_type_vector;
    typedef typename tuple_type<
        typename position_vector_type::iterator,
        typename id_vector_type::iterator,
        typename alive_vector_type::iterator,
        typename traits::template vector_type<typename TYPES::value_type>::type::iterator...
            >::type tuple_of_iterators_type;

    typedef typename tuple_type<
        typename position_vector_type::const_iterator,
        typename id_vector_type::const_iterator,
        typename alive_vector_type::const_iterator,
        typename traits::template vector_type<typename TYPES::value_type>::type::const_iterator...
            >::type tuple_of_const_iterators_type;

    typedef typename tuple_type< 
        const typename position::value_type&,
        const typename id::value_type&,
        const typename alive::value_type&,
        const typename TYPES::value_type&...
            >::type const_reference_data_type;

    typedef typename tuple_type< 
        typename position::value_type&,
        typename id::value_type&,
        typename alive::value_type&,
        typename TYPES::value_type&...
            >::type reference_data_type;

    typedef typename tuple_type< 
        typename position::value_type,
        typename id::value_type,
        typename alive::value_type,
        typename TYPES::value_type...
            >::type value_data_type;

    typedef std::tuple<
        position_vector_type,
        id_vector_type,
        alive_vector_type,
        typename traits::template vector_type<typename TYPES::value_type>::type...
            > vectors_data_type;



    typedef Aboria::value_type<value_data_type, mpl_type_vector> value_type;
    typedef Aboria::tuple_of_vectors_type<vectors_data_type, mpl_type_vector> data_type;
    typedef Aboria::reference_type<reference_data_type, mpl_type_vector> reference;
    typedef Aboria::reference_type<const_reference_data_type, mpl_type_vector> const_reference;

    typedef typename Aboria::zip_iterator<
        typename zip_type<tuple_of_iterators_type>::type, 
        tuple_of_iterators_type,reference,mpl_type_vector> iterator;
    typedef typename Aboria::zip_iterator<
        typename zip_type<tuple_of_const_iterators_type>::type, 
        tuple_of_iterators_type,const_reference,mpl_type_vector> const_iterator;

    /*
    typedef typename tuple_type<
        position_vector_type,
        id_vector_type,
        alive_vector_type,
        typename traits::template vector_type<typename TYPES::value_type>::type...
            >::type data_type;
            */

    

    const static size_t N = std::tuple_size<vectors_data_type>::value;

    template<std::size_t... I>
    static iterator begin_impl(data_type& data, index_sequence<I...>) {
        return boost::make_zip_iterator(boost::make_tuple(data.get<I>.begin()...));
    }

    template<std::size_t... I>
    static iterator end_impl(data_type& data, index_sequence<I...>) {
        return boost::make_zip_iterator(boost::make_tuple(data.get<I>.end()...));
    }

    template<std::size_t... I>
    static value_type index_impl(data_type& data, const size_t i, index_sequence<I...>) {
        return value_type(std::get<I>(data.get_tuple())[i]...);
    }

    template<std::size_t... I>
    static void clear_impl(data_type& data, index_sequence<I...>) {
        using expander = int[];
        (void)expander { 0, (std::get<I>(data.get_tuple()).clear())... };
    }

    template<std::size_t... I>
    static void push_back_impl(data_type& data, const value_type& val, index_sequence<I...>) {
        int dummy[] = { 0, (std::get<I>(data.get_tuple()).push_back(val.get_tuple().get<I>()),void(),0)... };
        static_cast<void>(dummy); // Avoid warning for unused variable.
    }

    template<std::size_t... I>
    static void pop_back_impl(data_type& data, index_sequence<I...>) {
        int dummy[] = { 0, (std::get<I>(data.get_tuple()).pop_back())... };
        static_cast<void>(dummy); // Avoid warning for unused variable.
    }

    template<std::size_t... I>
    static void insert_impl(data_type& data, iterator position, const value_type& val, index_sequence<I...>) {
        using expander = int[];
        (void)expander { 0, (std::get<I>(data.get_tuple()).insert(position,val.get_tuple().get<I>()))... };
    }

    template<std::size_t... I>
    static void insert_impl(data_type& data, iterator position, size_t n, const value_type& val, index_sequence<I...>) {
        using expander = int[];
        (void)expander { 0, (std::get<I>(data.get_tuple()).insert(position,n,val.get_tuple().get<I>()))... };
    }

    template<class InputIterator, std::size_t... I>
    static void insert_impl(data_type& data, iterator position, InputIterator first, InputIterator last, index_sequence<I...>) {
        using expander = int[];
        (void)expander { 0, (std::get<I>(data.get_tuple()).insert(position,first,last))... };
    }

    template<typename Indices = make_index_sequence<N>>
    static iterator begin(data_type& data) {
        return begin_impl(data, Indices());
    }

    template<typename Indices = make_index_sequence<N>>
    static iterator end(data_type& data) {
        return end_impl(data, Indices());
    }

    template<typename Indices = make_index_sequence<N>>
    static value_type index(data_type& data, const size_t i) {
        return index_impl(data, i, Indices());
    }

    template<typename Indices = make_index_sequence<N>>
    static void clear(data_type& data) {
        clear_impl(data, Indices());
    }

    template<typename Indices = make_index_sequence<N>>
    static void push_back(data_type& data, const value_type& val) {
        push_back_impl(data, val, Indices());
    }

    template<typename Indices = make_index_sequence<N>>
    static void pop_back(data_type& data) {
        pop_back_impl(data, Indices());
    }

    template<typename Indices = make_index_sequence<N>>
    static void insert(data_type& data, iterator pos, const value_type& val) {
        insert_impl(data, val, Indices());
    }

    template<typename Indices = make_index_sequence<N>>
    static void insert(data_type& data, iterator position, size_t n,  const value_type& val) {
        insert_impl(data, val, n, Indices());
    }

    template<class InputIterator, typename Indices = make_index_sequence<N>>
    static void insert(data_type& data, iterator pos, InputIterator first, InputIterator last) {
        insert_impl(data, pos, first, last, Indices());
    }



    typedef typename position_vector_type::size_type size_type; 
    typedef typename position_vector_type::difference_type difference_type; 
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
