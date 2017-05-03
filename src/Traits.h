
#ifndef TRAITS_H_ 
#define TRAITS_H_ 

#include "Vector.h"
#include "CudaInclude.h"
#include "Get.h"
#include <tuple>
#include <vector>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/lambda/lambda.hpp>
#include <algorithm>


namespace Aboria {

namespace mpl = boost::mpl;


struct default_traits {
    template <typename T>
    struct vector_type {
        typedef std::vector<T> type;
    };

};

template<template<typename,typename> class VECTOR>
struct Traits {};

template <>
struct Traits<std::vector>: public default_traits {};

#if defined(__CUDACC__)
template <>
struct Traits<thrust::device_vector>: public default_traits {
    template <typename T>
    struct vector_type {
        typedef thrust::device_vector<T> type;
    };


    template<typename ForwardIterator, typename InputIterator, typename OutputIterator>
    static void lower_bound(
            ForwardIterator first,
            ForwardIterator last,
            InputIterator values_first,
            InputIterator values_last,
            OutputIterator result) {
        thrust::lower_bound(first,last,values_first,values_last,result);
    }

    template<typename ForwardIterator, typename InputIterator, typename OutputIterator>
    static void upper_bound(
            ForwardIterator first,
            ForwardIterator last,
            InputIterator values_first,
            InputIterator values_last,
            OutputIterator result) {
        thrust::upper_bound(first,last,values_first,values_last,result);
    }

    template< class InputIt, class T, class BinaryOperation >
    static T reduce( 
        InputIt first, 
        InputIt last, T init,
        BinaryOperation op ) {
        return thrust::reduce(first,last,init,op);
    }

    template <class InputIterator, class OutputIterator, class UnaryOperation>
    static OutputIterator transform (
            InputIterator first, InputIterator last,
            OutputIterator result, UnaryOperation op) {
    }

    template <class ForwardIterator>
    static void sequence (ForwardIterator first, ForwardIterator last) {
    }

    template<typename ForwardIterator , typename UnaryOperation >
    static void tabulate (
            ForwardIterator first,
            ForwardIterator last,
            UnaryOperation  unary_op) {	

    }

    template<typename InputIterator, typename OutputIterator, typename UnaryFunction, 
        typename T, typename AssociativeOperator>
    static OutputIterator transform_exclusive_scan(
        InputIterator first, InputIterator last,
        OutputIterator result,
        UnaryFunction unary_op, T init, AssociativeOperator binary_op) {
    }


    template<typename InputIterator1, typename InputIterator2, 
        typename InputIterator3, typename RandomAccessIterator , typename Predicate >
    static void scatter_if(
            InputIterator1 first, InputIterator1 last,
            InputIterator2 map, InputIterator3 stencil,
            RandomAccessIterator output, Predicate pred) {
    }

    template<typename InputIterator1, typename InputIterator2, 
        typename OutputIterator, typename Predicate>
    static OutputIterator copy_if(
            InputIterator1 first, InputIterator1 last, 
            InputIterator2 stencil, OutputIterator result, Predicate pred) {
    }
};
#endif

template<typename ARG,unsigned int D, typename TRAITS>
struct TraitsCommon {
    typedef typename ARG::ERROR_FIRST_TEMPLATE_ARGUMENT_TO_PARTICLES_MUST_BE_A_STD_TUPLE_TYPE error; 
    };

template <typename traits, unsigned int D, typename ... TYPES>
struct TraitsCommon<std::tuple<TYPES...>,D,traits>:public traits {

    const static unsigned int dimension = D;
    typedef typename traits::template vector_type<Vector<double,D> >::type vector_double_d;
    typedef typename vector_double_d::iterator vector_double_d_iterator;
    typedef typename vector_double_d::const_iterator vector_double_d_const_iterator;
    typedef typename traits::template vector_type<Vector<int,D> >::type vector_int_d;
    typedef typename traits::template vector_type<Vector<unsigned int,D> >::type vector_unsigned_int_d;
    typedef typename vector_unsigned_int_d::iterator vector_unsigned_int_d_iterator;
    typedef typename vector_unsigned_int_d::const_iterator vector_unsigned_int_d_const_iterator;
    typedef typename traits::template vector_type<Vector<bool,D> >::type vector_bool_d;
    typedef typename traits::template vector_type<int>::type vector_int;
    typedef typename traits::template vector_type<unsigned int>::type vector_unsigned_int;
    typedef typename traits::template vector_type<size_t>::type vector_size_t;
    typedef typename vector_unsigned_int::iterator vector_unsigned_int_iterator;
    typedef typename vector_unsigned_int::const_iterator vector_unsigned_int_const_iterator;
    typedef typename traits::template vector_type<Vector<int,2> >::type vector_int2;

    typedef Vector<double,dimension> double_d;
    typedef Vector<int,dimension> int_d;
    typedef Vector<unsigned int,dimension> unsigned_int_d;
    typedef Vector<bool,dimension> bool_d;

    typedef position_d<dimension> position;
    typedef typename position::value_type position_value_type;
    typedef alive::value_type alive_value_type;
    typedef id::value_type id_value_type;
    typedef random::value_type random_value_type;
    typedef typename traits::template vector_type<position_value_type>::type position_vector_type;
    typedef typename traits::template vector_type<alive_value_type>::type alive_vector_type;
    typedef typename traits::template vector_type<id_value_type>::type id_vector_type;
    typedef typename traits::template vector_type<random_value_type>::type random_vector_type;

    typedef traits traits_type;
    typedef mpl::vector<position,id,alive,random,TYPES...> mpl_type_vector;

    /*
    typedef tuple_ns::tuple<
            typename position_vector_type::pointer,
            typename id_vector_type::pointer,
            typename alive_vector_type::pointer,
            typename traits::template vector_type<typename TYPES::value_type>::type::pointer...
            > tuple_of_pointers_type;

    typedef tuple_ns::tuple<
            typename position_vector_type::value_type*,
            typename id_vector_type::value_type*,
            typename alive_vector_type::value_type*,
            typename traits::template vector_type<typename TYPES::value_type>::type::value_type*...
            > tuple_of_raw_pointers_type;
            */

    typedef tuple_ns::tuple<
            typename position_vector_type::iterator,
            typename id_vector_type::iterator,
            typename alive_vector_type::iterator,
            typename random_vector_type::iterator,
            typename traits::template vector_type<typename TYPES::value_type>::type::iterator...
            > tuple_of_iterators_type;

    typedef tuple_ns::tuple<
            typename position_vector_type::const_iterator,
            typename id_vector_type::const_iterator,
            typename alive_vector_type::const_iterator,
            typename random_vector_type::const_iterator,
            typename traits::template vector_type<typename TYPES::value_type>::type::const_iterator...
            > tuple_of_const_iterators_type;


    typedef tuple_ns::tuple<
        position_vector_type,
        id_vector_type,
        alive_vector_type,
        random_vector_type,
        typename traits::template vector_type<typename TYPES::value_type>::type...
            > vectors_data_type;

    typedef typename Aboria::zip_iterator<tuple_of_iterators_type,mpl_type_vector> iterator;
    typedef typename Aboria::zip_iterator<tuple_of_const_iterators_type,mpl_type_vector> const_iterator;

    typedef typename iterator::reference reference;
    typedef typename iterator::value_type value_type;
    typedef typename iterator::pointer pointer;
    typedef typename iterator::tuple_raw_pointer raw_pointer;
    typedef typename iterator::tuple_raw_reference raw_reference;
    typedef typename const_iterator::reference const_reference;
    typedef Aboria::getter_type<vectors_data_type, mpl_type_vector> data_type;

    const static size_t N = tuple_ns::tuple_size<vectors_data_type>::value;

    template<std::size_t... I>
    static iterator begin_impl(data_type& data, detail::index_sequence<I...>) {
        return iterator(get_by_index<I>(data).begin()...);
    }

    template<std::size_t... I>
    static iterator end_impl(data_type& data, detail::index_sequence<I...>) {
        return iterator(get_by_index<I>(data).end()...);
    }

    template<std::size_t... I>
    static const_iterator cbegin_impl(const data_type& data, detail::index_sequence<I...>) {
        return const_iterator(get_by_index<I>(data).cbegin()...);
    }

    template<std::size_t... I>
    static const_iterator cend_impl(const data_type& data, detail::index_sequence<I...>) {
        return const_iterator(get_by_index<I>(data).cend()...);
    }

    template<std::size_t... I>
    static reference index_impl(data_type& data, const size_t i, detail::index_sequence<I...>, std::true_type) {
        return reference(get_by_index<I>(data)[i]...);
    }

    template<std::size_t... I>
    static reference index_impl(data_type& data, const size_t i, detail::index_sequence<I...>, std::false_type) {
        return reference(get_by_index<I>(data)[i]...);
    }

    template<std::size_t... I>
    static const_reference index_const_impl(const data_type& data, const size_t i, detail::index_sequence<I...>, std::true_type) {
        return const_reference(get_by_index<I>(data)[i]...);
    }

    template<std::size_t... I>
    static const_reference index_const_impl(const data_type& data, const size_t i, detail::index_sequence<I...>, std::false_type) {
        return const_reference(get_by_index<I>(data)[i]...);
    }

    template<std::size_t... I>
    static void clear_impl(data_type& data, detail::index_sequence<I...>) {
        int dummy[] = { 0, (get_by_index<I>(data).clear(),void(),0)... };
        static_cast<void>(dummy);
    }

    template<std::size_t... I>
    static void resize_impl(data_type& data, const size_t new_size, detail::index_sequence<I...>) {
        int dummy[] = { 0, (get_by_index<I>(data).resize(new_size),void(),0)... };
        static_cast<void>(dummy);
    }

    template<std::size_t... I>
    static void push_back_impl(data_type& data, const value_type& val, detail::index_sequence<I...>) {
        int dummy[] = { 0, (get_by_index<I>(data).push_back(get_by_index<I>(val)),void(),0)... };
        static_cast<void>(dummy); // Avoid warning for unused variable.
    }

    template<std::size_t... I>
    static void pop_back_impl(data_type& data, detail::index_sequence<I...>) {
        int dummy[] = { 0, (get_by_index<I>(data).pop_back(),void(),0)... };
        static_cast<void>(dummy); // Avoid warning for unused variable.
    }

    template<std::size_t... I>
    static void insert_impl(data_type& data, iterator position, const value_type& val, detail::index_sequence<I...>) {
        int dummy[] = { 0, (get_by_index<I>(data).insert(position,get_by_index<I>(val)))... };
        static_cast<void>(dummy);
    }

    template<std::size_t... I>
    static void insert_impl(data_type& data, iterator position, size_t n, const value_type& val, detail::index_sequence<I...>) {
        int dummy[] = { 0, (get_by_index<I>(data).insert(position,n,get_by_index<I>(val)))... };
        static_cast<void>(dummy);
    }

    template<class InputIterator, std::size_t... I>
    static void insert_impl(data_type& data, iterator position, InputIterator first, InputIterator last, detail::index_sequence<I...>) {
        int dummy[] = { 0, (get_by_index<I>(data).insert(position,first,last))... };
        static_cast<void>(dummy);
    }

    template<typename Indices = detail::make_index_sequence<N>>
    static iterator begin(data_type& data) {
        return begin_impl(data, Indices());
    }

    template<typename Indices = detail::make_index_sequence<N>>
    static iterator end(data_type& data) {
        return end_impl(data, Indices());
    }

    template<typename Indices = detail::make_index_sequence<N>>
    static const_iterator cbegin(const data_type& data) {
        return cbegin_impl(data, Indices());
    }

    template<typename Indices = detail::make_index_sequence<N>>
    static const_iterator cend(const data_type& data) {
        return cend_impl(data, Indices());
    }

    template<typename Indices = detail::make_index_sequence<N>>
    static reference index(data_type& data, const size_t i) {
        return index_impl(data, i, Indices(),std::is_reference<decltype(get<id>(data)[0])>());
    }

    template<typename Indices = detail::make_index_sequence<N>>
    static const_reference index(const data_type& data, const size_t i) {
        return index_const_impl(data, i, Indices(),std::is_reference<decltype(get<id>(data)[0])>());
    }

    template<typename Indices = detail::make_index_sequence<N>>
    static void clear(data_type& data) {
        clear_impl(data, Indices());
    }

    template<typename Indices = detail::make_index_sequence<N>>
    static void resize(data_type& data, const size_t new_size) {
        resize_impl(data, new_size, Indices());
    }

    template<typename Indices = detail::make_index_sequence<N>>
    static void push_back(data_type& data, const value_type& val) {
        push_back_impl(data, val, Indices());
    }

    template<typename Indices = detail::make_index_sequence<N>>
    static void pop_back(data_type& data) {
        pop_back_impl(data, Indices());
    }

    template<typename Indices = detail::make_index_sequence<N>>
    static void insert(data_type& data, iterator pos, const value_type& val) {
        insert_impl(data, val, Indices());
    }

    template<typename Indices = detail::make_index_sequence<N>>
    static void insert(data_type& data, iterator position, size_t n,  const value_type& val) {
        insert_impl(data, val, n, Indices());
    }

    template<class InputIterator, typename Indices = detail::make_index_sequence<N>>
    static void insert(data_type& data, iterator pos, InputIterator first, InputIterator last) {
        insert_impl(data, pos, first, last, Indices());
    }



    typedef typename position_vector_type::size_type size_type; 
    typedef typename position_vector_type::difference_type difference_type; 
};


}


#endif //TRAITS_H_
