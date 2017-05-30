
#ifndef ALGORITHMS_H_ 
#define ALGORITHMS_H_ 

#include "CudaInclude.h"
#include "Get.h"
#include "Traits.h"
#include <algorithm>



namespace Aboria {

namespace detail {

// reorders v such that v_new[i] == v[order[i]]
template< typename order_iterator, typename value_iterator >
void reorder_destructive( order_iterator order_begin, order_iterator order_end, value_iterator v )  {
    typedef typename std::iterator_traits< value_iterator >::value_type value_t;
    typedef typename std::iterator_traits< value_iterator >::reference reference;
    typedef typename std::iterator_traits< order_iterator >::value_type index_t;
    typedef typename std::iterator_traits< order_iterator >::difference_type diff_t;

    diff_t size = order_end - order_begin;
    
    size_t i, j, k;
    value_t temp;
    for(i = 0; i < size; i++){
        if(i != order_begin[i]){
            temp = v[i];
            k = i;
            while(i != (j = order_begin[k])){
                // every move places a value in it's final location
                // NOTE: need the static cast or assignment has no effect (TODO: why?)
                static_cast<reference>(v[k])  = v[j];
                order_begin[k] = k;
                k = j;
            }
            v[k] = temp;
            order_begin[k] = k;
        }
    }
}


template <typename value_type> 
struct iter_comp {
    bool operator()(const value_type& t1, const value_type& t2) { return get<0>(t1.get_tuple()) < get<0>(t2.get_tuple()); }
};

template <typename T> 
struct lower_bound_impl {
    const T& values_first,values_last;
    lower_bound_impl(const T& values_first, const T& values_last):
        values_first(values_first),values_last(values_last) {}
    typedef typename std::iterator_traits<T>::value_type value_type;
    typedef typename std::iterator_traits<T>::difference_type difference_type;
    difference_type operator()(const value_type & search) { 
        return std::distance(values_first,std::lower_bound(values_first,values_last,search)); 
    }
};

template <typename T> 
struct upper_bound_impl {
    const T& values_first,values_last;
    upper_bound_impl(const T& values_first, const T& values_last):
        values_first(values_first),values_last(values_last) {}
    typedef typename std::iterator_traits<T>::value_type value_type;
    typedef typename std::iterator_traits<T>::difference_type difference_type;
    difference_type operator()(const value_type & search) { 
        return std::distance(values_first,std::upper_bound(values_first,values_last,search)); 
    }
};

#ifdef __aboria_use_thrust_algorithms__
template <typename T>
using counting_iterator = thrust::counting_iterator<T>;

static const thrust::detail::functional::placeholder<0>::type _1;
static const thrust::detail::functional::placeholder<1>::type _2;
static const thrust::detail::functional::placeholder<2>::type _3;

using thrust::make_transform_iterator;
using thrust::make_zip_iterator;
using thrust::make_tuple;
#else

template <typename T>
using counting_iterator = boost::counting_iterator<T>;

template <class UnaryFunction, class Iterator>
using transform_iterator = boost::transform_iterator<UnaryFunction, Iterator>;

const boost::lambda::placeholder1_type _1;
const boost::lambda::placeholder2_type _2;
const boost::lambda::placeholder3_type _3;

using boost::make_transform_iterator;
using boost::make_zip_iterator;
using boost::make_tuple;

//template <class UnaryFunction, class Iterator>
//transform_iterator<UnaryFunction, Iterator>
//make_transform_iterator(Iterator&& it, UnaryFunction&& fun) {
//    return boost::make_transform_iterator(
//            std::forward<Iterator>(it), 
//            std::forward<UnaryFunction>(fun));
//}
#endif



template< class InputIt, class UnaryFunction >
#ifdef __aboria_use_thrust_algorithms__
InputIt for_each( InputIt first, InputIt last, UnaryFunction f ) {
    return thrust::for_each(first,last,f);
#else
UnaryFunction for_each( InputIt first, InputIt last, UnaryFunction f ) {
    return std::for_each(first,last,f);
#endif
}

template<typename T1, typename T2>
void sort_by_key(T1 start_keys,
        T1 end_keys,
        T2 start_data) {

#ifdef __aboria_use_thrust_algorithms__
    thrust::sort_by_key(start_keys,end_keys,start_data);
#else
    typedef zip_iterator<tuple_ns::tuple<T1,T2>,mpl::vector<>> pair_zip_type;
    typedef typename pair_zip_type::reference reference;
    typedef typename pair_zip_type::value_type value_type;

    std::sort(
            pair_zip_type(start_keys,start_data),
            pair_zip_type(end_keys,start_data+std::distance(start_keys,end_keys)),
            detail::iter_comp<value_type>());
#endif
}

template<typename ForwardIterator, typename InputIterator, typename OutputIterator>
void lower_bound(
        ForwardIterator first,
        ForwardIterator last,
        InputIterator values_first,
        InputIterator values_last,
        OutputIterator result) {

#ifdef __aboria_use_thrust_algorithms__
    thrust::lower_bound(first,last,values_first,values_last,result);
#else
    std::transform(values_first,values_last,result,
            detail::lower_bound_impl<ForwardIterator>(first,last));
#endif
}

template<typename ForwardIterator, typename InputIterator, typename OutputIterator>
void upper_bound(
        ForwardIterator first,
        ForwardIterator last,
        InputIterator values_first,
        InputIterator values_last,
        OutputIterator result) {

#ifdef __aboria_use_thrust_algorithms__
    thrust::upper_bound(first,last,values_first,values_last,result);
#else
    std::transform(values_first,values_last,result,
            detail::upper_bound_impl<ForwardIterator>(first,last));
#endif
}

template< class InputIt, class T, class BinaryOperation >
T reduce( 
    InputIt first, 
    InputIt last, T init,
    BinaryOperation op ) {

#ifdef __aboria_use_thrust_algorithms__
    return thrust::reduce(first,last,init,op);
#else
    return std::accumulate(first,last,init,op);
#endif
}

template <class InputIterator, class OutputIterator, class UnaryOperation>
OutputIterator transform (
        InputIterator first, InputIterator last,
        OutputIterator result, UnaryOperation op) {


#ifdef __aboria_use_thrust_algorithms__
    return thrust::transform(first,last,result,op);
#else
    return std::transform(first,last,result,op);
#endif
}


template <class ForwardIterator>
void sequence (ForwardIterator first, ForwardIterator last) {


#ifdef __aboria_use_thrust_algorithms__
    thrust::sequence(first,last);
#else
    counting_iterator<unsigned int> count(0);
    std::transform(first,last,count,first,
        [](const typename std::iterator_traits<ForwardIterator>::reference, const unsigned int i) {
            return i;
        });
#endif
}

template<typename ForwardIterator , typename UnaryOperation >
void tabulate (
        ForwardIterator first,
        ForwardIterator last,
        UnaryOperation  unary_op) {	


#ifdef __aboria_use_thrust_algorithms__
    thrust::tabulate(first,last,unary_op);
#else
    counting_iterator<unsigned int> count(0);
    std::transform(first,last,count,first,
        [&unary_op](const typename std::iterator_traits<ForwardIterator>::reference, const unsigned int i) {
            return unary_op(i);
        });
#endif

}

template<typename InputIterator , typename OutputIterator >
OutputIterator copy(InputIterator first, InputIterator last, OutputIterator result) {
#ifdef __aboria_use_thrust_algorithms__
    return thrust::copy(first,last,result);
#else
    return std::copy(first,last,result);
#endif
}

template<typename InputIterator, typename OutputIterator, typename UnaryFunction, 
    typename T, typename AssociativeOperator>
OutputIterator transform_exclusive_scan(
    InputIterator first, InputIterator last,
    OutputIterator result,
    UnaryFunction unary_op, T init, AssociativeOperator binary_op) {

#ifdef __aboria_use_thrust_algorithms__
    return thrust::transform_exclusive_scan(first,last,result,unary_op,init,binary_op);
#else
    const size_t n = last-first;
    result[0] = init;
    for (int i=1; i<n; ++i) {
        result[i] = binary_op(result[i-1],unary_op(first[i-1]));
    }
    return result + n;
#endif
}


template<typename InputIterator1, typename InputIterator2, 
    typename InputIterator3, typename RandomAccessIterator , typename Predicate >
void scatter_if(
        InputIterator1 first, InputIterator1 last,
        InputIterator2 map, InputIterator3 stencil,
        RandomAccessIterator output, Predicate pred) {

#ifdef __aboria_use_thrust_algorithms__
    scatter_if(first,last,map,stencil,output,pred);
#else
    const size_t n = last-first;
    for (int i=0; i<n; ++i) {
        if (pred(stencil[i])) {
            output[map[i]] = first[i];
        }
    }
#endif
}

template<typename InputIterator1, typename InputIterator2, 
    typename OutputIterator, typename Predicate>
OutputIterator copy_if(
        InputIterator1 first, InputIterator1 last, 
        InputIterator2 stencil, OutputIterator result, Predicate pred) {

#ifdef __aboria_use_thrust_algorithms__
    return copy_if(first,last,stencil,result,pred);
#else
    const size_t n = last-first;
    for (int i=0; i<n; ++i) {
        if (pred(stencil[i])) {
            *result = first[i];
            ++result;
        }
    }
    return result;
#endif
}

}
}

#endif
