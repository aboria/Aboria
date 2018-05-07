
#ifndef ALGORITHMS_H_
#define ALGORITHMS_H_

#include "../CudaInclude.h"
#include "Get.h"
#include "Traits.h"
#include <algorithm>
#include <omp.h>
#include <random>

namespace Aboria {

namespace detail {

// TODO: this is a complicated mass of preprocessor, can I improve?
#ifdef HAVE_THRUST
template <typename Traits>
size_t concurrent_processes(
    typename std::enable_if<
        std::is_same<typename Traits::vector_unsigned_int,
                     thrust::device_vector<unsigned int>>::value>::type * = 0) {
  // if using GPU just return "lots"
#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
  return 9999;
#else
  return omp_get_max_threads();
#endif
}

template <typename Traits>
size_t concurrent_processes(
    typename std::enable_if<
        !std::is_same<typename Traits::vector_unsigned_int,
                      thrust::device_vector<unsigned int>>::value>::type * =
        0) {
#ifdef HAVE_OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

#else // no thrust

template <typename Traits> size_t concurrent_processes() {
#ifdef HAVE_OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

#endif

template <typename T> struct is_std_iterator {
#ifdef HAVE_THRUST
  typedef std::integral_constant<bool, false>
      /*
      !std::is_convertible<typename
      std::iterator_traits<T>::iterator_category,
                           thrust::random_access_device_iterator_tag>::value
      std::is_same<typename std::iterator_traits<T>::iterator_category,
                   std::random_access_iterator_tag>::value ||
          std::is_convertible<
              typename std::iterator_traits<T>::iterator_category,
              boost::random_access_traversal_tag>::value>
                         */

      type;
  // typedef std::integral_constant<bool,false> type;
#else
  typedef std::integral_constant<bool, true> type;
#endif
};

template <typename T> struct lower_bound_impl {
  const T &values_first, values_last;
  lower_bound_impl(const T &values_first, const T &values_last)
      : values_first(values_first), values_last(values_last) {}
  typedef typename std::iterator_traits<T>::value_type value_type;
  typedef typename std::iterator_traits<T>::difference_type difference_type;
  difference_type operator()(const value_type &search) {
    return std::distance(values_first,
                         std::lower_bound(values_first, values_last, search));
  }
};

template <typename T> struct upper_bound_impl {
  const T &values_first, values_last;
  upper_bound_impl(const T &values_first, const T &values_last)
      : values_first(values_first), values_last(values_last) {}
  typedef typename std::iterator_traits<T>::value_type value_type;
  typedef typename std::iterator_traits<T>::difference_type difference_type;
  difference_type operator()(const value_type &search) {
    return std::distance(values_first,
                         std::upper_bound(values_first, values_last, search));
  }
};

template <class ForwardIt, class T>
void fill(ForwardIt first, ForwardIt last, const T &value, std::true_type) {
  std::fill(first, last, value);
}

#ifdef HAVE_THRUST
template <class ForwardIt, class T>
void fill(ForwardIt first, ForwardIt last, const T &value, std::false_type) {
  thrust::fill(first, last, value);
}
#endif

template <class ForwardIt, class T>
void fill(ForwardIt first, ForwardIt last, const T &value) {
  fill(first, last, value, typename is_std_iterator<ForwardIt>::type());
}

template <class InputIt, class UnaryFunction>
UnaryFunction for_each(InputIt first, InputIt last, UnaryFunction f,
                       std::true_type) {
  return std::for_each(first, last, f);
}

#ifdef HAVE_THRUST
template <class InputIt, class UnaryFunction>
UnaryFunction for_each(InputIt first, InputIt last, UnaryFunction f,
                       std::false_type) {
  thrust::for_each(first, last, f);
  return f;
}
#endif

template <class InputIt, class UnaryFunction>
UnaryFunction for_each(InputIt first, InputIt last, UnaryFunction f) {
  return for_each(first, last, f, typename is_std_iterator<InputIt>::type());
}

template <typename RandomIt>
void sort(RandomIt start, RandomIt end, std::true_type) {
  std::sort(start, end);
}

#ifdef HAVE_THRUST
template <typename RandomIt>
void sort(RandomIt start, RandomIt end, std::false_type) {
  thrust::sort(start, end);
}
#endif

template <typename RandomIt> void sort(RandomIt start, RandomIt end) {
  detail::sort(start, end, typename is_std_iterator<RandomIt>::type());
}

template <typename RandomIt, typename StrictWeakOrdering>
void sort(RandomIt start, RandomIt end, StrictWeakOrdering comp,
          std::true_type) {
  std::sort(start, end, comp);
}

#ifdef HAVE_THRUST
template <typename RandomIt, typename StrictWeakOrdering>
void sort(RandomIt start, RandomIt end, StrictWeakOrdering comp,
          std::false_type) {
  thrust::sort(start, end, comp);
}
#endif

template <typename RandomIt, typename StrictWeakOrdering>
void sort(RandomIt start, RandomIt end, StrictWeakOrdering comp) {
  detail::sort(start, end, comp, typename is_std_iterator<RandomIt>::type());
}

template <typename T1, typename T2>
void sort_by_key(T1 start_keys, T1 end_keys, T2 start_data, std::true_type) {
  typedef zip_iterator<std::tuple<T1, T2>, mpl::vector<>> pair_zip_type;

  std::sort(
      pair_zip_type(start_keys, start_data),
      pair_zip_type(end_keys, start_data + std::distance(start_keys, end_keys)),
      [](auto a, auto b) {
        return std::get<0>(a.get_tuple()) < std::get<0>(b.get_tuple());
      });
  /*
std::sort(
  boost::make_zip_iterator(boost::make_tuple(start_keys, start_data)),
  boost::make_zip_iterator(boost::make_tuple(
      end_keys, start_data + std::distance(start_keys, end_keys))),
  [](auto a, auto b) { return a.template get<0>() < b.template get<0>(); });
  */
}

#ifdef HAVE_THRUST
template <typename T1, typename T2>
void sort_by_key(T1 start_keys, T1 end_keys, T2 start_data, std::false_type) {
  thrust::sort_by_key(start_keys, end_keys, start_data);
}
#endif

// TODO: only works for random access iterators
template <typename T1, typename T2>
void sort_by_key(T1 start_keys, T1 end_keys, T2 start_data) {
  // TODO: how to check its generically a std iterator as opposed to a
  // thrust::iterator
  detail::sort_by_key(start_keys, end_keys, start_data,
                      typename is_std_iterator<T1>::type());
}

template <typename ForwardIterator, typename InputIterator,
          typename OutputIterator>
void lower_bound(ForwardIterator first, ForwardIterator last,
                 InputIterator values_first, InputIterator values_last,
                 OutputIterator result, std::true_type) {
  std::transform(values_first, values_last, result,
                 detail::lower_bound_impl<ForwardIterator>(first, last));
}

#ifdef HAVE_THRUST
template <typename ForwardIterator, typename InputIterator,
          typename OutputIterator>
void lower_bound(ForwardIterator first, ForwardIterator last,
                 InputIterator values_first, InputIterator values_last,
                 OutputIterator result, std::false_type) {
  thrust::lower_bound(first, last, values_first, values_last, result);
}
#endif

template <typename ForwardIterator, typename InputIterator,
          typename OutputIterator>
void lower_bound(ForwardIterator first, ForwardIterator last,
                 InputIterator values_first, InputIterator values_last,
                 OutputIterator result) {

  detail::lower_bound(first, last, values_first, values_last, result,
                      typename is_std_iterator<ForwardIterator>::type());
}

template <typename ForwardIterator, typename LessThanComparable>
ForwardIterator lower_bound(ForwardIterator first, ForwardIterator last,
                            const LessThanComparable &value, std::true_type) {

  return std::lower_bound(first, last, value);
}

#ifdef HAVE_THRUST
template <typename ForwardIterator, typename LessThanComparable>
ForwardIterator lower_bound(ForwardIterator first, ForwardIterator last,
                            const LessThanComparable &value, std::false_type) {

  return thrust::lower_bound(first, last, value);
}
#endif

template <typename ForwardIterator, typename LessThanComparable>
ForwardIterator lower_bound(ForwardIterator first, ForwardIterator last,
                            const LessThanComparable &value) {

  return detail::lower_bound(first, last, value,
                             typename is_std_iterator<ForwardIterator>::type());
}

template <typename ForwardIterator, typename InputIterator,
          typename OutputIterator>
void upper_bound(ForwardIterator first, ForwardIterator last,
                 InputIterator values_first, InputIterator values_last,
                 OutputIterator result, std::true_type) {

  std::transform(values_first, values_last, result,
                 detail::upper_bound_impl<ForwardIterator>(first, last));
}

#ifdef HAVE_THRUST
template <typename ForwardIterator, typename InputIterator,
          typename OutputIterator>
void upper_bound(ForwardIterator first, ForwardIterator last,
                 InputIterator values_first, InputIterator values_last,
                 OutputIterator result, std::false_type) {
  thrust::upper_bound(first, last, values_first, values_last, result);
}
#endif

template <typename ForwardIterator, typename InputIterator,
          typename OutputIterator>
void upper_bound(ForwardIterator first, ForwardIterator last,
                 InputIterator values_first, InputIterator values_last,
                 OutputIterator result) {

  detail::upper_bound(first, last, values_first, values_last, result,
                      typename is_std_iterator<ForwardIterator>::type());
}

template <class InputIt, class T, class BinaryOperation>
T reduce(InputIt first, InputIt last, T init, BinaryOperation op,
         std::true_type) {

  return std::accumulate(first, last, init, op);
}

#ifdef HAVE_THRUST
template <class InputIt, class T, class BinaryOperation>
T reduce(InputIt first, InputIt last, T init, BinaryOperation op,
         std::false_type) {

  return thrust::reduce(first, last, init, op);
}
#endif

template <class InputIt, class T, class BinaryOperation>
T reduce(InputIt first, InputIt last, T init, BinaryOperation op) {

  detail::reduce(first, last, init, op,
                 typename is_std_iterator<InputIt>::type());
}

template <class InputIterator, class OutputIterator, class UnaryOperation>
OutputIterator transform(InputIterator first, InputIterator last,
                         OutputIterator result, UnaryOperation op,
                         std::true_type) {
  return std::transform(first, last, result, op);
}

#ifdef HAVE_THRUST
template <class InputIterator, class OutputIterator, class UnaryOperation>
OutputIterator transform(InputIterator first, InputIterator last,
                         OutputIterator result, UnaryOperation op,
                         std::false_type) {
  return thrust::transform(first, last, result, op);
}
#endif

template <class InputIterator, class OutputIterator, class UnaryOperation>
OutputIterator transform(InputIterator first, InputIterator last,
                         OutputIterator result, UnaryOperation op) {
  return detail::transform(first, last, result, op,
                           typename is_std_iterator<OutputIterator>::type());
}

template <class ForwardIterator>
void sequence(ForwardIterator first, ForwardIterator last, std::true_type) {
  boost::counting_iterator<unsigned int> count(0);
  std::transform(
      first, last, count, first,
      [](const typename std::iterator_traits<ForwardIterator>::reference,
         const unsigned int i) { return i; });
}

template <class ForwardIterator, typename T>
void sequence(ForwardIterator first, ForwardIterator last, T init,
              std::true_type) {
  boost::counting_iterator<unsigned int> count(init);
  std::transform(
      first, last, count, first,
      [](const typename std::iterator_traits<ForwardIterator>::reference,
         const unsigned int i) { return i; });
}

#ifdef HAVE_THRUST
template <class ForwardIterator>
void sequence(ForwardIterator first, ForwardIterator last, std::false_type) {
  thrust::sequence(first, last);
}

template <class ForwardIterator, typename T>
void sequence(ForwardIterator first, ForwardIterator last, T init,
              std::false_type) {
  thrust::sequence(first, last, init);
}
#endif

template <class ForwardIterator>
void sequence(ForwardIterator first, ForwardIterator last) {
  detail::sequence(first, last,
                   typename is_std_iterator<ForwardIterator>::type());
}

template <class ForwardIterator, typename T>
void sequence(ForwardIterator first, ForwardIterator last, T init) {
  detail::sequence(first, last, init,
                   typename is_std_iterator<ForwardIterator>::type());
}

template <typename ForwardIterator, typename UnaryOperation>
void tabulate(ForwardIterator first, ForwardIterator last,
              UnaryOperation unary_op, std::true_type) {

  boost::counting_iterator<unsigned int> count(0);
  std::transform(
      first, last, count, first,
      [&unary_op](
          const typename std::iterator_traits<ForwardIterator>::reference,
          const unsigned int i) { return unary_op(i); });
}

#ifdef HAVE_THRUST
template <typename ForwardIterator, typename UnaryOperation>
void tabulate(ForwardIterator first, ForwardIterator last,
              UnaryOperation unary_op, std::false_type) {
  thrust::tabulate(first, last, unary_op);
}
#endif

template <typename ForwardIterator, typename UnaryOperation>
void tabulate(ForwardIterator first, ForwardIterator last,
              UnaryOperation unary_op) {
  detail::tabulate(first, last, unary_op,
                   typename is_std_iterator<ForwardIterator>::type());
}

template <class ForwardIt, class UnaryPredicate>
ForwardIt partition(ForwardIt first, ForwardIt last, UnaryPredicate p,
                    std::true_type) {
  return std::partition(first, last, p);
}

#ifdef HAVE_THRUST
template <class ForwardIt, class UnaryPredicate>
ForwardIt partition(ForwardIt first, ForwardIt last, UnaryPredicate p,
                    std::false_type) {
  return thrust::partition(first, last, p);
}
#endif

template <class ForwardIt, class UnaryPredicate>
ForwardIt partition(ForwardIt first, ForwardIt last, UnaryPredicate p) {
  return detail::partition(first, last, p,
                           typename is_std_iterator<ForwardIt>::type());
}

template <class ForwardIt, class UnaryPredicate>
ForwardIt stable_partition(ForwardIt first, ForwardIt last, UnaryPredicate p,
                           std::true_type) {
  return std::stable_partition(first, last, p);
}

#ifdef HAVE_THRUST
template <class ForwardIt, class UnaryPredicate>
ForwardIt stable_partition(ForwardIt first, ForwardIt last, UnaryPredicate p,
                           std::false_type) {
  return thrust::stable_partition(first, last, p);
}
#endif

template <class ForwardIt, class UnaryPredicate>
ForwardIt stable_partition(ForwardIt first, ForwardIt last, UnaryPredicate p) {
  return detail::stable_partition(first, last, p,
                                  typename is_std_iterator<ForwardIt>::type());
}

template <class ForwardIt>
ForwardIt unique(ForwardIt first, ForwardIt last, std::true_type) {
  return std::unique(first, last);
}

#ifdef HAVE_THRUST
template <class ForwardIt>
ForwardIt unique(ForwardIt first, ForwardIt last, std::false_type) {
  return thrust::unique(first, last);
}
#endif

template <class ForwardIt> ForwardIt unique(ForwardIt first, ForwardIt last) {
  return detail::unique(first, last,
                        typename is_std_iterator<ForwardIt>::type());
}

template <typename InputIterator, typename OutputIterator>
OutputIterator copy(InputIterator first, InputIterator last,
                    OutputIterator result, std::true_type) {
  return std::copy(first, last, result);
}

#ifdef HAVE_THRUST
template <typename InputIterator, typename OutputIterator>
OutputIterator copy(InputIterator first, InputIterator last,
                    OutputIterator result, std::false_type) {
  return thrust::copy(first, last, result);
}
#endif

template <typename InputIterator, typename OutputIterator>
OutputIterator copy(InputIterator first, InputIterator last,
                    OutputIterator result) {
  return detail::copy(first, last, result,
                      typename is_std_iterator<InputIterator>::type());
}

template <typename InputIterator, typename OutputIterator,
          typename UnaryFunction, typename T, typename AssociativeOperator>
OutputIterator
transform_exclusive_scan(InputIterator first, InputIterator last,
                         OutputIterator result, UnaryFunction unary_op, T init,
                         AssociativeOperator binary_op, std::true_type) {

  const size_t n = last - first;
  result[0] = init;
  for (size_t i = 1; i < n; ++i) {
    result[i] = binary_op(result[i - 1], unary_op(first[i - 1]));
  }
  return result + n;
}

#ifdef HAVE_THRUST
template <typename InputIterator, typename OutputIterator,
          typename UnaryFunction, typename T, typename AssociativeOperator>
OutputIterator
transform_exclusive_scan(InputIterator first, InputIterator last,
                         OutputIterator result, UnaryFunction unary_op, T init,
                         AssociativeOperator binary_op, std::false_type) {

  return thrust::transform_exclusive_scan(first, last, result, unary_op, init,
                                          binary_op);
}
#endif

template <typename InputIterator, typename OutputIterator,
          typename UnaryFunction, typename T, typename AssociativeOperator>
OutputIterator transform_exclusive_scan(InputIterator first, InputIterator last,
                                        OutputIterator result,
                                        UnaryFunction unary_op, T init,
                                        AssociativeOperator binary_op) {

  return detail::transform_exclusive_scan(
      first, last, result, unary_op, init, binary_op,
      typename is_std_iterator<InputIterator>::type());
}

template <class InputIt, class OutputIt>
OutputIt inclusive_scan(InputIt first, InputIt last, OutputIt d_first,
                        std::true_type) {
#if __cplusplus >= 201703L
  // C++17 code here
  return std::inclusive_scan(first, last, d_first);
#else
  return std::partial_sum(first, last, d_first);
#endif
}

#ifdef HAVE_THRUST
template <class InputIt, class OutputIt>
OutputIt inclusive_scan(InputIt first, InputIt last, OutputIt d_first,
                        std::false_type) {
  return thrust::inclusive_scan(first, last, d_first);
}
#endif

template <class InputIt, class OutputIt>
OutputIt inclusive_scan(InputIt first, InputIt last, OutputIt d_first) {
  return detail::inclusive_scan(first, last, d_first,
                                typename is_std_iterator<InputIt>::type());
}

template <class InputIt, class OutputIt, class T>
OutputIt exclusive_scan(InputIt first, InputIt last, OutputIt d_first, T init,
                        std::true_type) {
#if __cplusplus >= 201703L
  // C++17 code here
  return std::exclusive_scan(first, last, d_first,init);
#else
  if (first != last) {
    *d_first++ = init;
    for (; first != last - 1; ++first, ++d_first) {
      *d_first = *(d_first - 1) + *first;
    }
  }
  return d_first;
#endif
}

#ifdef HAVE_THRUST
template <class InputIt, class OutputIt, class T>
OutputIt exclusive_scan(InputIt first, InputIt last, OutputIt d_first, T init,
                        std::false_type) {
  return thrust::exclusive_scan(first, last, d_first,init);
}
#endif

template <class InputIt, class OutputIt, class T>
OutputIt exclusive_scan(InputIt first, InputIt last, OutputIt d_first, T init) {
  return detail::exclusive_scan(first, last, d_first, init,
                                typename is_std_iterator<InputIt>::type());
}

template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator, typename T>
OutputIterator
exclusive_scan_by_key(InputIterator1 first1, InputIterator1 last1,
                      InputIterator2 first2, OutputIterator result, T init,
                      std::true_type) {
  for (; first1 != last1; ++first2) {
    *result++ = init;
    ++first1;
    for (; *first1 == *(first1 - 1) && first1 != last1;
         ++first1, ++first2, ++result)
      *result = *(result - 1) + *first2;
  }
  return result;
}

#ifdef HAVE_THRUST
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator, typename T>
OutputIterator
exclusive_scan_by_key(InputIterator1 first1, InputIterator1 last1,
                      InputIterator2 first2, OutputIterator result, T init,
                      std::false_type) {
  return thrust::exclusive_scan_by_key(first1, last1, first2, result, init);
}
#endif

template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator, typename T>
OutputIterator
exclusive_scan_by_key(InputIterator1 first1, InputIterator1 last1,
                      InputIterator2 first2, OutputIterator result, T init) {
  return detail::exclusive_scan_by_key(
      first1, last1, first2, result, init,
      typename is_std_iterator<InputIterator1>::type());
}

template <typename InputIterator1, typename InputIterator2,
          typename RandomAccessIterator>
void scatter(InputIterator1 first, InputIterator1 last, InputIterator2 map,
             RandomAccessIterator output, std::true_type) {
  const size_t n = last - first;
  for (size_t i = 0; i < n; ++i) {
    output[map[i]] = first[i];
  }
}

#ifdef HAVE_THRUST
template <typename InputIterator1, typename InputIterator2,
          typename RandomAccessIterator>
void scatter(InputIterator1 first, InputIterator1 last, InputIterator2 map,
             RandomAccessIterator output, std::false_type) {
  thrust::scatter(first, last, map, output);
}
#endif

template <typename InputIterator1, typename InputIterator2,
          typename RandomAccessIterator>
void scatter(InputIterator1 first, InputIterator1 last, InputIterator2 map,
             RandomAccessIterator output) {
  detail::scatter(first, last, map, output,
                  typename is_std_iterator<RandomAccessIterator>::type());
}

template <typename InputIterator1, typename InputIterator2,
          typename InputIterator3, typename RandomAccessIterator,
          typename Predicate>
void scatter_if(InputIterator1 first, InputIterator1 last, InputIterator2 map,
                InputIterator3 stencil, RandomAccessIterator output,
                Predicate pred, std::true_type) {

  const size_t n = last - first;
  for (size_t i = 0; i < n; ++i) {
    if (pred(stencil[i])) {
      output[map[i]] = first[i];
    }
  }
}

#ifdef HAVE_THRUST
template <typename InputIterator1, typename InputIterator2,
          typename InputIterator3, typename RandomAccessIterator,
          typename Predicate>
void scatter_if(InputIterator1 first, InputIterator1 last, InputIterator2 map,
                InputIterator3 stencil, RandomAccessIterator output,
                Predicate pred, std::false_type) {

  thrust::scatter_if(first, last, map, stencil, output, pred);
}
#endif

template <typename InputIterator1, typename InputIterator2,
          typename InputIterator3, typename RandomAccessIterator,
          typename Predicate>
void scatter_if(InputIterator1 first, InputIterator1 last, InputIterator2 map,
                InputIterator3 stencil, RandomAccessIterator output,
                Predicate pred) {

  detail::scatter_if(first, last, map, stencil, output, pred,
                     typename is_std_iterator<RandomAccessIterator>::type());
}

template <typename InputIterator1, typename InputIterator2,
          typename InputIterator3, typename RandomAccessIterator>
void scatter_if(InputIterator1 first, InputIterator1 last, InputIterator2 map,
                InputIterator3 stencil, RandomAccessIterator output,
                std::true_type) {

  const size_t n = last - first;
  for (size_t i = 0; i < n; ++i) {
    if (stencil[i]) {
      output[map[i]] = first[i];
    }
  }
}

#ifdef HAVE_THRUST
template <typename InputIterator1, typename InputIterator2,
          typename InputIterator3, typename RandomAccessIterator>
void scatter_if(InputIterator1 first, InputIterator1 last, InputIterator2 map,
                InputIterator3 stencil, RandomAccessIterator output,
                std::false_type) {

  thrust::scatter_if(first, last, map, stencil, output);
}
#endif

template <typename InputIterator1, typename InputIterator2,
          typename InputIterator3, typename RandomAccessIterator>
void scatter_if(InputIterator1 first, InputIterator1 last, InputIterator2 map,
                InputIterator3 stencil, RandomAccessIterator output) {

  detail::scatter_if(first, last, map, stencil, output,
                     typename is_std_iterator<RandomAccessIterator>::type());
}

template <typename InputIterator, typename RandomAccessIterator,
          typename OutputIterator>
void gather(InputIterator map_first, InputIterator map_last,
            RandomAccessIterator input_first, OutputIterator result,
            std::true_type) {
  std::transform(map_first, map_last, result,
                 [&input_first](typename InputIterator::value_type const &i) {
                   return input_first[i];
                 });
}

#ifdef HAVE_THRUST
template <typename InputIterator, typename RandomAccessIterator,
          typename OutputIterator>
void gather(InputIterator map_first, InputIterator map_last,
            RandomAccessIterator input_first, OutputIterator result,
            std::false_type) {
  thrust::gather(map_first, map_last, input_first, result);
}
#endif

template <typename InputIterator, typename RandomAccessIterator,
          typename OutputIterator>
void gather(InputIterator map_first, InputIterator map_last,
            RandomAccessIterator input_first, OutputIterator result) {
  detail::gather(map_first, map_last, input_first, result,
                 typename is_std_iterator<RandomAccessIterator>::type());
}

template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator, typename Predicate>
OutputIterator copy_if(InputIterator1 first, InputIterator1 last,
                       InputIterator2 stencil, OutputIterator result,
                       Predicate pred, std::true_type) {
  // TODO: need to parallise this....
  const size_t n = last - first;
  for (size_t i = 0; i < n; ++i) {
    if (pred(stencil[i])) {
      *result = first[i];
      ++result;
    }
  }
  return result;
}

#ifdef HAVE_THRUST
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator, typename Predicate>
OutputIterator copy_if(InputIterator1 first, InputIterator1 last,
                       InputIterator2 stencil, OutputIterator result,
                       Predicate pred, std::false_type) {

  return thrust::copy_if(first, last, stencil, result, pred);
}
#endif

template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator, typename Predicate>
OutputIterator copy_if(InputIterator1 first, InputIterator1 last,
                       InputIterator2 stencil, OutputIterator result,
                       Predicate pred) {

  return detail::copy_if(first, last, stencil, result, pred,
                         typename is_std_iterator<InputIterator1>::type());
}

template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator>
OutputIterator copy_if(InputIterator1 first, InputIterator1 last,
                       InputIterator2 stencil, OutputIterator result) {

  return detail::copy_if(first, last, stencil, result,
                         [](auto i) { return static_cast<bool>(i); },
                         typename is_std_iterator<InputIterator1>::type());
}

} // namespace detail
} // namespace Aboria

#endif
