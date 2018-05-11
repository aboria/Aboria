
#ifndef GET_H_
#define GET_H_

#include <algorithm>
#include <boost/mpl/find.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/vector.hpp>
#include <tuple>
#include <type_traits>

#include "detail/GetterType.h"

namespace Aboria {

///
/// @brief Generic type that can be used by ::get<variable>()
///
/// @tparam Tuple a `std::tuple` or `thrust::tuple` of values
/// @tparam MplVector a boost::mpl typelist of variable types
///
template <typename Tuple, typename MplVector> struct getter_type;

///
/// @brief Generic iterator that zips multiple iterators together
///
/// @tparam iterator_tuple_type `std::tuple` or `thrust::tuple` of individual
///         iterators
/// @tparam mpl_vector_type a `boost::mpl` typelist of variable types
///
template <typename iterator_tuple_type, typename mpl_vector_type>
class zip_iterator;

///
/// @brief Generic pointer that zips multiple pointers together
///
/// @tparam TupleType std::tuple or thrust::tuple of individual
///         iterators
/// @tparam MplVector a boost::mpl typelist
///
template <typename TupleType, typename MplVector> struct zip_pointer;

} // namespace Aboria

#include "Log.h"
#include "detail/Get.h"

namespace Aboria {

template <typename Tuple, typename MplVector> struct getter_type {};

///
/// @brief specialisation of @ref getter_type for `std::tuple`
///
/// This needs to be all host only code, use `thrust::tuple` specialisation
/// for device code
///
/// @TODO: api is not good, consider doing an iterator_facade type thing
///
template <typename MplVector, typename... Types>
struct getter_type<std::tuple<Types...>, MplVector> {
  typedef std::tuple<Types...> tuple_type;
  typedef MplVector mpl_vector_type;

  typedef typename detail::getter_helper<tuple_type>::tuple_reference
      tuple_reference;
  template <typename T>
  using elem_by_type = detail::get_elem_by_type<T, mpl_vector_type>;
  template <typename T>
  using return_type = typename detail::getter_helper<
      tuple_type>::template return_type<elem_by_type<T>::index>;

  getter_type() {}

  explicit getter_type(const tuple_type &data) : data(data) {}

  getter_type(const getter_type &other) : data(other.data) {}

  getter_type(getter_type &&other) : data(std::move(other.data)) {}

  template <typename T = tuple_reference,
            typename = typename std::enable_if<
                !std::is_same<T, tuple_type>::value>::type>
  getter_type(const getter_type<tuple_reference, mpl_vector_type> &other)
      : data(other.data) {}

  template <typename T = tuple_reference,
            typename = typename std::enable_if<
                !std::is_same<T, tuple_type>::value>::type>
  getter_type(getter_type<tuple_reference, mpl_vector_type> &&other)
      : data(std::move(other.data)) {}

  template <typename T, typename = typename std::enable_if<
                            (!std::is_same<T, tuple_type>::value) &&
                            (!std::is_same<T, tuple_reference>::value)>::type>
  getter_type(const getter_type<T, mpl_vector_type> &other)
      : data(other.data) {}

  template <typename T, typename = typename std::enable_if<
                            (!std::is_same<T, tuple_type>::value) &&
                            (!std::is_same<T, tuple_reference>::value)>::type>
  getter_type(getter_type<T, mpl_vector_type> &&other)
      : data(std::move(other.data)) {}

  template <typename T1, typename T2, typename... T3>
  getter_type(T1 &&arg1, T2 &&arg2, T3 &&... args)
      : data(std::forward<T1>(arg1), std::forward<T2>(arg2),
             std::forward<T3>(args)...) {}

  template <typename T> int throw_away(const T &in) { return 0; }

  template <typename Tuple, std::size_t... I>
  void copy_impl(const Tuple &other_data, detail::index_sequence<I...>) {
    int dummy[] = {0,
                   throw_away(std::get<I>(data) = std::get<I>(other_data))...};
    static_cast<void>(dummy);
  }

  getter_type &operator=(const getter_type &other) {
    // copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
    data = other.data;
    return *this;
  }
  getter_type &operator=(getter_type &&other) {
    // copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
    data = other.data;
    return *this;
  }

  template <typename T = tuple_reference,
            typename = typename std::enable_if<
                !std::is_same<T, tuple_type>::value>::type>
  getter_type &
  operator=(const getter_type<tuple_reference, mpl_vector_type> &other) {
    // copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
    data = other.data;
    return *this;
  }

  template <typename T = tuple_reference,
            typename = typename std::enable_if<
                !std::is_same<T, tuple_type>::value>::type>
  getter_type &
  operator=(getter_type<tuple_reference, mpl_vector_type> &&other) {
    // copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
    data = other.data;
    return *this;
  }

  template <typename T, typename = typename std::enable_if<
                            (!std::is_same<T, tuple_reference>::value) &&
                            (!std::is_same<T, tuple_type>::value)>::type>
  getter_type &operator=(const getter_type<T, mpl_vector_type> &other) {
    // copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
    data = other.data;
    return *this;
  }

  template <typename T, typename = typename std::enable_if<
                            (!std::is_same<T, tuple_reference>::value) &&
                            (!std::is_same<T, tuple_type>::value)>::type>
  getter_type &operator=(getter_type<T, mpl_vector_type> &&other) {
    // copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
    data = other.data;
    return *this;
  }

#if defined(__aboria_use_thrust_algorithms__) || defined(__CUDACC__)
  template <typename T1, typename T2, typename PointerType,
            typename DerivedType>
  getter_type &
  operator=(const thrust::reference<getter_type<T1, T2>, PointerType,
                                    DerivedType> &other) {
    data = static_cast<getter_type<T1, T2>>(other).data;
    return *this;
  }
#endif

  template <typename T1, typename T2>
  bool operator==(const getter_type<T1, T2> &other) {
    return data == other.data;
  }

  void swap(getter_type &other) { data.swap(other.data); }

  template <typename tuple_type2, std::size_t... I>
  void swap_via_tie(tuple_type2 &tuple, detail::index_sequence<I...>) {
    tuple_type tmp = std::tie(std::get<I>(tuple)...);
    data.swap(tmp);
  }

  CUDA_HOST_DEVICE
  const tuple_type &get_tuple() const {
#if defined(__CUDA_ARCH__)
    ERROR_CUDA("Cannot use `getter_type<std::tuple>` in device code");
#endif
    return data;
  }

  CUDA_HOST_DEVICE
  tuple_type &get_tuple() {
#if defined(__CUDA_ARCH__)
    ERROR_CUDA("Cannot use `getter_type<std::tuple>` in device code");
#endif
    return data;
  }

  tuple_type data;
};

#ifdef HAVE_THRUST
///
/// @brief specialisation of @ref getter_type for `thrust::tuple`
///
/// Use device code only in this specialisation
/// @TODO: api is not good, consider doing an iterator_facade type thing
///
template <typename MplVector, typename TT1, typename TT2, typename TT3,
          typename TT4, typename TT5, typename TT6, typename TT7, typename TT8,
          typename TT9>
struct getter_type<thrust::tuple<TT1, TT2, TT3, TT4, TT5, TT6, TT7, TT8, TT9>,
                   MplVector> {
  typedef thrust::tuple<TT1, TT2, TT3, TT4, TT5, TT6, TT7, TT8, TT9> tuple_type;
  typedef MplVector mpl_vector_type;

  typedef typename detail::getter_helper<tuple_type>::tuple_reference
      tuple_reference;
  typedef typename detail::getter_helper<tuple_type>::tuple_device_reference
      tuple_device_reference;
  typedef typename detail::getter_helper<tuple_type>::index_type index_type;
  template <typename T>
  using elem_by_type = detail::get_elem_by_type<T, mpl_vector_type>;
  template <typename T>
  using return_type = typename detail::getter_helper<
      tuple_type>::template return_type<elem_by_type<T>::index>;

  // typedef typename detail::zip_helper<tuple_type>::pointer pointer;

  CUDA_HOST_DEVICE
  getter_type() {}

  CUDA_HOST_DEVICE
  explicit getter_type(const tuple_type &data) : data(data) {}

  CUDA_HOST_DEVICE
  getter_type(const getter_type &other) : data(other.data) {}

  template <typename T = tuple_reference,
            typename = typename std::enable_if<
                !std::is_same<T, tuple_type>::value>::type>
  CUDA_HOST_DEVICE
  getter_type(const getter_type<tuple_reference, mpl_vector_type> &other)
      : data(other.data) {}

  template <typename T = tuple_reference,
            typename = typename std::enable_if<
                !std::is_same<T, tuple_type>::value>::type>
  CUDA_HOST_DEVICE
  getter_type(getter_type<tuple_reference, mpl_vector_type> &&other)
      : data(std::move(other.data)) {}

  template <typename T = tuple_device_reference,
            typename = typename std::enable_if<
                !std::is_same<T, tuple_type>::value>::type>
  CUDA_HOST_DEVICE
  getter_type(const getter_type<tuple_device_reference, mpl_vector_type> &other)
      : data(detail::getter_helper<tuple_type>::raw_reference_cast(
            other.data, index_type())) {}

  template <typename T = tuple_device_reference,
            typename = typename std::enable_if<
                !std::is_same<T, tuple_type>::value>::type>
  CUDA_HOST_DEVICE
  getter_type(getter_type<tuple_device_reference, mpl_vector_type> &&other)
      : data(detail::getter_helper<tuple_type>::raw_reference_cast(
            std::move(other.data), index_type())) {}

  template <typename T, typename = typename std::enable_if<
                            (!std::is_same<T, tuple_type>::value) &&
                            (!std::is_same<T, tuple_reference>::value)>::type>
  CUDA_HOST_DEVICE getter_type(const getter_type<T, mpl_vector_type> &other)
      : data(other.data) {}

  template <typename T, typename = typename std::enable_if<
                            (!std::is_same<T, tuple_type>::value) &&
                            (!std::is_same<T, tuple_reference>::value)>::type>
  CUDA_HOST_DEVICE getter_type(getter_type<T, mpl_vector_type> &&other)
      : data(std::move(other.data)) {}

  template <typename T1, typename T2, typename... T3>
  CUDA_HOST_DEVICE getter_type(T1 &&arg1, T2 &&arg2, T3 &&... args)
      : data(std::forward<T1>(arg1), std::forward<T2>(arg2),
             std::forward<T3>(args)...) {}

  template <typename T> CUDA_HOST_DEVICE int throw_away(const T &in) {
    return 0;
  }

  template <typename Tuple, std::size_t... I>
  CUDA_HOST_DEVICE void copy_impl(const Tuple &other_data,
                                  detail::index_sequence<I...>) {
    int dummy[] = {
        0, throw_away(thrust::get<I>(data) = thrust::get<I>(other_data))...};
    static_cast<void>(dummy);
  }

  CUDA_HOST_DEVICE
  getter_type &operator=(const getter_type &other) {
    // copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
    data = other.data;
    return *this;
  }

  CUDA_HOST_DEVICE
  getter_type &operator=(getter_type &&other) {
    // copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
    data = other.data;
    return *this;
  }

  template <typename T = tuple_reference,
            typename = typename std::enable_if<
                !std::is_same<T, tuple_type>::value>::type>
  CUDA_HOST_DEVICE getter_type &
  operator=(const getter_type<tuple_reference, mpl_vector_type> &other) {
    // copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
    data = other.data;
    return *this;
  }

  template <typename T = tuple_reference,
            typename = typename std::enable_if<
                !std::is_same<T, tuple_type>::value>::type>
  CUDA_HOST_DEVICE getter_type &
  operator=(getter_type<tuple_reference, mpl_vector_type> &&other) {
    // copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
    data = other.data;
    return *this;
  }

  template <typename T, typename = typename std::enable_if<
                            (!std::is_same<T, tuple_reference>::value) &&
                            (!std::is_same<T, tuple_type>::value)>::type>
  CUDA_HOST_DEVICE getter_type &
  operator=(const getter_type<T, mpl_vector_type> &other) {
    // copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
    data = other.data;
    return *this;
  }

  template <typename T, typename = typename std::enable_if<
                            (!std::is_same<T, tuple_reference>::value) &&
                            (!std::is_same<T, tuple_type>::value)>::type>
  CUDA_HOST_DEVICE getter_type &
  operator=(getter_type<T, mpl_vector_type> &&other) {
    // copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
    data = other.data;
    return *this;
  }

  template <typename T, typename Tag>
  CUDA_HOST_DEVICE getter_type &
  operator=(const thrust::reference<
            getter_type<T, mpl_vector_type>,
            thrust::pointer<getter_type<T, mpl_vector_type>, Tag>> &other) {
    data = getter_type<T, mpl_vector_type>(other).data;
    return *this;
  }

  template <typename T1, typename T2>
  CUDA_HOST_DEVICE bool operator==(const getter_type<T1, T2> &other) {
    return data == other.data;
  }

  CUDA_HOST_DEVICE
  void swap(getter_type &other) { data.swap(other.data); }

  template <typename tuple_type2, std::size_t... I>
  CUDA_HOST_DEVICE void swap_via_tie(tuple_type2 &tuple,
                                     detail::index_sequence<I...>) {
    tuple_type tmp = std::tie(thrust::get<I>(tuple)...);
    data.swap(tmp);
  }

  CUDA_HOST_DEVICE
  const tuple_type &get_tuple() const { return data; }

  CUDA_HOST_DEVICE
  tuple_type &get_tuple() { return data; }

  tuple_type data;
};
#endif

///
/// @brief specialisation of @ref zip_pointer for `std::tuple`
///
/// Use host only code in this specialisation
///
template <typename MplVector, typename... Types>
struct zip_pointer<std::tuple<Types *...>, MplVector> {
  typedef std::tuple<Types *...> tuple_type;
  typedef MplVector mpl_vector_type;

  template <typename T>
  using elem_by_type = detail::get_elem_by_type<T, mpl_vector_type>;
  template <typename T>
  using return_type = typename detail::getter_helper<
      tuple_type>::template return_type<elem_by_type<T>::index>;

  typedef getter_type<typename detail::zip_helper<tuple_type>::tuple_reference,
                      mpl_vector_type>
      reference;
  typedef
      typename detail::zip_helper<tuple_type>::difference_type difference_type;
  typedef typename detail::zip_helper<tuple_type>::index_type index_type;

  zip_pointer() {}

  explicit zip_pointer(const tuple_type &data) : data(data) {}

  zip_pointer(const zip_pointer &other) : data(other.data) {}

  zip_pointer(zip_pointer &&other) : data(std::move(other.data)) {}

  template <typename tuple_type2,
            typename = typename std::enable_if<
                std::is_convertible<tuple_type2, tuple_type>::value>::type>
  zip_pointer(const zip_pointer<tuple_type2, mpl_vector_type> &other)
      : data(other.data) {}

  template <typename tuple_type2,
            typename = typename std::enable_if<
                std::is_convertible<tuple_type2, tuple_type>::value>::type>
  zip_pointer(const zip_pointer<tuple_type2, mpl_vector_type> &&other)
      : data(std::move(other.data)) {}

  template <typename T1, typename T2, typename... T3>
  zip_pointer(T1 &&arg1, T2 &&arg2, T3 &&... args)
      : data(std::forward<T1>(arg1), std::forward<T2>(arg2),
             std::forward<T3>(args)...) {}

  zip_pointer &operator=(const zip_pointer &other) {
    data = other.data;
    return *this;
  }
  zip_pointer &operator=(zip_pointer &&other) {
    data = std::move(other.data);
    return *this;
  }
  template <typename T1, typename T2>
  zip_pointer &operator=(const zip_pointer<T1, T2> &other) {
    data = other.data;
    return *this;
  }
  template <typename T1, typename T2>
  zip_pointer &operator=(zip_pointer<T1, T2> &&other) {
    data = std::move(other.data);
    return *this;
  }

  void swap(zip_pointer &other) { data.swap(other.data); }

  CUDA_HOST_DEVICE
  const tuple_type &get_tuple() const {
#if defined(__CUDA_ARCH__)
    ERROR_CUDA("Cannot use `zip_pointer` with `std::tuple` in device code");
#endif
    return data;
  }

  CUDA_HOST_DEVICE
  tuple_type &get_tuple() {
#if defined(__CUDA_ARCH__)
    ERROR_CUDA("Cannot use `zip_pointer` with `std::tuple` in device code");
#endif
    return data;
  }

  bool operator==(const zip_pointer &other) const { return equal(other); }

  zip_pointer &operator++() {
    advance(1);
    return *this;
  }

  zip_pointer operator++(int) {
    zip_pointer temp(*this);
    advance(1);
    return temp; // return saved state
  }

  zip_pointer operator+(const difference_type &n) const {
    zip_pointer ret(*this);
    ret.advance(n);
    return ret;
  }

  zip_pointer operator-(const difference_type &n) const {
    zip_pointer ret(*this);
    ret.advance(-n);
    return ret;
  }

  zip_pointer &operator+=(const difference_type &n) {
    advance(n);
    return *this;
  }

  difference_type operator-(const zip_pointer &other) const {
    return other.distance_to(*this);
  }

  bool operator<(const zip_pointer &other) const {
    return distance_to(other) > 0;
  }

  bool operator<=(const zip_pointer &other) const {
    return distance_to(other) >= 0;
  }

  bool operator>(const zip_pointer &other) const {
    return distance_to(other) < 0;
  }

  bool operator>=(const zip_pointer &other) const {
    return distance_to(other) <= 0;
  }

  zip_pointer &operator--() {
    decrement();
    return *this;
  }

  zip_pointer operator--(int) {
    zip_pointer temp(*this);
    decrement();
    return temp; // return saved state
  }

  reference operator*() const { return dereference(); }

  void increment() {
    detail::zip_helper<tuple_type>::increment_impl(data, index_type());
  }

  void decrement() {
    detail::zip_helper<tuple_type>::decrement_impl(data, index_type());
  }

  bool equal(zip_pointer const &other) const {
    return std::get<0>(other.data) == std::get<0>(data);
  }

  reference dereference() const {
    return reference(
        detail::zip_helper<tuple_type>::make_reference(data, index_type()));
  }

  difference_type distance_to(zip_pointer const &other) const {
    return std::get<0>(other.data) - std::get<0>(data);
  }

  void advance(difference_type n) {
    detail::zip_helper<tuple_type>::advance_impl(data, n, index_type());
  }

  tuple_type data;
};

#ifdef HAVE_THRUST
///
/// @brief specialisation of @ref zip_pointer for `thrust::tuple`
///
/// Use device only code in this specialisation
///
template <typename MplVector, typename TT1, typename TT2, typename TT3,
          typename TT4, typename TT5, typename TT6, typename TT7, typename TT8,
          typename TT9>
struct zip_pointer<thrust::tuple<TT1, TT2, TT3, TT4, TT5, TT6, TT7, TT8, TT9>,
                   MplVector> {
  typedef thrust::tuple<TT1, TT2, TT3, TT4, TT5, TT6, TT7, TT8, TT9> tuple_type;
  typedef MplVector mpl_vector_type;

  template <typename T>
  using elem_by_type = detail::get_elem_by_type<T, mpl_vector_type>;
  template <typename T>
  using return_type = typename detail::getter_helper<
      tuple_type>::template return_type<elem_by_type<T>::index>;

  typedef getter_type<typename detail::zip_helper<tuple_type>::tuple_reference,
                      mpl_vector_type>
      reference;
  typedef
      typename detail::zip_helper<tuple_type>::difference_type difference_type;
  typedef typename detail::zip_helper<tuple_type>::index_type index_type;

  CUDA_HOST_DEVICE
  zip_pointer() {}
  CUDA_HOST_DEVICE
  explicit zip_pointer(const tuple_type &data) : data(data) {}
  CUDA_HOST_DEVICE
  zip_pointer(const zip_pointer &other) : data(other.data) {}
  CUDA_HOST_DEVICE
  zip_pointer(zip_pointer &&other) : data(std::move(other.data)) {}

  template <typename tuple_type2,
            typename = typename std::enable_if<
                std::is_convertible<tuple_type2, tuple_type>::value>::type>
  CUDA_HOST_DEVICE
  zip_pointer(const zip_pointer<tuple_type2, mpl_vector_type> &other)
      : data(other.data) {}

  template <typename tuple_type2,
            typename = typename std::enable_if<
                std::is_convertible<tuple_type2, tuple_type>::value>::type>
  CUDA_HOST_DEVICE
  zip_pointer(const zip_pointer<tuple_type2, mpl_vector_type> &&other)
      : data(std::move(other.data)) {}

  template <typename T1, typename T2, typename... T3>
  CUDA_HOST_DEVICE zip_pointer(T1 &&arg1, T2 &&arg2, T3 &&... args)
      : data(std::forward<T1>(arg1), std::forward<T2>(arg2),
             std::forward<T3>(args)...) {}

  CUDA_HOST_DEVICE
  zip_pointer &operator=(const zip_pointer &other) {
    data = other.data;
    return *this;
  }
  CUDA_HOST_DEVICE
  zip_pointer &operator=(zip_pointer &&other) {
    data = std::move(other.data);
    return *this;
  }
  template <typename T1, typename T2>
  CUDA_HOST_DEVICE zip_pointer &operator=(const zip_pointer<T1, T2> &other) {
    data = other.data;
    return *this;
  }
  template <typename T1, typename T2>
  CUDA_HOST_DEVICE zip_pointer &operator=(zip_pointer<T1, T2> &&other) {
    data = std::move(other.data);
    return *this;
  }

  CUDA_HOST_DEVICE
  void swap(zip_pointer &other) { data.swap(other.data); }

  CUDA_HOST_DEVICE
  const tuple_type &get_tuple() const { return data; }
  CUDA_HOST_DEVICE
  tuple_type &get_tuple() { return data; }

  CUDA_HOST_DEVICE
  bool operator==(const zip_pointer &other) const { return equal(other); }

  CUDA_HOST_DEVICE
  zip_pointer &operator++() {
    advance(1);
    return *this;
  }

  CUDA_HOST_DEVICE
  zip_pointer operator++(int) {
    zip_pointer temp(*this);
    advance(1);
    return temp; // return saved state
  }

  CUDA_HOST_DEVICE
  zip_pointer operator+(const difference_type &n) const {
    zip_pointer ret(*this);
    ret.advance(n);
    return ret;
  }

  zip_pointer &operator+=(const difference_type &n) {
    advance(n);
    return *this;
  }

  CUDA_HOST_DEVICE
  zip_pointer operator-(const difference_type &n) const {
    zip_pointer ret(*this);
    ret.advance(-n);
    return ret;
  }

  CUDA_HOST_DEVICE
  difference_type operator-(const zip_pointer &other) const {
    return other.distance_to(*this);
  }

  CUDA_HOST_DEVICE
  zip_pointer &operator--() {
    decrement();
    return *this;
  }

  CUDA_HOST_DEVICE
  zip_pointer operator--(int) {
    zip_pointer temp(*this);
    decrement();
    return temp; // return saved state
  }

  CUDA_HOST_DEVICE
  bool operator<(const zip_pointer &other) const {
    return distance_to(other) > 0;
  }

  CUDA_HOST_DEVICE
  bool operator<=(const zip_pointer &other) const {
    return distance_to(other) >= 0;
  }

  CUDA_HOST_DEVICE
  bool operator>(const zip_pointer &other) const {
    return distance_to(other) < 0;
  }

  CUDA_HOST_DEVICE
  bool operator>=(const zip_pointer &other) const {
    return distance_to(other) <= 0;
  }

  CUDA_HOST_DEVICE
  reference operator*() const { return dereference(); }

  CUDA_HOST_DEVICE
  void increment() {
    detail::zip_helper<tuple_type>::increment_impl(data, index_type());
  }
  CUDA_HOST_DEVICE
  void decrement() {
    detail::zip_helper<tuple_type>::decrement_impl(data, index_type());
  }

  CUDA_HOST_DEVICE
  bool equal(zip_pointer const &other) const {
    return thrust::get<0>(other.data) == thrust::get<0>(data);
  }

  CUDA_HOST_DEVICE
  reference dereference() const {
    return reference(
        detail::zip_helper<tuple_type>::make_reference(data, index_type()));
  }

  CUDA_HOST_DEVICE
  difference_type distance_to(zip_pointer const &other) const {
    return thrust::get<0>(other.data) - thrust::get<0>(data);
  }

  CUDA_HOST_DEVICE
  void advance(difference_type n) {
    detail::zip_helper<tuple_type>::advance_impl(data, n, index_type());
  }

  tuple_type data;
};
#endif

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
    detail::getter_helper<tuple_type>::is_reference::value &&
    !detail::getter_helper<tuple_type2>::is_reference::value>::type
swap(getter_type<tuple_type, mpl_vector_type> x,
     getter_type<tuple_type2, mpl_vector_type> &y) {
  x.swap_via_tie(y.data, detail::getter_helper<tuple_type2>::index_type());
}

template <typename tuple_type, typename tuple_type2, typename mpl_vector_type>
typename std::enable_if<
    !detail::getter_helper<tuple_type>::is_reference::value &&
    detail::getter_helper<tuple_type2>::is_reference::value>::type
swap(getter_type<tuple_type, mpl_vector_type> &x,
     getter_type<tuple_type2, mpl_vector_type> y) {
  y.swap_via_tie(x.data, detail::getter_helper<tuple_type2>::index_type());
}

template <typename tuple_type, typename mpl_vector_type>
void swap(getter_type<tuple_type, mpl_vector_type> x,
          getter_type<tuple_type, mpl_vector_type> y) {
  y.swap(x);
}

/*
template <typename tuple_type, typename tuple_type2, typename
mpl_vector_type> typename std::enable_if<
    detail::getter_helper<tuple_type>::is_reference::value
    && detail::getter_helper<tuple_type2>::is_reference::value
    >::type
swap(getter_type<tuple_type,mpl_vector_type> x,
     getter_type<tuple_type2,mpl_vector_type> y) {
    y.swap(x);
}

template <typename tuple_type, typename tuple_type2, typename
mpl_vector_type> typename std::enable_if<
    !detail::getter_helper<tuple_type>::is_reference::value
    && !detail::getter_helper<tuple_type2>::is_reference::value
    >::type
swap(getter_type<tuple_type,mpl_vector_type>& x,
     getter_type<tuple_type2,mpl_vector_type>& y) {
    y.swap(x);
}
*/

template <typename iterator_tuple_type, typename mpl_vector_type>
class zip_iterator {};

///
/// @brief specialisation of @ref zip_iterator for `std::tuple`
///
/// Use host only code in this specialisation
///
template <typename mpl_vector_type, typename... Types>
class zip_iterator<std::tuple<Types...>, mpl_vector_type>
    : public detail::zip_iterator_base<std::tuple<Types...>,
                                       mpl_vector_type>::type {
  typedef std::tuple<Types...> iterator_tuple_type;

public:
  typedef iterator_tuple_type tuple_type;
  typedef typename detail::zip_iterator_base<
      tuple_type, mpl_vector_type>::value_type value_type;
  typedef
      typename detail::zip_iterator_base<tuple_type, mpl_vector_type>::reference
          reference;

  typedef typename detail::zip_helper<iterator_tuple_type>::difference_type
      difference_type;
  typedef typename detail::zip_helper<iterator_tuple_type>::iterator_category
      iterator_category;

  typedef getter_type<
      typename detail::zip_helper<iterator_tuple_type>::tuple_pointer,
      mpl_vector_type>
      pointer;

  // Note: Can't call it raw_pointer or thrust thinks it is a trivial
  // iterator!
  typedef zip_pointer<
      typename detail::zip_helper<iterator_tuple_type>::tuple_raw_pointer,
      mpl_vector_type>
      getter_raw_pointer;
  typedef getter_type<
      typename detail::zip_helper<iterator_tuple_type>::tuple_raw_reference,
      mpl_vector_type>
      getter_raw_reference;

  template <typename T>
  using elem_by_type = detail::get_elem_by_type<T, mpl_vector_type>;

  template <typename T> struct return_type {
    static const size_t N = elem_by_type<T>::index;
    typedef const typename detail::zip_helper<
        iterator_tuple_type>::template tuple_element<N>::type type;
  };

  zip_iterator() {}

  explicit zip_iterator(iterator_tuple_type iter) : iter(iter) {}

  template <typename... T> explicit zip_iterator(T... args) : iter(args...) {}

  CUDA_HOST_DEVICE
  const iterator_tuple_type &get_tuple() const {
#if defined(__CUDA_ARCH__)
    ERROR_CUDA("Cannot use `zip_iterator` in device code");
#endif
    return iter;
  }

  CUDA_HOST_DEVICE
  iterator_tuple_type &get_tuple() {
#if defined(__CUDA_ARCH__)
    ERROR_CUDA("Cannot use `zip_iterator` in device code");
#endif
    return iter;
  }

private:
  typedef
      typename detail::zip_helper<iterator_tuple_type>::index_type index_type;

  void increment() {
    detail::zip_helper<iterator_tuple_type>::increment_impl(iter, index_type());
  }

  void decrement() {
    detail::zip_helper<iterator_tuple_type>::decrement_impl(iter, index_type());
  }

  bool equal(zip_iterator const &other) const {
    return detail::get_impl<0>(other.iter) == detail::get_impl<0>(iter);
  }

  reference dereference() const {
    return reference(detail::zip_helper<iterator_tuple_type>::make_reference(
        iter, index_type()));
  }

  difference_type distance_to(zip_iterator const &other) const {
#if defined(__CUDA_ARCH__)
    ERROR_CUDA("Cannot use `zip_iterator` in device code");
#endif
    return detail::get_impl<0>(other.iter) - detail::get_impl<0>(iter);
  }

  void advance(difference_type n) {
    detail::zip_helper<iterator_tuple_type>::advance_impl(iter, n,
                                                          index_type());
  }

  iterator_tuple_type iter;
  friend typename detail::zip_helper<iterator_tuple_type>::iterator_core_access;
};

#ifdef HAVE_THRUST
///
/// @brief specialisation of @ref zip_iterator for `thrust::tuple`
///
/// Use device only code in this specialisation
///
template <typename mpl_vector_type, typename... Types>
class zip_iterator<thrust::tuple<Types...>, mpl_vector_type>
    : public detail::zip_iterator_base<thrust::tuple<Types...>,
                                       mpl_vector_type>::type {
  typedef thrust::tuple<Types...> iterator_tuple_type;

public:
  typedef iterator_tuple_type tuple_type;
  typedef typename detail::zip_iterator_base<
      tuple_type, mpl_vector_type>::value_type value_type;
  typedef
      typename detail::zip_iterator_base<tuple_type, mpl_vector_type>::reference
          reference;

  typedef typename detail::zip_helper<iterator_tuple_type>::difference_type
      difference_type;
  typedef typename detail::zip_helper<iterator_tuple_type>::iterator_category
      iterator_category;

  typedef getter_type<
      typename detail::zip_helper<iterator_tuple_type>::tuple_pointer,
      mpl_vector_type>
      pointer;

  // Note: Can't call it raw_pointer or thrust thinks it is a trivial
  // iterator!
  typedef zip_pointer<
      typename detail::zip_helper<iterator_tuple_type>::tuple_raw_pointer,
      mpl_vector_type>
      getter_raw_pointer;
  typedef getter_type<
      typename detail::zip_helper<iterator_tuple_type>::tuple_raw_reference,
      mpl_vector_type>
      getter_raw_reference;

  template <typename T>
  using elem_by_type = detail::get_elem_by_type<T, mpl_vector_type>;

  template <typename T> struct return_type {
    static const size_t N = elem_by_type<T>::index;
    typedef const typename detail::zip_helper<
        iterator_tuple_type>::template tuple_element<N>::type type;
  };

  CUDA_HOST_DEVICE
  zip_iterator() {}

  CUDA_HOST_DEVICE
  explicit zip_iterator(iterator_tuple_type iter) : iter(iter) {}

  template <typename... T>
  CUDA_HOST_DEVICE explicit zip_iterator(T... args) : iter(args...) {}

  CUDA_HOST_DEVICE
  const iterator_tuple_type &get_tuple() const { return iter; }

  CUDA_HOST_DEVICE
  iterator_tuple_type &get_tuple() { return iter; }

private:
  typedef
      typename detail::zip_helper<iterator_tuple_type>::index_type index_type;

  CUDA_HOST_DEVICE
  void increment() {
    detail::zip_helper<iterator_tuple_type>::increment_impl(iter, index_type());
  }

  CUDA_HOST_DEVICE
  void decrement() {
    detail::zip_helper<iterator_tuple_type>::decrement_impl(iter, index_type());
  }

  CUDA_HOST_DEVICE
  bool equal(zip_iterator const &other) const {
    return detail::get_impl<0>(other.iter) == detail::get_impl<0>(iter);
  }

  CUDA_HOST_DEVICE
  reference dereference() const {
    return reference(detail::zip_helper<iterator_tuple_type>::make_reference(
        iter, index_type()));
  }

  CUDA_HOST_DEVICE
  difference_type distance_to(zip_iterator const &other) const {
    return thrust::get<0>(other.iter) - thrust::get<0>(iter);
  }

  CUDA_HOST_DEVICE
  void advance(difference_type n) {
    detail::zip_helper<iterator_tuple_type>::advance_impl(iter, n,
                                                          index_type());
  }

  iterator_tuple_type iter;
  friend typename detail::zip_helper<iterator_tuple_type>::iterator_core_access;
};
#endif

///
/// @brief convert an generic iterator to a @ref raw_pointer
///
/// @tparam Iterator the iterator type to be converted. Can be a normal STL
/// iterator
///         or a @zip_iterator
/// @param arg  the actual iterator object to be converted
/// @return decltype(detail::iterator_to_raw_pointer(arg,typename
/// detail::is_zip_iterator<Iterator>::type()))
///
template <typename Iterator>
auto iterator_to_raw_pointer(const Iterator &arg)
    -> decltype(detail::iterator_to_raw_pointer(
        arg, typename detail::is_zip_iterator<Iterator>::type()))

{
  return detail::iterator_to_raw_pointer(
      arg, typename detail::is_zip_iterator<Iterator>::type());
}

/*
template <typename Iterator>
auto pointer_to_raw_pointer(const Iterator& arg) ->
    decltype(detail::pointer_to_raw_pointer(arg,
                typename detail::is_zip_iterator<Iterator>::type()
                           ))

{
    return detail::pointer_to_raw_pointer(arg,
            typename detail::is_zip_iterator<Iterator>::type()
                           );
}
*/

//
// Particle getters/setters
//
//

/// get a variable from a particle \p arg
/// \param T Variable type
/// \param arg the particle
/// \return a const reference to a T::value_type holding the variable data
///

/*
template<typename T, typename ValueType>
CUDA_HOST_DEVICE
auto get(const ValueType& arg) ->
    decltype(detail::get_impl<ValueType::template
elem_by_type<T>::index>(arg.get_tuple()))
{
    //std::cout << "get const reference" << std::endl;
    return detail::get_impl<ValueType::template
elem_by_type<T>::index>(arg.get_tuple());
    //return arg.template get<T>();
}
*/

///
/// @brief get the value of a variable from a const @ref getter_type,
///        @ref zip_iterator, @ref zip_pointer, or @ref Particles
///
/// @tparam T the variable type to get, @see ABORIA_VARIABLE
/// @tparam ValueType the @ref getter_type, @ref zip_iterator, @ref
/// zip_pointer,
///                   or @ref Particles
/// @param arg the object with type ValueType
/// @return CUDA_HOST_DEVICE typename ValueType::template return_type<T>
/// ::type const&
///
template <typename T, typename ValueType>
CUDA_HOST_DEVICE typename ValueType::template return_type<T>::type const &
get(const ValueType &arg) {
  // std::cout << "get reference" << std::endl;
  return detail::get_impl<ValueType::template elem_by_type<T>::index>(
      arg.get_tuple());
  // return arg.template get<T>();
}

///
/// @brief get the value of a variable from a @ref getter_type,
///        @ref zip_iterator, @ref zip_pointer, or @ref Particles
///
/// @tparam T the variable type to get, @see ABORIA_VARIABLE
/// @tparam ValueType the @ref getter_type, @ref zip_iterator, @ref
/// zip_pointer,
///                   or @ref Particles
/// @param arg the object with type ValueType
/// @return CUDA_HOST_DEVICE typename ValueType::template return_type<T>
/// ::type const&
///
template <typename T, typename ValueType>
CUDA_HOST_DEVICE typename ValueType::template return_type<T>::type &
get(ValueType &arg) {
  // std::cout << "get reference" << std::endl;
  return detail::get_impl<ValueType::template elem_by_type<T>::index>(
      arg.get_tuple());
  // return arg.template get<T>();
}

///
/// @brief get the value of a variable from a rvalue @ref getter_type,
///        @ref zip_iterator, @ref zip_pointer, or @ref Particles
///
/// @tparam T the variable type to get, @see ABORIA_VARIABLE
/// @tparam ValueType the @ref getter_type, @ref zip_iterator, @ref
/// zip_pointer,
///                   or @ref Particles
/// @param arg the object with type ValueType
/// @return CUDA_HOST_DEVICE typename ValueType::template return_type<T>
/// ::type const&
///
template <typename T, typename ValueType>
CUDA_HOST_DEVICE typename ValueType::template return_type<T>::type &
get(ValueType &&arg) {
  // std::cout << "get reference" << std::endl;
  return detail::get_impl<ValueType::template elem_by_type<T>::index>(
      arg.get_tuple());
  // return arg.template get<T>();
}

///
/// @brief get the value of an indexed variable from a const @ref
/// getter_type,
///        @ref zip_iterator, @ref zip_pointer, or @ref Particles
///
/// @tparam N the index of the variable to get. The order of variables is
/// set to
///           (position,id,alive,generator,user_variable1,user_variable2,...)
/// @tparam ValueType the @ref getter_type, @ref zip_iterator, @ref
/// zip_pointer,
///                   or @ref Particles
/// @param arg the object with type ValueType
/// @return CUDA_HOST_DEVICE const typename detail::getter_helper<typename
/// ValueType::tuple_type> ::template return_type<N> ::type&
///
template <unsigned int N, typename ValueType>
CUDA_HOST_DEVICE const typename detail::getter_helper<
    typename ValueType::tuple_type>::template return_type<N>::type &
get_by_index(const ValueType &arg) {
  return detail::get_impl<N>(arg.get_tuple());
}

///
/// @brief get the value of an indexed variable from a lvalue @ref
/// getter_type,
///        @ref zip_iterator, @ref zip_pointer, or @ref Particles
///
/// @tparam N the index of the variable to get. The order of variables is
/// set to
///           (position,id,alive,generator,user_variable1,user_variable2,...)
/// @tparam ValueType the @ref getter_type, @ref zip_iterator, @ref
/// zip_pointer,
///                   or @ref Particles
/// @param arg the object with type ValueType
/// @return CUDA_HOST_DEVICE const typename detail::getter_helper<typename
/// ValueType::tuple_type> ::template return_type<N> ::type&
///
template <unsigned int N, typename ValueType>
CUDA_HOST_DEVICE typename detail::getter_helper<
    typename ValueType::tuple_type>::template return_type<N>::type &
get_by_index(ValueType &arg) {
  return detail::get_impl<N>(arg.get_tuple());
}

///
/// @brief get the value of an indexed variable from a rvalue @ref
/// getter_type,
///        @ref zip_iterator, @ref zip_pointer, or @ref Particles
///
/// @tparam N the index of the variable to get. The order of variables is
/// set to
///           (position,id,alive,generator,user_variable1,user_variable2,...)
/// @tparam ValueType the @ref getter_type, @ref zip_iterator, @ref
/// zip_pointer,
///                   or @ref Particles
/// @param arg the object with type ValueType
/// @return CUDA_HOST_DEVICE const typename detail::getter_helper<typename
/// ValueType::tuple_type> ::template return_type<N> ::type&
///
template <unsigned int N, typename ValueType>
CUDA_HOST_DEVICE typename detail::getter_helper<
    typename ValueType::tuple_type>::template return_type<N>::type &
get_by_index(ValueType &&arg) {
  return detail::get_impl<N>(arg.get_tuple());
}

} // namespace Aboria

#endif // GET_H_
