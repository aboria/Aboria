
#ifndef GETTERTYPE_DETAIL_H_ 
#define GETTERTYPE_DETAIL_H_ 

#include "Helpers.h"
#include "Log.h"

namespace Aboria {
namespace detail {



template <bool Pointers, typename MplVector, typename Tuple> 
struct getter_type_base{};

/*
 * specialisation for std::tuple. This needs to be host only
 * TODO: api is not good, consider doing an iterator_facade type thing
 */
template <typename MplVector, typename ... Types> 
struct getter_type_base<false, MplVector, std::tuple<Types...>>{
    typedef std::tuple<Types...> tuple_type;
    typedef MplVector mpl_vector_type;
 
    typedef typename detail::getter_helper<tuple_type>::tuple_reference tuple_reference;
    template <typename T>
    using elem_by_type = detail::get_elem_by_type<T,mpl_vector_type>;
    template <typename T>
    using return_type = typename detail::getter_helper<tuple_type>::template return_type<elem_by_type<T>::index>;

    //typedef typename detail::zip_helper<tuple_type>::pointer pointer;

    getter_type_base() {}

    explicit getter_type_base(const tuple_type& data):data(data) {}

    getter_type_base(const getter_type_base& other):data(other.data) {}

    getter_type_base(getter_type_base&& other):data(std::move(other.data)) {}

    template <typename T=tuple_reference, typename = typename
    std::enable_if<
        !std::is_same<T,tuple_type>::value
        >::type>
    getter_type_base(const getter_type_base<false,mpl_vector_type,tuple_reference>& other):
        data(other.data) {}

    template <typename T=tuple_reference, typename = typename
    std::enable_if<
        !std::is_same<T,tuple_type>::value
        >::type>
    getter_type_base(getter_type_base<false,mpl_vector_type, tuple_reference>&& other):
        data(std::move(other.data)) {}

    template <typename T, typename = typename
    std::enable_if<
        (!std::is_same<T,tuple_type>::value)
        && (!std::is_same<T,tuple_reference>::value)
        >::type>
    getter_type_base(const getter_type_base<false,mpl_vector_type,T>& other):
        data(other.data) {}

    template <typename T, typename = typename
    std::enable_if<
        (!std::is_same<T,tuple_type>::value)
        && (!std::is_same<T,tuple_reference>::value)
        >::type>
    getter_type_base(getter_type_base<false,mpl_vector_type,T>&& other):
        data(std::move(other.data)) {}

    
    template <typename T1, typename T2, typename... T3>
    getter_type_base(T1&& arg1, T2&& arg2, T3&&... args):data(std::forward<T1>(arg1),std::forward<T2>(arg2),std::forward<T3>(args)...) {}

    template <typename T>
    int throw_away(const T& in) {
        return 0;
    }

    template<typename Tuple, std::size_t... I>
    void copy_impl(const Tuple& other_data, detail::index_sequence<I...>) {
        int dummy[] = { 0, throw_away(std::get<I>(data) = std::get<I>(other_data))... };
        static_cast<void>(dummy);
    }

    getter_type_base& operator=( const getter_type_base& other ) {
        //copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
        data = other.data;
        return *this;
    }
    getter_type_base& operator=( getter_type_base&& other ) {
        //copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
        data = other.data;
        return *this;
    }

    template <typename T=tuple_reference, typename = typename
    std::enable_if<
        !std::is_same<T,tuple_type>::value
        >::type>
    getter_type_base& operator=( const getter_type_base<false,mpl_vector_type,tuple_reference>& other) {
        //copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
        data = other.data;
        return *this;
    }

    template <typename T=tuple_reference, typename = typename
    std::enable_if<
        !std::is_same<T,tuple_type>::value
        >::type>
    getter_type_base& operator=(getter_type_base<false,mpl_vector_type,tuple_reference>&& other) {
        //copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
        data = other.data;
        return *this;
    }

    template <typename T, typename = typename
    std::enable_if<
        (!std::is_same<T,tuple_reference>::value)
        && (!std::is_same<T,tuple_type>::value)
        >::type>
    getter_type_base& operator=( const getter_type_base<false,mpl_vector_type,T>& other) {
        //copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
        data = other.data;
        return *this;
    }

    template <typename T, typename = typename
    std::enable_if<
        (!std::is_same<T,tuple_reference>::value)
        && (!std::is_same<T,tuple_type>::value)
        >::type>
    getter_type_base& operator=( getter_type_base<false,mpl_vector_type,T>&& other) {
        //copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
        data = other.data;
        return *this;
    }

#if defined(__aboria_use_thrust_algorithms__) || defined(__CUDACC__)
    template <typename T1, typename T2, typename PointerType, typename DerivedType> 
    getter_type_base& operator=( const thrust::reference<getter_type_base<false,T1,T2>,PointerType,DerivedType>& other) {
        data = static_cast<getter_type_base<false,T1,T2>>(other).data;
        return *this;
    }
#endif

    template <typename T1, typename T2> 
    bool operator==( const getter_type_base<false,T1,T2>& other) {
        return data == other.data;
    }
    
    void swap(getter_type_base &other) {
        data.swap(other.data);
    }

    template <typename tuple_type2,std::size_t... I>
    void swap_via_tie(tuple_type2 &tuple, detail::index_sequence<I...>) {
        tuple_type tmp = std::tie(std::get<I>(tuple)...);
        data.swap(tmp);
    }

    CUDA_HOST_DEVICE
    const tuple_type & get_tuple() const { 
        #if defined(__CUDA_ARCH__)
        ERROR_CUDA("Cannot use `getter_type_base<std::tuple>` in device code");
        #endif
        return data;
    }

    CUDA_HOST_DEVICE
    tuple_type & get_tuple() { 
        #if defined(__CUDA_ARCH__)
        ERROR_CUDA("Cannot use `getter_type_base<std::tuple>` in device code");
        #endif
        return data; 
    }


    tuple_type data;
};

/*
 * specialisation for thrust::tuple.
 * TODO: api is not good, consider doing an iterator_facade type thing
 */
#ifdef __aboria_have_thrust__ 
template <typename MplVector, typename TT1, typename TT2, typename TT3, typename TT4, typename TT5, typename TT6, typename TT7, typename TT8, typename TT9> 
struct getter_type_base<false, MplVector, thrust::tuple<TT1,TT2,TT3,TT4,TT5,TT6,TT7,TT8,TT9>>{
    typedef thrust::tuple<TT1,TT2,TT3,TT4,TT5,TT6,TT7,TT8,TT9> tuple_type;
    typedef MplVector mpl_vector_type;
 
    typedef typename detail::getter_helper<tuple_type>::tuple_reference tuple_reference;
    typedef typename detail::getter_helper<tuple_type>::tuple_device_reference tuple_device_reference;
    typedef typename detail::getter_helper<tuple_type>::index_type index_type;
    template <typename T>
    using elem_by_type = detail::get_elem_by_type<T,mpl_vector_type>;
    template <typename T>
    using return_type = typename detail::getter_helper<tuple_type>::template return_type<elem_by_type<T>::index>;

    //typedef typename detail::zip_helper<tuple_type>::pointer pointer;

    CUDA_HOST_DEVICE
    getter_type_base() {}

    CUDA_HOST_DEVICE
    explicit getter_type_base(const tuple_type& data):data(data) {}

    CUDA_HOST_DEVICE
    getter_type_base(const getter_type_base& other):data(other.data) {}

    template <typename T=tuple_reference, typename = typename
    std::enable_if<
        !std::is_same<T,tuple_type>::value
        >::type>
    CUDA_HOST_DEVICE
    getter_type_base(const getter_type_base<false,mpl_vector_type,tuple_reference>& other):
        data(other.data) {}

    template <typename T=tuple_reference, typename = typename
    std::enable_if<
        !std::is_same<T,tuple_type>::value
        >::type>
    CUDA_HOST_DEVICE
    getter_type_base(getter_type_base<false,mpl_vector_type,tuple_reference>&& other):
        data(std::move(other.data)) {}

    template <typename T=tuple_device_reference, typename = typename
    std::enable_if<
        !std::is_same<T,tuple_type>::value
        >::type>
    CUDA_HOST_DEVICE
    getter_type_base(const getter_type_base<false,mpl_vector_type,tuple_device_reference>& other):
        data(detail::getter_helper<tuple_type>::raw_reference_cast(other.data,index_type())) {}

    template <typename T=tuple_device_reference, typename = typename
    std::enable_if<
        !std::is_same<T,tuple_type>::value
        >::type>
    CUDA_HOST_DEVICE
    getter_type_base(getter_type_base<false,mpl_vector_type,tuple_device_reference>&& other):
        data(detail::getter_helper<tuple_type>::raw_reference_cast(std::move(other.data),index_type())) {}

    template <typename T, typename = typename
    std::enable_if<
        (!std::is_same<T,tuple_type>::value)
        && (!std::is_same<T,tuple_reference>::value)
        >::type>
    CUDA_HOST_DEVICE
    getter_type_base(const getter_type_base<false,mpl_vector_type,T>& other):
        data(other.data) {}

    template <typename T, typename = typename
    std::enable_if<
        (!std::is_same<T,tuple_type>::value)
        && (!std::is_same<T,tuple_reference>::value)
        >::type>
    CUDA_HOST_DEVICE
    getter_type_base(getter_type_base<false,mpl_vector_type,T>&& other):
        data(std::move(other.data)) {}

    
    template <typename T1, typename T2, typename... T3>
    CUDA_HOST_DEVICE
    getter_type_base(T1&& arg1, T2&& arg2, T3&&... args):data(std::forward<T1>(arg1),std::forward<T2>(arg2),std::forward<T3>(args)...) {}

    template <typename T>
    CUDA_HOST_DEVICE
    int throw_away(const T& in) {
        return 0;
    }

    template<typename Tuple, std::size_t... I>
    CUDA_HOST_DEVICE
    void copy_impl(const Tuple& other_data, detail::index_sequence<I...>) {
        int dummy[] = { 0, throw_away(thrust::get<I>(data) = thrust::get<I>(other_data))... };
        static_cast<void>(dummy);
    }

    CUDA_HOST_DEVICE
    getter_type_base& operator=( const getter_type_base& other ) {
        //copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
        data = other.data;
        return *this;
    }

    CUDA_HOST_DEVICE
    getter_type_base& operator=( getter_type_base&& other ) {
        //copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
        data = other.data;
        return *this;
    }

    template <typename T=tuple_reference, typename = typename
    std::enable_if<
        !std::is_same<T,tuple_type>::value
        >::type>
    CUDA_HOST_DEVICE
    getter_type_base& operator=( const getter_type_base<false,mpl_vector_type,tuple_reference>& other) {
        //copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
        data = other.data;
        return *this;
    }

    template <typename T=tuple_reference, typename = typename
    std::enable_if<
        !std::is_same<T,tuple_type>::value
        >::type>
    CUDA_HOST_DEVICE
    getter_type_base& operator=(getter_type_base<false,mpl_vector_type,tuple_reference>&& other) {
        //copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
        data = other.data;
        return *this;
    }

    template <typename T, typename = typename
    std::enable_if<
        (!std::is_same<T,tuple_reference>::value)
        && (!std::is_same<T,tuple_type>::value)
        >::type>
    CUDA_HOST_DEVICE
    getter_type_base& operator=( const getter_type_base<false,mpl_vector_type,T>& other) {
        //copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
        data = other.data;
        return *this;
    }

    template <typename T, typename = typename
    std::enable_if<
        (!std::is_same<T,tuple_reference>::value)
        && (!std::is_same<T,tuple_type>::value)
        >::type>
    CUDA_HOST_DEVICE
    getter_type_base& operator=( getter_type_base<false,mpl_vector_type,T>&& other) {
        //copy_impl(other.data,detail::make_index_sequence<std::tuple_size<tuple_type>::value>());
        data = other.data;
        return *this;
    }

    template <typename T1, typename T2> 
    CUDA_HOST_DEVICE
    bool operator==( const getter_type_base<false,T1,T2>& other) {
        return data == other.data;
    }
    
    CUDA_HOST_DEVICE
    void swap(getter_type_base &other) {
        data.swap(other.data);
    }

    template <typename tuple_type2,std::size_t... I>
    CUDA_HOST_DEVICE
    void swap_via_tie(tuple_type2 &tuple, detail::index_sequence<I...>) {
        tuple_type tmp = std::tie(thrust::get<I>(tuple)...);
        data.swap(tmp);
    }

    CUDA_HOST_DEVICE
    const tuple_type & get_tuple() const { return data; }

    CUDA_HOST_DEVICE
    tuple_type & get_tuple() { return data; }

    tuple_type data;
};
#endif

/*
 * specialisation for tuple of pointers. Note hack to detect tuple of pointers.. could be better?
 * TODO: api is not good, consider doing an iterator_facade type thing
 */
template <typename MplVector, typename ... Types> 
struct getter_type_base<true, MplVector, std::tuple<Types*...>>{
    typedef std::tuple<Types*...> tuple_type;
    typedef MplVector mpl_vector_type;

    template <typename T>
    using elem_by_type = detail::get_elem_by_type<T,mpl_vector_type>;
    template <typename T>
    using return_type = typename detail::getter_helper<tuple_type>::template return_type<elem_by_type<T>::index>;

    typedef getter_type_base<false,mpl_vector_type,
            typename detail::zip_helper<tuple_type>::tuple_reference> reference;
    typedef typename detail::zip_helper<tuple_type>::difference_type difference_type;
    typedef typename detail::zip_helper<tuple_type>::index_type index_type;

    getter_type_base() {}

    explicit getter_type_base(const tuple_type& data):data(data) {}

    getter_type_base(const getter_type_base& other):data(other.data) {}

    getter_type_base(getter_type_base&& other):data(std::move(other.data)) {}

    template <typename tuple_type2, typename = typename
    std::enable_if<std::is_convertible<tuple_type2,tuple_type>::value>::type>
    getter_type_base(const getter_type_base<true,mpl_vector_type,tuple_type2>& other):data(other.data) {}

    template <typename tuple_type2, typename = typename
    std::enable_if<std::is_convertible<tuple_type2,tuple_type>::value>::type>
    getter_type_base(const getter_type_base<true,mpl_vector_type,tuple_type2>&& other):data(std::move(other.data)) {}

    template <typename T1, typename T2, typename... T3>
    getter_type_base(T1&& arg1, T2&& arg2, T3&&... args):data(std::forward<T1>(arg1),std::forward<T2>(arg2),std::forward<T3>(args)...) {}

    getter_type_base& operator=( const getter_type_base& other ) {
        data = other.data;
        return *this;
    }
    getter_type_base& operator=( getter_type_base&& other ) {
        data = std::move(other.data);
        return *this;
    }
    template <typename T1, typename T2> 
    getter_type_base& operator=( const getter_type_base<true,T1,T2>& other) {
        data = other.data;
        return *this;
    }
    template <typename T1, typename T2> 
    getter_type_base& operator=( getter_type_base<true,T1,T2>&& other) {
        data = std::move(other.data);
        return *this;
    }
    
    void swap(getter_type_base &other) {
        data.swap(other.data);
    }

    const tuple_type & get_tuple() const { return data; }
    CUDA_HOST_DEVICE
    tuple_type & get_tuple() { return data; }

    
    bool operator==(const getter_type_base& other) const {
        return equal(other);
    }

    getter_type_base& operator++() {
        advance(1);
        return *this;
    }
    
    getter_type_base operator++(int) {
        getter_type_base temp(*this);
        advance(1);
        return temp; // return saved state
    }

    getter_type_base operator+(const difference_type& n) const {
        getter_type_base ret(*this);
        ret.advance(n);
        return ret;
    }

    difference_type operator-(const getter_type_base& other) const {
        return other.distance_to(*this);
    }

    getter_type_base& operator--() {
        decrement();
        return *this;
    }
    
    getter_type_base operator--(int) {
        getter_type_base temp(*this);
        decrement();
        return temp; // return saved state
    }

    reference operator*() const {
        return dereference();
    }

    void increment() { 
        detail::zip_helper<tuple_type>::increment_impl(data,index_type()); 
    }

    void decrement() { 
        detail::zip_helper<tuple_type>::decrement_impl(data,index_type()); 
    }

    bool equal(getter_type_base const& other) const { 
        return std::get<0>(other.data) == std::get<0>(data);
    }

    reference dereference() const { 
        return reference(
                detail::zip_helper<tuple_type>::make_reference(
                    data,index_type())
                ); 
    }

    difference_type distance_to(getter_type_base const& other) const { 
        return std::get<0>(other.data) - std::get<0>(data);
    }

    void advance(difference_type n) { 
        detail::zip_helper<tuple_type>::advance_impl(data,n,index_type()); 
    }

    tuple_type data;
};

#ifdef __aboria_have_thrust__ 
template <typename MplVector, typename TT1, typename TT2, typename TT3, typename TT4, typename TT5, typename TT6, typename TT7, typename TT8, typename TT9> 
struct getter_type_base<true, MplVector, thrust::tuple<TT1,TT2,TT3,TT4,TT5,TT6,TT7,TT8,TT9>>{
    typedef thrust::tuple<TT1,TT2,TT3,TT4,TT5,TT6,TT7,TT8,TT9> tuple_type;
    typedef MplVector mpl_vector_type;

    template <typename T>
    using elem_by_type = detail::get_elem_by_type<T,mpl_vector_type>;
    template <typename T>
    using return_type = typename detail::getter_helper<tuple_type>::template return_type<elem_by_type<T>::index>;

    typedef getter_type_base<false,mpl_vector_type,
            typename detail::zip_helper<tuple_type>::tuple_reference> reference;
    typedef typename detail::zip_helper<tuple_type>::difference_type difference_type;
    typedef typename detail::zip_helper<tuple_type>::index_type index_type;


    CUDA_HOST_DEVICE
    getter_type_base() {}
    CUDA_HOST_DEVICE
    explicit getter_type_base(const tuple_type& data):data(data) {}
    CUDA_HOST_DEVICE
    getter_type_base(const getter_type_base& other):data(other.data) {}
    CUDA_HOST_DEVICE
    getter_type_base(getter_type_base&& other):data(std::move(other.data)) {}

    template <typename tuple_type2, typename = typename
    std::enable_if<std::is_convertible<tuple_type2,tuple_type>::value>::type>
    CUDA_HOST_DEVICE
    getter_type_base(const getter_type_base<true,mpl_vector_type,tuple_type2>& other):data(other.data) {}

    template <typename tuple_type2, typename = typename
    std::enable_if<std::is_convertible<tuple_type2,tuple_type>::value>::type>
    CUDA_HOST_DEVICE
    getter_type_base(const getter_type_base<true,mpl_vector_type,tuple_type2>&& other):data(std::move(other.data)) {}

    template <typename T1, typename T2, typename... T3>
    CUDA_HOST_DEVICE
    getter_type_base(T1&& arg1, T2&& arg2, T3&&... args):data(std::forward<T1>(arg1),std::forward<T2>(arg2),std::forward<T3>(args)...) {}

    CUDA_HOST_DEVICE
    getter_type_base& operator=( const getter_type_base& other ) {
        data = other.data;
        return *this;
    }
    CUDA_HOST_DEVICE
    getter_type_base& operator=( getter_type_base&& other ) {
        data = std::move(other.data);
        return *this;
    }
    template <typename T1, typename T2> 
    CUDA_HOST_DEVICE
    getter_type_base& operator=( const getter_type_base<true,T1,T2>& other) {
        data = other.data;
        return *this;
    }
    template <typename T1, typename T2> 
    CUDA_HOST_DEVICE
    getter_type_base& operator=( getter_type_base<true,T1,T2>&& other) {
        data = std::move(other.data);
        return *this;
    }
    
    CUDA_HOST_DEVICE
    void swap(getter_type_base &other) {
        data.swap(other.data);
    }

    CUDA_HOST_DEVICE
    const tuple_type & get_tuple() const { return data; }
    CUDA_HOST_DEVICE
    tuple_type & get_tuple() { return data; }

    
    CUDA_HOST_DEVICE
    bool operator==(const getter_type_base& other) const {
        return equal(other);
    }

    CUDA_HOST_DEVICE
    getter_type_base& operator++() {
        advance(1);
        return *this;
    }
    
    CUDA_HOST_DEVICE
    getter_type_base operator++(int) {
        getter_type_base temp(*this);
        advance(1);
        return temp; // return saved state
    }


    CUDA_HOST_DEVICE
    getter_type_base operator+(const difference_type& n) const {
        getter_type_base ret(*this);
        ret.advance(n);
        return ret;
    }

    CUDA_HOST_DEVICE
    difference_type operator-(const getter_type_base& other) const {
        return other.distance_to(*this);
    }

    CUDA_HOST_DEVICE
    getter_type_base& operator--() {
        decrement();
        return *this;
    }
    
    CUDA_HOST_DEVICE
    getter_type_base operator--(int) {
        getter_type_base temp(*this);
        decrement();
        return temp; // return saved state
    }


    CUDA_HOST_DEVICE
    reference operator*() const {
        return dereference();
    }

    CUDA_HOST_DEVICE
    void increment() { 
        detail::zip_helper<tuple_type>::increment_impl(data,index_type()); 
    }
    CUDA_HOST_DEVICE
    void decrement() { 
        detail::zip_helper<tuple_type>::decrement_impl(data,index_type()); 
    }

    CUDA_HOST_DEVICE
    bool equal(getter_type_base const& other) const { 
        return thrust::get<0>(other.data) == thrust::get<0>(data);
    }

    CUDA_HOST_DEVICE
    reference dereference() const { 
        return reference(
                detail::zip_helper<tuple_type>::make_reference(
                    data,index_type())
                ); 
    }

    CUDA_HOST_DEVICE
    difference_type distance_to(getter_type_base const& other) const { 
        return thrust::get<0>(other.data) - thrust::get<0>(data);
    }

    CUDA_HOST_DEVICE
    void advance(difference_type n) { 
        detail::zip_helper<tuple_type>::advance_impl(data,n,index_type()); 
    }

    tuple_type data;
};
#endif


template<typename Tuple>
struct pointer_check {};

template<typename FirstType, typename ... Types>
struct pointer_check<std::tuple<FirstType,Types...>> {
    static const bool value = std::is_pointer<FirstType>::value;
};

#ifdef __aboria_have_thrust__ 
template<typename FirstType, typename ... Types>
struct pointer_check<thrust::tuple<FirstType,Types...>> {
    static const bool value = std::is_pointer<FirstType>::value;
};
#endif

}
}

#endif
