#ifndef CUDA_INCLUDE_H_ 
#define CUDA_INCLUDE_H_ 

#if defined(__aboria_use_thrust_algorithms__) || defined(__CUDACC__)
    #include <thrust/device_vector.h>
    #include <thrust/host_vector.h>
    #include <thrust/sort.h>
    #include <thrust/binary_search.h>
    #include <thrust/iterator/iterator_facade.h>
    #include <thrust/sequence.h>
    #include <thrust/transform_scan.h>
    #include <nppdefs.h>
#endif


#if defined(__CUDACC__)
    #define CUDA_HOST_DEVICE __host__ __device__

    #if not defined(__aboria_use_thrust_algorithms__)
        #define __aboria_use_thrust_algorithms__
    #endif

    #define __aboria_hd_warning_disable__ \
    #pragma hd_warning_disable

    // if compiling with cuda compiler use cuda's tuple and iterator_facade
    namespace tuple_ns = thrust;
    namespace iterator_facade_ns = thrust;
#else
    #define CUDA_HOST_DEVICE 
    #define __aboria_hd_warning_disable__
    
    #include <boost/iterator/iterator_facade.hpp>

    // if compiling with non-cuda compiler use std tuple and boost's iterator_facade
    namespace tuple_ns = std;
    namespace iterator_facade_ns = boost;
#endif

#endif //CUDA_INCLUDE_H_
