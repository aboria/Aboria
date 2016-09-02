#ifndef CUDA_INCLUDE_H_ 
#define CUDA_INCLUDE_H_ 

#if defined(__CUDACC__)
    #include <thrust/device_vector.h>
    #include <thrust/host_vector.h>
    #include <thrust/sort.h>
    #include <thrust/binary_search.h>
    #include <thrust/iterator/iterator_facade.h>
    #include <thrust/sequence.h>
    #include <thrust/transform_scan.h>
    #include <nppdefs.h>
    #define CUDA_HOST_DEVICE __host__ __device__

    #define __aboria_hd_warning_disable__ \
    #pragma hd_warning_disable

    namespace tuple_ns = thrust;
#else
    #define CUDA_HOST_DEVICE 
    #define __aboria_hd_warning_disable__

    namespace tuple_ns = std;
#endif

#endif //CUDA_INCLUDE_H_
