#ifndef CUDA_INCLUDE_H_ 
#define CUDA_INCLUDE_H_ 

#ifdef HAVE_THRUST
    #include <thrust/device_vector.h>
    #include <thrust/sort.h>
    #include <thrust/binary_search.h>
    #include <thrust/sequence.h>
    #include <thrust/transform_scan.h>
    #include <nppdefs.h>
    #define CUDA_HOST_DEVICE __host__ __device__
#else
    #define CUDA_HOST_DEVICE 
#endif

#endif //CUDA_INCLUDE_H_
