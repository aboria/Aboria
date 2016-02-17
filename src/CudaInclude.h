#ifdef HAVE_THRUST
    #include <thrust/device_vector.h>
    #define CUDA_HOST_DEVICE __host__ __device__
#else
    #define CUDA_HOST_DEVICE 
#endif

