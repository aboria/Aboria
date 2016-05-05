#ifdef HAVE_THRUST
    #include <thrust/device_vector.h>
    #define CUDA_HOST_DEVICE __host__ __device__
    namespace astd = thrust;
    namespace extra_iterators = thrust;
    namespace lambda = thrust::placeholders;
#else
    #include <boost/lambda/lambda.hpp>
    #include <boost/iterator/counting_iterator.hpp>
    #define CUDA_HOST_DEVICE 
    namespace astd = std;
    namespace extra_iterators = boost;
    namespace lambda = boost::lambda;
#endif


