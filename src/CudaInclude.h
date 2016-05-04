#ifdef HAVE_THRUST
    #include <thrust/device_vector.h>
    #define CUDA_HOST_DEVICE __host__ __device__
    using namespace trust;
    using namespace trust::placeholders;
#else
    #include <boost/lambda/lambda.hpp>
    #include <boost/iterator/counting_iterator.hpp>
    #define CUDA_HOST_DEVICE 
    using namespace std;
    using namespace boost;
    using namespace boost::lambda;
#endif


