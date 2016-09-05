
#ifndef ALGORITHMS_H_ 
#define ALGORITHMS_H_ 

#include "CudaInclude.h"
#include "Get.h"
#include "Traits.h"
#include <algorithm>



namespace Aboria {

template< class InputIt, class UnaryFunction >
#if defined(__CUDACC__)
InputIt for_each( InputIt first, InputIt last, UnaryFunction f ) {
    return thrust::for_each(first,last,f);
#else
UnaryFunction for_each( InputIt first, InputIt last, UnaryFunction f ) {
    return std::for_each(first,last,f);
#endif
}

}

#endif
