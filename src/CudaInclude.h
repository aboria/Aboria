#ifndef CUDA_INCLUDE_H_
#define CUDA_INCLUDE_H_

#if defined(HAVE_THRUST)


#if defined(__CUDACC__)

#undef THRUST_DEVICE_SYSTEM
#define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_CUDA
#undef THRUST_HOST_SYSTEM
#define THRUST_HOST_SYSTEM THRUST_HOST_SYSTEM_CPP

#else

#undef THRUST_DEVICE_SYSTEM
#define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_OMP
#undef THRUST_HOST_SYSTEM
#define THRUST_HOST_SYSTEM THRUST_HOST_SYSTEM_OMP

#endif

#include <nppdefs.h>
#include <thrust/binary_search.h>
#include <thrust/device_vector.h>
#include <thrust/gather.h>
#include <thrust/host_vector.h>
#include <thrust/iterator/iterator_categories.h>
#include <thrust/iterator/iterator_facade.h>
#include <thrust/partition.h>
#include <thrust/random.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/transform_scan.h>
#include <thrust/tuple.h>
#include <thrust/unique.h>
#endif

#if defined(__CUDACC__)
#define CUDA_HOST_DEVICE __host__ __device__
#define CUDA_DEVICE __device__

#if not defined(__aboria_use_thrust_algorithms__)
#define __aboria_use_thrust_algorithms__
#endif

#define __aboria_hd_warning_disable__ #pragma hd_warning_disable

#define ABORIA_HOST_DEVICE_IGNORE_WARN #pragma hd_warning_disable

#else
#define CUDA_HOST_DEVICE
#define CUDA_DEVICE
#define __aboria_hd_warning_disable__
#define ABORIA_HOST_DEVICE_IGNORE_WARN

#include <boost/iterator/iterator_facade.hpp>

#endif

#endif // CUDA_INCLUDE_H_
