
#ifndef GETTERTYPE_DETAIL_H_ 
#define GETTERTYPE_DETAIL_H_ 

#include "Helpers.h"
#include "Log.h"

namespace Aboria {
namespace detail {




/*
template<typename Tuple>
struct pointer_check {};

template<typename FirstType, typename ... Types>
struct pointer_check<std::tuple<FirstType,Types...>> {
    static const bool value = std::is_pointer<FirstType>::value;
};

#ifdef HAVE_THRUST 
template<typename FirstType, typename ... Types>
struct pointer_check<thrust::tuple<FirstType,Types...>> {
    static const bool value = std::is_pointer<FirstType>::value;
};
#endif
*/

}
}

#endif
