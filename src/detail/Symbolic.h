/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Aboria.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/



#ifndef SYMBOLIC_DETAIL_H_
#define SYMBOLIC_DETAIL_H_


#include <boost/mpl/bool.hpp>
#include <boost/mpl/assert.hpp>
#include "boost/mpl/and.hpp"
#include <boost/mpl/greater.hpp>
#include <boost/mpl/equal.hpp>
#include <boost/fusion/include/cons.hpp>
#include <boost/fusion/include/pair.hpp>
//#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/remove_if.hpp>
#include <boost/fusion/include/make_list.hpp>
#include <boost/fusion/include/make_map.hpp>
#include <boost/proto/core.hpp>
#include <boost/proto/context.hpp>
#include <boost/proto/traits.hpp>
#include <boost/proto/transform.hpp>
#include <boost/proto/debug.hpp>
#include <boost/utility/enable_if.hpp>
//#include <boost/type_traits/ice.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/result_of.hpp>
//#include <type_traits>
#include <tuple>
#include <map>

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;
namespace proto = boost::proto;
using proto::_;
using proto::N;


namespace Aboria {

    template<unsigned int I, typename P>
    struct Label; 

    namespace detail {

        template<typename Expr>
        struct SymbolicExpr;

        template<typename Expr>
        struct LabelExpr;

        template<typename Expr>
        struct GeometryExpr;

        // forward declare here so we can use the nice eval functions defined in Symbolic.h....
        template<typename labels_type=fusion::nil, typename dx_type=fusion::nil>
        struct EvalCtx;
    }
}

#include "detail/Terminal.h"
#include "detail/Grammars.h"



#include "detail/Domains.h"
#include "detail/Expressions.h"
#include "Vector.h"
#include "Get.h"
#include "Random.h"




#endif /* SYMBOLIC_H_ */
