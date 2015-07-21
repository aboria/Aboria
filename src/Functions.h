/*
 * Functions.h
 * 
 * Copyright 2015 Martin Robinson
 *
 * This file is part of Aboria.
 *
 * Aboria is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Aboria is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Aboria.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 5 Feb 2015
 *      Author: robinsonm
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include "Symbolic.h"

namespace mpl = boost::mpl;
namespace proto = boost::proto;

namespace Aboria {



#define ABORIA_UNARY_FUNCTION(function_name, function_to_wrap,  domain) \
    template<typename Expr> \
	typename proto::result_of::make_expr<\
	proto::tag::function \
	, domain \
	, function_to_wrap     \
	, Expr const &    \
	>::type const   \
	function_name(Expr const &arg)  \
	{ \
		return proto::make_expr<proto::tag::function, domain>( \
				function_to_wrap()   \
				, boost::ref(arg)\
		);  \
	}\

#define ABORIA_BINARY_FUNCTION(function_name, function_to_wrap,  domain) \
    template<typename Expr1, typename Expr2> \
	typename proto::result_of::make_expr<\
	proto::tag::function \
	, domain \
	, function_to_wrap     \
	, Expr1 const &    \
	, Expr2 const &    \
	>::type const   \
	function_name(Expr1 const &arg1, Expr2 const &arg2)  \
	{ \
		return proto::make_expr<proto::tag::function, domain>( \
				function_to_wrap()   \
				, boost::ref(arg1) \
				, boost::ref(arg2) \
		);  \
	}\

#define ABORIA_TERNARY_FUNCTION(function_name, function_to_wrap,  domain) \
    template<typename Expr1, typename Expr2, typename Expr3> \
	typename proto::result_of::make_expr<\
	proto::tag::function \
	, domain \
	, function_to_wrap     \
	, Expr1 const &    \
	, Expr2 const &    \
	, Expr3 const &    \
	>::type const   \
	function_name(Expr1 const &arg1, Expr2 const &arg2, Expr3 const &arg3)  \
	{ \
		return proto::make_expr<proto::tag::function, domain>( \
				function_to_wrap()   \
				, boost::ref(arg1) \
				, boost::ref(arg2) \
				, boost::ref(arg3) \
		);  \
	}\

#define ABORIA_TAGGED_FUNCTION(name,tag) \
    template<typename LABEL, typename CONDITIONAL, typename ARG, typename INIT> \
	typename proto::result_of::make_expr< \
	tag , \
	DataVectorDomain, \
	LABEL const &, \
	CONDITIONAL const &, \
	ARG const &, \
	INIT const & \
	>::type const \
	name(LABEL const & label,CONDITIONAL const & conditional, ARG const & arg, INIT const & init) \
	{ \
		return proto::make_expr< tag , DataVectorDomain>( \
				boost::ref(label), \
				boost::ref(conditional), \
				boost::ref(arg), \
				boost::ref(init) \
		); \
		} \




	struct norm_fun
	{
		typedef double result_type;

		result_type operator()(const Vect3d& vector) const
		{
			return vector.norm();
		}
	};

    ABORIA_UNARY_FUNCTION(norm_, norm_fun, DataVectorDomain);


    template< typename T >
	struct reflect_fun
	{
		typedef Vect3d result_type;

		result_type operator()(const Vect3d& from, const Vect3d& to, const T geometry) const
		{
            Vect3d result = to;
            reflect_once(from,result,geometry);
			return result;
		}
	};

    ABORIA_TERNARY_FUNCTION(reflect_, reflect_fun<Expr3>, DataVectorDomain);


}
#endif /* FUNCTIONS_H_ */
