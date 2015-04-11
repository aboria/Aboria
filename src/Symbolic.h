/*
 * Symbolic.h
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

#ifndef SYMBOLIC_H_
#define SYMBOLIC_H_

#include "DataVector.h"
#include "Vector.h"


#include <boost/mpl/bool.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/proto/core.hpp>
#include <boost/proto/context.hpp>
#include <boost/proto/traits.hpp>
#include <boost/proto/transform.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/ice.hpp>
#include <type_traits>
#include <tuple>
//#include <boost/fusion/include/vector.hpp>
//#include <boost/fusion/include/as_vector.hpp>
//#include <boost/fusion/include/joint_view.hpp>
//#include <boost/fusion/include/single_view.hpp>

namespace mpl = boost::mpl;
namespace proto = boost::proto;
using proto::_;
using proto::N;





namespace Aboria {

template<typename I>
struct label_ {};


namespace tag {
		struct sum_;
	}




template<typename Expr>
struct DataVectorExpr;

template<typename Expr>
struct LabelExpr;

struct null {
};

template<typename Value>
struct get_particles_ref : proto::callable
{

	typedef null result_type;

	result_type operator ()(const Value& arg)
	{
		return null();
	}
};

template<typename I, typename DataType>
struct get_particles_ref<DataVector<I, DataType> > : proto::callable
{

	typedef const DataType &result_type;

	result_type operator ()(const DataVector<I,DataType> &arr)
	{
		return arr.get_particles();
	}
};

template<typename DataType1, typename DataType2>
struct particles_ref_accumulate : proto::callable
{
	typedef DataType1 result_type;

	result_type operator ()(DataType1 existing, DataType2 to_add) {
		BOOST_MPL_ASSERT(( boost::is_same< DataType1,DataType2 > ));
		if (existing != to_add) {
			throw std::runtime_error("Expression not valid: conflicting Particles data structures");
		}
		return existing;
	}
};



template<typename ToAdd>
struct particles_ref_accumulate<ToAdd,null> : proto::callable
{
	typedef ToAdd result_type;

	result_type operator ()(ToAdd to_add, null existing) {
		return to_add;
	}
};

template<typename Existing>
struct particles_ref_accumulate<null,Existing> : proto::callable
{
	typedef Existing result_type;

	result_type operator ()(null to_add, Existing existing) {
		return existing;
	}
};

template<>
struct particles_ref_accumulate<null,null> : proto::callable
{
	typedef null result_type;

	result_type operator ()(null to_add, null existing) {
		return null();
	}
};



// A grammar which matches all the assignment operators,
// so we can easily disable them.
struct AssignOps
  : proto::switch_<struct AssignOpsCases>
{};

// Here are the cases used by the switch_ above.
struct AssignOpsCases
{
    template<typename Tag, int D = 0> struct case_  : proto::not_<_> {};

    template<int D> struct case_< proto::tag::plus_assign, D >         : _ {};
    template<int D> struct case_< proto::tag::minus_assign, D >        : _ {};
    template<int D> struct case_< proto::tag::multiplies_assign, D >   : _ {};
    template<int D> struct case_< proto::tag::divides_assign, D >      : _ {};
    template<int D> struct case_< proto::tag::modulus_assign, D >      : _ {};
    template<int D> struct case_< proto::tag::shift_left_assign, D >   : _ {};
    template<int D> struct case_< proto::tag::shift_right_assign, D >  : _ {};
    template<int D> struct case_< proto::tag::bitwise_and_assign, D >  : _ {};
    template<int D> struct case_< proto::tag::bitwise_or_assign, D >   : _ {};
    template<int D> struct case_< proto::tag::bitwise_xor_assign, D >  : _ {};
};


// This grammar describes which DataVector expressions
// are allowed;#include <tuple>
//struct DataVectorGrammar
//  : proto::or_<
//  	  proto::when<proto::terminal, get_particles_ref<proto::_value>(proto::_value)>
//  	  , proto::when<proto::sum_<DataVectorGrammar,DataVectorGrammar>, null()>
//  	  , proto::when<proto::and_<
//              	  	  proto::nary_expr<_, proto::vararg<DataVectorGrammar> >
//            		, proto::not_<AssignOps> >
//  	  	  	  	  	, proto::fold<_, void *(), particles_ref_accumulate<get_particles_ptrs,proto::_state>(get_particles_ptrs, proto::_state) > >
//    >
//{};

// This grammar describes which DataVector expressions
// are allowed;
struct DataVectorGrammar
  : proto::or_<
  	  proto::terminal<_>
  	  , proto::and_<
              	  	  proto::nary_expr<_, proto::vararg<DataVectorGrammar> >
            		, proto::not_<AssignOps>
  	  	  	  	   >
  >

{};

  struct LabelGrammar
    : proto::or_<
        proto::when< proto::terminal<Particles<_> >
                    , proto::_value >
        , proto::terminal<label_<_> >
        , proto::when< 
                    _ 
    				, LabelGrammar(proto::_right) >
      >
  {};



//struct DataVectorGrammar
//		: proto::or_<
//
//		  // DataVectorTerminals return their value
//		  proto::when< proto::terminal< DataVector<_,_> >
//				, fusion::single_view<proto::_value>(proto::_value) >
//
//		  // Any other terminals return nothing ...
//		  , proto::when< proto::terminal<_>
//				, fusion::nil() >
//
//		  , proto::when<proto::sum_<DataVectorGrammar,DataVectorGrammar>, fusion::nil()>
//
//		  // For any non-terminals, concat all children values
//		  , proto::when< proto::nary_expr<_, proto::vararg<_> >
//		  	  , proto::fold<_, fusion::nil()
//					, fusion::joint_view<DataVectorGrammer,boost::add_const<proto::_state> > (DataVectorGrammer, proto::_state)
//					>
//		>
//		{};

//struct ParticlesConsistentCtx
//: proto::callable_context< ParticlesConsistentCtx, proto::null_context >
//{
//	ParticlesConsistentCtx(void *particles) {boost/type_traits/ice.hpp
//		particles_ptrs.push_back(particles);
//	}
//
//	typedef void * result_type;
//	template<int I, typename DataType>
//	void operator ()(proto::tag::terminal, const DataVector<I,DataType> &arr)
//	{
//		DataType *this_particles_ptr = &(arr.get_particles());
//		for (void *i: particles_ptrs) {
//			if (this_particles_ptr != i) {
//				throw std::runtime_error("Expression not valid: Particles data structure " <<
//						this_particles_ptr << " not in list of valid pointers");
//			}
//		}
//	}
//
//	template<typename Expr>
//	void operator ()(tag::sum, Expr arr)
//	{
//		DataType *this_particles_ptr = &(arr.get_particles());
//		for (void *i: particles_ptrs) {
//			if (this_particles_ptr != i) {
//				throw std::runtime_error("Expression not valid: Particles data structure " <<
//						this_particles_ptr << " not in list of valid pointers");
//			}
//		}
//	}
//
//	void *particles_ptr;
//};

  	    struct dx_ {};


		// Here is an evaluation context that indexes into a lazy vector
		// expression, and combines the result.
		template<typename ParticlesType1, typename ParticlesType2>
		struct TwoParticleCtx
		{
			TwoParticleCtx(const Vect3d& dx, const typename ParticlesType1::value_type& particle1, const typename ParticlesType1::value_type& particle2)
			: dx_(dx),particle1_(particle1),particle2_(particle2)
			{}

			template<
			typename Expr
			// defaulted template parameters, so we can
			// specialize on the expressions that need
			// special handling.
			, typename Tag = typename proto::tag_of<Expr>::type
			, typename Arg0 = typename Expr::proto_child0
			>
			struct eval: proto::default_eval<Expr, TwoParticleCtx const>
			{};


			template<typename Expr, typename I, typename ParticlesTypeOther>
			struct eval<Expr, proto::tag::terminal, DataVector<I,ParticlesTypeOther> >
			{
				BOOST_MPL_ASSERT(( boost::type_traits::ice_or<
										boost::is_same< ParticlesTypeOther,ParticlesType1 >::value
									  , boost::is_same< ParticlesTypeOther,ParticlesType2 >::value
										>));
				typedef typename Elem<I::value,ParticlesTypeOther>::type result_type;

				result_type operator ()(Expr &expr, TwoParticleCtx const &ctx) const
				{
					return Elem<I::value,ParticlesTypeOther>::get(ctx.particle1_);
				}
			};


			// Handle dx terminals here...
			template<typename Expr>
			struct eval<Expr, proto::tag::terminal, dx_ >
			{
				typedef const Vect3d& result_type;

				result_type operator ()(Expr &expr, TwoParticleCtx const &ctx) const
				{
					return ctx.dx_;
				}
			};

			const typename ParticlesType1::value_type& particle1_;
			const typename ParticlesType2::value_type& particle2_;
			const Vect3d& dx_;

		};

// Here is an evaluation context that indexes into a lazy vector
// expression, and combines the result.
template<typename ParticlesType>
struct ParticleCtx
{
	ParticleCtx(const typename ParticlesType::value_type& particle)
      : particle_(particle)
    {}

	template<
	        typename Expr
	        // defaulted template parameters, so we can
	        // specialize on the expressions that need
	        // special handling.
	      , typename Tag = typename proto::tag_of<Expr>::type
	      , typename Arg0 = typename Expr::proto_child0
	    >
	struct eval: proto::default_eval<Expr, ParticleCtx const>
						{};



			// Handle vector terminals here...
			template<typename Expr, typename I, typename ParticlesType2>
			struct eval<Expr, proto::tag::terminal, DataVector<I,ParticlesType2> >
			{
				typedef typename Elem<I::value,ParticlesType2>::type result_type;

				BOOST_MPL_ASSERT(( boost::is_same< ParticlesType,ParticlesType2 > ));

				result_type operator ()(Expr &expr, ParticleCtx const &ctx) const
				{
					return Elem<I::value,ParticlesType2>::get(ctx.particle_);
				}
			};



			// Handle sums here...
			template<typename Expr, typename Arg0>
			struct eval<Expr, tag::sum_, Arg0 >
			{
				typedef typename proto::result_of::child_c<Expr,0>::type child0_type;
				typedef typename proto::result_of::child_c<Expr,1>::type child1_type;
				typedef typename proto::result_of::child_c<Expr,2>::type child2_type;

				//BOOST_MPL_ASSERT(( proto::matches< child0_type, LabelGrammar<1> > ));

				typedef typename boost::result_of<LabelGrammar(child0_type)>::type particles_type_ref;
                typedef typename std::remove_reference<particles_type_ref>::type particles_type;
				//typedef typename Expr::proto_child0::proto_child1::proto_value particles_type;
				typedef typename proto::result_of::eval<child1_type const, TwoParticleCtx<ParticlesType,particles_type> const>::type conditional_type;

				BOOST_MPL_ASSERT(( boost::is_same<conditional_type,bool > ));

				typedef typename proto::result_of::eval<child2_type const, TwoParticleCtx<ParticlesType,particles_type> const>::type result_type_const_ref;
				typedef typename std::remove_const<typename std::remove_reference<result_type_const_ref>::type>::type  result_type;


				result_type operator ()(Expr &expr, ParticleCtx const &ctx) const
				{
					particles_type_ref particlesb = LabelGrammar()(proto::child_c<0>(expr));
					child1_type conditional = proto::child_c<1>(expr);
					child2_type arg = proto::child_c<2>(expr);
					result_type sum = 0;
                    std::cout << "doing sum for particle "<<ctx.particle_.get_id()<<std::endl;
					for (auto i: particlesb.get_neighbours(ctx.particle_.get_position())) {
                        std::cout << "doing neighbour "<<std::get<0>(i).get_id()<<std::endl;
						TwoParticleCtx<ParticlesType,particles_type> ctx2(std::get<1>(i),std::get<0>(i),ctx.particle_);
						if (proto::eval(conditional,ctx2)) {
                            std::cout <<"conditional is true"<<std::endl;
							sum += proto::eval(arg,ctx2);
						} else {
                            std::cout <<"conditional is true"<<std::endl;
                        }
					}
					return sum;
				}
				};

	const typename ParticlesType::value_type& particle_;
};






// Tell proto how to generate expressions in the DataVectorDomain
struct DataVectorDomain
		: proto::domain<proto::generator<DataVectorExpr>, DataVectorGrammar >
		{};

// Declare that phoenix_domain is a sub-domain of spirit_domain
 struct LabelDomain 
        : proto::domain<proto::generator<LabelExpr>, LabelGrammar, DataVectorDomain>
        {};


	template< typename Expr >
	struct norm_fun
	{
		typedef double result_type;

		double operator()(const Vect3d& vector) const
		{
			return vector.norm();
		}
	};

	template<typename Expr>
	typename proto::result_of::make_expr<
	proto::tag::function  // Tag type
	, DataVectorDomain
	, norm_fun< Expr >        // First child (by value)
	, Expr const &
	>::type const
	norm_(Expr const &arg)
	{
		return proto::make_expr<proto::tag::function, DataVectorDomain>(
				norm_fun<Expr>()    // First child (by value)
				, boost::ref(arg)
		);
	}



	template<typename LABEL, typename CONDITIONAL, typename ARG>
	typename proto::result_of::make_expr<
	tag::sum_,
	DataVectorDomain,
	LABEL const &,
	CONDITIONAL const &,
	ARG const &
	>::type const
	sum_(LABEL const & label,CONDITIONAL const & conditional, ARG const & arg)
	{
		return proto::make_expr<tag::sum_, DataVectorDomain>(
				boost::ref(label),
				boost::ref(conditional),
				boost::ref(arg)
		);
		}

// Here is DataVectorExpr, which extends a proto expr type by
// giving it an operator [] which uses the ParticleCtx
// to evaluate an expression with a given index.
template<typename Expr>
struct DataVectorExpr
		: proto::extends<Expr, DataVectorExpr<Expr>, DataVectorDomain>
{
			explicit DataVectorExpr(Expr const &expr)
			: proto::extends<Expr, DataVectorExpr<Expr>, DataVectorDomain>(expr)
			  {}

			// Use the ParticleCtx to implement subscripting
			// of a DataVector expression tree.
			template<typename ParticleType>
			typename proto::result_of::eval<Expr const, ParticleCtx<ParticleType> const>::type
			eval( const typename ParticleType::value_type& particle) const
			{
				ParticleCtx<ParticleType> ctx(particle);
				return proto::eval(*this, ctx);
			}

			template<typename ParticleType1, typename ParticleType2>
			typename proto::result_of::eval<Expr const, TwoParticleCtx<ParticleType1,ParticleType2> const>::type
			eval( const Vect3d& dx, const typename ParticleType1::value_type& particle1,  const typename ParticleType2::value_type& particle2) const
			{
				TwoParticleCtx<ParticleType1,ParticleType2> ctx(dx, particle1, particle2);
				return proto::eval(*this, ctx);
			}
};

template<typename Expr>
struct LabelExpr: proto::extends<Expr, LabelExpr<Expr>, LabelDomain>
{
	explicit LabelExpr(Expr const &expr)
		: proto::extends<Expr, LabelExpr<Expr>, LabelDomain>(expr)
		{}


    BOOST_PROTO_EXTENDS_USING_ASSIGN(LabelExpr)
};


template<unsigned int I>
struct Label 
    : LabelExpr<typename proto::terminal<label_<mpl::int_<I> > >::type> {

	typedef typename proto::terminal<label_<mpl::int_<I> > >::type expr_type;
    typedef label_<mpl::int_<I> > data_type;

	explicit Label()
		: LabelExpr<expr_type>( expr_type::make(data_type()))
	  	{}


    BOOST_PROTO_EXTENDS_USING_ASSIGN(Label)
};




struct Dx
    : proto::terminal<dx_>::type {};


template<typename I, typename ParticlesType>
struct DataVectorSymbolic
	: DataVectorExpr<typename proto::terminal<DataVector<I,ParticlesType> >::type> {

	typedef typename proto::terminal<DataVector<I,ParticlesType> >::type expr_type;
	typedef typename ParticlesType::value_type particle_type;

	explicit DataVectorSymbolic(ParticlesType &p)
	: DataVectorExpr<expr_type>( expr_type::make(DataVector<I,ParticlesType>(p)) )
	  {}


	template< typename Expr >
	DataVectorSymbolic &operator =(Expr const & expr) {
        BOOST_MPL_ASSERT_NOT(( boost::is_same<I,mpl::int_<ID> > ));
		return this->assign(proto::as_expr<DataVectorDomain>(expr));
	}

private:

	template< typename Expr >
	DataVectorSymbolic &assign(Expr const & expr)
	{
		ParticlesType &particles = proto::value(*this).get_particles();
		std::for_each(particles.begin(),particles.end(),[&expr](particle_type& i) {
			Elem<I::value,ParticlesType>::set(i,expr.template eval<ParticlesType>(i));
		});
        if (I::value == POSITION) {
            particles.update_positions();
        }
		return *this;
	}
};



template<int I, typename ParticlesType>
DataVectorSymbolic<mpl::int_<I>,ParticlesType> get_vector(ParticlesType &p) {
	return DataVectorSymbolic<mpl::int_<I>,ParticlesType>(p);
}

}



#endif /* SYMBOLIC_H_ */
