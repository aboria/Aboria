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

#include "Vector.h"
#include "Get.h"
#include "Random.h"


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
#include <boost/type_traits/ice.hpp>
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


        ////////////////////////
        /// Terminal Classes ///
        ////////////////////////

        template<typename I, typename P>
        struct label {
            typedef P particles_type;
            typedef I depth;

            label(P& p):p(p) {}
            P& get_particles() const { return p; }

            P& p;
        };


        template<typename T>
        struct symbolic {
            typedef T variable_type;
            typedef typename T::value_type value_type;

            std::vector<value_type>& get_buffer(const void *address) { return map[address]; }

            std::map<const void *,std::vector<value_type> > map;
        };

        /*
        template<typename I>
        struct unknown {
            typedef typename I::value value;
        };
        */

        template <typename T>
        struct accumulate {
            typedef T functor_type;
            typedef typename T::result_type init_type;
            accumulate():init(0) {};
            accumulate(const T& functor):functor(functor),init(0) {};
            void set_init(const init_type& arg) {
                init = arg;
            }
            T functor;
            init_type init;
        };

        template <typename T,unsigned int N>
        struct vector {
            typedef Vector<T,N> result_type;

            template <typename ...Types>
            result_type operator()(const Types... args) const {
                return result_type(args...);
            }
        };

        template <typename L1, typename L2>
        struct dx {
            typedef L1 label_a_type;
            typedef L2 label_b_type;
            label_a_type la;
            label_b_type lb;
            dx(label_a_type &la, label_b_type &lb):la(la),lb(lb) {};
            const label_a_type & get_label_a() const { return la; }
            const label_b_type & get_label_b() const { return lb; }
        };

        struct normal {
            //typedef std::mt19937 generator_type;
            normal() {};
            normal(uint32_t seed): generator(seed) {};
            double operator()() {
                std::normal_distribution<double> normal_distribution;
                return normal_distribution(generator);
            }
            /*
            double operator()(const size_t& id) const {
                std::normal_distribution<double> normal_distribution;
                generator_type gen(id);
                return normal_distribution(gen);
            }
            */
            double operator()(generator_type& gen) const {
                std::normal_distribution<double> normal_distribution;
                return normal_distribution(gen);
            }
            generator_type generator;

        };

        struct uniform {
            //typedef std::mt19937 generator_type;
            uniform() {};
            uniform(uint32_t seed): generator(seed) {};
            double operator()() {
                std::uniform_real_distribution<double> uniform_distribution;
                return uniform_distribution(generator);
            }
            /*
            double operator()(const size_t& id) const {
                std::uniform_real_distribution<double> uniform_distribution;
                generator_type gen(id);
                return uniform_distribution(gen);
            }
            */
            double operator()(generator_type& gen) const {
                std::uniform_real_distribution<double> normal_distribution;
                return normal_distribution(gen);
            }
            generator_type generator;
        };

        template <typename T>
        struct geometry {
            typedef T result_type;

            template <typename... A>
            result_type operator()(const A&... args) const {
                return T(args...);
            }
        };


        template <typename T>
        struct geometries {
            typedef T result_type;

            //result_type operator()(const Vect3d& centre, const double radius, const bool in) const
            template <typename... A>
            result_type operator()(const A&... args) const {
                return T(args...);
            }
        };




        ////////////////
        /// Grammars ///
        ////////////////


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


        struct SymbolicGrammar
            : proto::or_<
                proto::terminal<_>
                , proto::and_<
                    proto::nary_expr<_, proto::vararg<SymbolicGrammar> >
                    , proto::not_<AssignOps>
                >
              >

        {};

        struct get_particles: proto::callable {
            template<typename Sig>
                struct result;

            template<typename This, typename LabelType>
                struct result<This(const LabelType&)>
                : boost::add_reference<typename LabelType::particles_type>
                {};

            /*
            template<typename This, typename LabelType>
                struct result<This(LabelType)>
                : boost::add_reference<typename LabelType::particles_type>
                {};
                */

            template<typename LabelType>
                typename LabelType::particles_type &operator()(const LabelType& label) {
                    return label.get_particles();
                }
        };



        struct LabelGrammar
            : proto::when< proto::terminal<label<_,_> >
              , get_particles(proto::_value) >
        {};

        struct SubscriptGrammar
            : proto::when< proto::terminal<label<_,_> >
              , proto::_value >
        {};


        struct GeometryGrammar
            : proto::function< proto::terminal< geometry<_> >, SymbolicGrammar, SymbolicGrammar, SymbolicGrammar>
        {}; 

        struct GeometriesGrammar
            : proto::function< proto::terminal< geometries<_> >, LabelGrammar, SymbolicGrammar>
        {}; 


        //TODO: should make a "conditional" grammar to match 2nd arguement
        struct AccumulateGrammar
            : proto::function< proto::terminal< accumulate<_> >, LabelGrammar, SymbolicGrammar, SymbolicGrammar>
        {}; 


        struct remove_label: proto::callable {
            template<typename Sig>
            struct result;

            template<typename This, typename T1, typename T2>
            struct result<This(T1, T2)>:
                fusion::result_of::as_list<
                    typename fusion::result_of::remove_if<typename std::remove_reference<T2>::type,
                                boost::is_same<
                                    //typename std::remove_reference<typename std::remove_const<T1>::type>::type &,
                                    T1,
                                    mpl::_>
                                >::type
                > 
            {};

            /*
            template<typename This, typename T1, typename T2>
            struct result<This(T1, const T2&)>:
                fusion::result_of::as_list<
                    typename fusion::result_of::remove_if<T2,boost::is_same<T1,mpl::_>>::type
                > 
            {};
            */

            template<typename T1, typename T2>
            typename result<remove_label(const T1&, const T2&)>::type
            operator()(const T1& label, const T2& state) {
                return fusion::as_list(fusion::remove_if<boost::is_same<const T1&,mpl::_>>(state));
            }
        };

        struct push_back_if_new: proto::callable {
            template<typename Sig>
            struct result;

            /*
            template<typename This, typename T1, typename T2>
            struct result<This(const T1&, const T2&)> {
                typedef fusion::cons<const T1&, 
                        typename boost::result_of<remove_label(const T1&,const T2&)>::type> type;
            };
            */

            template<typename This, typename T1, typename T2>
            struct result<This(T1&, T2)> {
                typedef fusion::cons<const T1&, 
                        typename boost::result_of<remove_label(const T1&,const T2&)>::type> type;
            };

            /*
            template<typename This, typename T1, typename T2>
            struct result<This(T1&, const T2&)> {
                typedef fusion::cons<const T1&, 
                        typename boost::result_of<remove_label(const T1&,const T2&)>::type> type;
            };
            */

            /*
            template<typename This, typename T1, typename T2>
            struct result<This(T1&, const T2&)>:
                fusion::result_of::as_list<
                    typename fusion::result_of::push_back<
                        typename boost::result_of<remove_label(T1&,const T2&)>::type, T1>::type 
                > 
            {};
            */

            template<typename T1, typename T2>
            typename result<push_back_if_new(const T1&, const T2&)>::type
            operator()(const T1& label, const T2& state) {
                return typename result<push_back_if_new(const T1&, const T2&)>::type(label,remove_label()(label,state));
            }
        };


        struct get_a_from_dx: proto::callable {
            template<typename Sig>
                struct result;

            template<typename This, typename DxType>
            struct result<This(DxType&)> {
                typedef const typename DxType::label_a_type& type;
            };

            template<typename DxType>
            const typename DxType::label_a_type & operator()(const DxType& dx) {
                return dx.get_label_a();
            }
        };

        struct get_b_from_dx: proto::callable {
            template<typename Sig>
                struct result;

            template<typename This, typename DxType>
            struct result<This(DxType&)> {
                typedef const typename DxType::label_b_type& type;
            };

            template<typename DxType>
            const typename DxType::label_b_type & operator()(const DxType& dx) {
                return dx.get_label_b();
            }
        };

        struct get_labels: 
            proto::or_<
                
                // Put all label terminals at the head of the
                // list that we're building in the "state" parameter
                proto::when<
                    proto::terminal< label<_,_> >
                  , push_back_if_new(proto::_value,proto::_state)
                >
                , proto::when< 
                    proto::terminal< dx<_,_> >
                    , push_back_if_new(get_a_from_dx(proto::_value),
                            push_back_if_new(get_b_from_dx(proto::_value),proto::_state))
                >
                //don't add other terminals
                , proto::when<
                    proto::terminal<_>
                  , proto::_state
                >
                , proto::when<
                    proto::function< proto::terminal< accumulate<_> >, _, _, _>
                    , remove_label(proto::_value(proto::_child1), 
                                    get_labels(proto::_child2,
                                            get_labels(proto::_child3,proto::_state)))
                >
                , proto::otherwise< 
                    proto::fold<_, proto::_state, get_labels> 
                >
           >
        {};

    }
	struct norm_fun;
    namespace detail {

        struct norm_dx:
            proto::function< proto::terminal<Aboria::norm_fun>, proto::terminal<dx<_,_>>> 
        {};

        struct range_if_expr:
            proto::or_<
                proto::less<norm_dx,SymbolicGrammar>
                , proto::less_equal<norm_dx,SymbolicGrammar>
                , proto::greater<SymbolicGrammar,norm_dx>
                , proto::greater_equal<SymbolicGrammar,norm_dx>
                , proto::nary_expr<_, proto::vararg<range_if_expr> >
            >
        {};

        namespace result_of {
            template <typename Expr>
            struct get_labels {
                typedef typename boost::result_of<Aboria::detail::get_labels(Expr,fusion::nil_)>::type get_labels_result;
                typedef typename std::remove_reference<get_labels_result>::type type;
            };

        }


        template<typename Expr>
        struct is_const: 
            mpl::equal_to<
                typename fusion::result_of::size<
                    typename result_of::get_labels<Expr>::type
                >::type
                , mpl::int_<0>
            >
        {};

        template<typename Expr>
        struct is_univariate: 
            mpl::equal_to<
                typename fusion::result_of::size<
                    typename result_of::get_labels<Expr>::type
                >::type
                , mpl::int_<1>
            >
        {};

        template<typename Expr>
        struct is_bivariate: 
            mpl::equal_to<
                typename fusion::result_of::size<
                    typename result_of::get_labels<Expr>::type
                >::type
                , mpl::int_<2>
            >
        {};

        template<typename Expr,unsigned int I>
        struct get_label_c {
            typedef typename result_of::get_labels<Expr>::type labels_type;
            typedef typename fusion::result_of::value_at_c<get_labels,I>::type type;
        };
            

        ////////////////
        /// Contexts ///
        ////////////////

        // Here is an evaluation context that indexes into a lazy vector
        // expression, and combines the result.
        template<typename labels_type=fusion::nil_, typename dx_type=fusion::nil_>
        struct EvalCtx {
            typedef typename fusion::result_of::size<labels_type>::type size_type;
            typedef typename fusion::result_of::size<dx_type>::type dx_size_type;
            static constexpr int dx_size = size_type::value*(size_type::value-1)/2;

            //BOOST_MPL_ASSERT_MSG(dx_size_type::value==dx_size,DX_SIZE_NOT_CONSISTENT_WITH_LABELS_SIZE,(dx_size,dx_size_type));
            static_assert(dx_size_type::value==dx_size,"dx size not consitent with labels_size");
            
            EvalCtx(labels_type labels=fusion::nil_(), dx_type dx=fusion::nil())
                : m_labels(labels),m_dx(dx)
            {}

            template<
                typename Expr
                // defaulted template parameters, so we can
                // specialize on the expressions that need
                // special handling.
                , typename Tag = typename proto::tag_of<Expr>::type
                , typename Enable = void
                >
            struct eval: proto::default_eval<Expr, EvalCtx const>
            {};


            // Handle normal and uniform terminals here...
            template<typename Expr>
            struct eval<Expr, proto::tag::terminal,
            typename boost::enable_if<
                mpl::and_<
                    mpl::or_<
                        proto::matches<Expr, proto::terminal<normal>>,
                        proto::matches<Expr, proto::terminal<uniform>>
                            >,
                    mpl::greater<size_type,mpl::int_<0>>
                >>::type> {

                typedef double result_type;

                result_type operator ()(Expr &expr, EvalCtx const &ctx) const
                {
                    //TODO: get better (parallel) random number generator
                    //return const_cast<typename ParticlesType::value_type&>(ctx.particle_).rand_normal();
                    return proto::value(expr)(
                            const_cast<generator_type&>(get<random>(fusion::front(ctx.m_labels).second)));
                            //get<id>(fusion::front(ctx.m_labels).second));
                }
            };


            // Handle subscripts here...
            template<typename Expr>
            struct eval<Expr, proto::tag::subscript,
            typename boost::enable_if<
                mpl::and_<
                    proto::matches<typename proto::result_of::child_c<Expr,1>::type,
                        proto::terminal<label<_,_>>>, 
                    mpl::greater<size_type,mpl::int_<0>>
                >>::type> {

                typedef typename proto::result_of::child_c<Expr, 1>::type child1_type;
                typedef typename proto::result_of::value<child1_type>::type label_type;

                /*
                typedef typename label_type::particles_type particles_type;
                typedef typename particles_type::const_reference particle_ref;
                typedef typename fusion::pair<label_type,particle_ref> search_type;

                */
                static_assert(fusion::result_of::has_key<labels_type,label_type>::value,
                        "label not in evaluation context");

                typedef typename proto::result_of::child_c<Expr,0>::type child0_type;
                typedef typename proto::result_of::value<child0_type>::type symbolic_type;
                typedef typename symbolic_type::variable_type variable_type; 
                typedef const typename variable_type::value_type & result_type;

                result_type operator ()(Expr &expr, EvalCtx const &ctx) const {
                    return get<variable_type>(fusion::at_key<label_type>(ctx.m_labels));
                }
            };

            // Handle dx terminals here...
            template<typename Expr>
            struct eval<Expr, proto::tag::terminal, 
            typename boost::enable_if<
                mpl::and_<
                    proto::matches<Expr, proto::terminal<dx<_,_>>>,
                    mpl::equal<size_type,mpl::int_<2>>
            >>::type> {
                typedef typename fusion::result_of::front<const dx_type>::type result_type;
                typedef typename proto::result_of::value<Expr>::type expr_dx;
                typedef typename expr_dx::label_a_type expr_label_a_type;
                typedef typename expr_dx::label_b_type expr_label_b_type;

                /*
                BOOST_MPL_ASSERT_MSG((fusion::result_of::has_key<labels_type,expr_label_a_type>::value),ASDFASDFASDF,(expr_dx,labels_type,expr_label_a_type));
                BOOST_MPL_ASSERT_MSG((fusion::result_of::has_key<labels_type,expr_label_b_type>::value),ASDFASDFASDF,(expr_dx));
                */

                static_assert(fusion::result_of::has_key<labels_type,expr_label_a_type>::value,
                        "dx label a not in evaluation context");
                static_assert(fusion::result_of::has_key<labels_type,expr_label_b_type>::value,
                        "dx label b not in evaluation context");

                result_type operator ()(Expr &expr, EvalCtx const &ctx) const {
                    return fusion::front(ctx.m_dx);
                }
            };

            template <typename result_type,
                     typename label_type,
                     typename if_expr_type, 
                     typename expr_type,
                     typename accumulate_type>
            static
            typename boost::enable_if<
                mpl::not_<proto::matches<if_expr_type,range_if_expr>>
            ,result_type>::type
            sum_impl(label_type label, 
                    if_expr_type& if_expr, 
                    expr_type& expr, 
                    accumulate_type& accum,
                    const EvalCtx& ctx,mpl::int_<0>) { //note: using tag dispatching here cause I couldn't figure out how to do this via enable_if....

                result_type sum = accum.init;
                for (auto i: label.get_particles()) {
                    auto new_labels = fusion::make_map<label_type>(i);
                    EvalCtx<decltype(new_labels),decltype(ctx.m_dx)> const new_ctx(new_labels,ctx.m_dx);

                    if (proto::eval(if_expr,new_ctx)) {
                        sum = accum.functor(sum,proto::eval(expr,new_ctx));
                    }
                }
                return sum;
            }

            template <typename result_type,
                     typename label_b_type,
                     typename if_expr_type, 
                     typename expr_type,
                     typename accumulate_type>
            static
            typename boost::enable_if<
                    mpl::not_<proto::matches<if_expr_type,range_if_expr>>
            ,result_type>::type
            sum_impl(label_b_type label, 
                    if_expr_type& if_expr, 
                    expr_type& expr, 
                    accumulate_type& accum,
                    const EvalCtx& ctx, mpl::int_<1>) {

                typedef typename label_b_type::particles_type particles_b_type;

                result_type sum = accum.init;
                const particles_b_type& particlesb = label.get_particles(); //copies particles!!!!
                auto ai = fusion::front(ctx.m_labels).second;

                typedef typename particles_b_type::position position;
                typedef typename position::value_type double_d;
                typedef typename std::remove_reference<
                    typename fusion::result_of::at_c<labels_type,0>::type>::type::first_type label_a_type;



                ASSERT(!particlesb.get_periodic().any(),"periodic does not work with dense");
                const size_t nb = particlesb.size();
                for (size_t i=0; i<nb; ++i) {
                    typename particles_b_type::const_reference bi = particlesb[i];

                    auto new_labels = fusion::make_map<label_a_type,label_b_type>(
                                        ai,bi);

                    const double_d dx = get<position>(bi)-get<position>(ai);


                    EvalCtx<decltype(new_labels),fusion::list<const double_d &>> const new_ctx(
                            new_labels,
                            fusion::make_list(boost::cref(dx)));

                    if (proto::eval(if_expr,new_ctx)) {
                        sum = accum.functor(sum,proto::eval(expr,new_ctx));
                    }
                }
                return sum;
            }

            template <typename result_type,
                     typename label_b_type,
                     typename if_expr_type, 
                     typename expr_type,
                     typename accumulate_type,typename dummy=size_type>
            static
            typename boost::enable_if<
                    proto::matches<if_expr_type,range_if_expr>
            ,result_type>::type
            sum_impl(label_b_type label,
                    if_expr_type& if_expr, 
                    expr_type& expr, 
                    accumulate_type& accum,
                    const EvalCtx& ctx, mpl::int_<1>) {


                typedef typename label_b_type::particles_type particles_b_type;
                result_type sum = accum.init;
                const particles_b_type& particlesb = label.get_particles(); //copies particles!!!!
                auto ai = fusion::front(ctx.m_labels).second;

                typedef typename std::remove_reference<
                    typename fusion::result_of::at_c<labels_type,0>::type>::type::first_type label_a_type;
                typedef typename label_a_type::particles_type particles_a_type;
                typedef typename particles_a_type::position position;

                for (auto i: particlesb.get_neighbours(get<position>(ai))) {
                    auto new_labels = fusion::make_map<label_a_type,label_b_type>(
                                        ai,std::get<0>(i));

                    auto new_dx = fusion::make_list(boost::cref(std::get<1>(i)));

                    EvalCtx<decltype(new_labels),decltype(new_dx)> const new_ctx(new_labels,new_dx);

                    if (proto::eval(if_expr,new_ctx)) {
                        sum = accum.functor(sum,proto::eval(expr,new_ctx));
                        
                    }
                }
                return sum;
            }

            template<typename Expr>
            struct eval<Expr, proto::tag::function,
            typename boost::enable_if<
                    proto::matches<Expr,AccumulateGrammar>
                >::type> {

                typedef typename proto::result_of::child_c<Expr,0>::type child0_type;
                typedef typename proto::result_of::value<child0_type>::type functor_terminal_type;
                typedef typename functor_terminal_type::functor_type functor_type;
                typedef typename functor_type::result_type result_type;

                result_type operator ()(Expr &expr, EvalCtx const &ctx) const {
                    return sum_impl<result_type>(proto::value(proto::child_c<1>(expr)),
                            proto::child_c<2>(expr),
                            proto::child_c<3>(expr),
                            proto::value(proto::child_c<0>(expr)),ctx,size_type());
                }
            };

            labels_type m_labels;
            dx_type m_dx;
    };

        ////////////////
        /// Domains  ///
        ////////////////



        // Tell proto how to generate expressions in the SymbolicDomain
        struct SymbolicDomain
            : proto::domain<proto::generator<SymbolicExpr>, SymbolicGrammar > {};

        // Declare that phoenix_domain is a sub-domain of spirit_domain
        struct LabelDomain 
            : proto::domain<proto::generator<LabelExpr>, LabelGrammar, SymbolicDomain> {};

        struct GeometryDomain 
            : proto::domain<proto::generator<GeometryExpr>, GeometryGrammar> {};


        ///////////////////
        /// Expressions ///
        ///////////////////

        template<typename Expr>
        struct SymbolicExpr
            : proto::extends<Expr, SymbolicExpr<Expr>, SymbolicDomain> {

            explicit SymbolicExpr(Expr const &expr)
                : proto::extends<Expr, SymbolicExpr<Expr>, SymbolicDomain>(expr)
            {}
        };

        template<typename Expr, typename Enable=void>
        struct symbolic_helper {};

        template<typename Expr>
        struct symbolic_helper<Expr, typename boost::enable_if<is_const<
                    typename proto::result_of::as_expr<
                                        Expr,detail::SymbolicDomain>::type
                    >>::type> {

            typedef typename proto::result_of::as_expr<
                            Expr,detail::SymbolicDomain>::type expr_type;

            typedef typename std::remove_const<typename proto::result_of::deep_copy<expr_type>::type>::type deep_copy_type;

            typedef EvalCtx<> const_context_type;
            typedef typename proto::result_of::eval<expr_type, const_context_type const>::type result;

            typedef typename std::remove_cv<typename std::remove_reference<
                        result
                        >::type>::type result_base_type;
     
        };

        template<typename Expr>
        struct symbolic_helper<Expr, typename boost::enable_if<is_univariate<
                    typename proto::result_of::as_expr<
                                        Expr,detail::SymbolicDomain>::type
                    >>::type> {

            typedef typename proto::result_of::as_expr<
                            Expr,detail::SymbolicDomain>::type expr_type;

            typedef typename std::remove_const<typename proto::result_of::deep_copy<expr_type>::type>::type deep_copy_type;
 
            typedef typename fusion::result_of::at_c<typename std::result_of<get_labels(expr_type,fusion::nil_)>::type,0>::type label_a_type_ref;
            typedef typename std::remove_const<
                typename std::remove_reference<label_a_type_ref>::type>::type label_a_type;
            typedef typename label_a_type::particles_type particles_a_type;
            typedef typename particles_a_type::double_d double_d;
            typedef typename particles_a_type::const_reference particle_a_reference;
            typedef EvalCtx<fusion::map<fusion::pair<label_a_type,particle_a_reference>>> univariate_context_type;
            typedef typename proto::result_of::eval<expr_type, univariate_context_type const>::type result;

            typedef typename std::remove_cv<typename std::remove_reference<
                        result
                        >::type>::type result_base_type;
            

        };

        template<typename Expr>
        struct symbolic_helper<Expr, typename boost::enable_if<is_bivariate<
                    typename proto::result_of::as_expr<
                                        Expr,detail::SymbolicDomain>::type
                    >>::type> {
            typedef typename proto::result_of::as_expr<
                            Expr,detail::SymbolicDomain>::type expr_type;

            typedef typename std::remove_const<typename proto::result_of::deep_copy<expr_type>::type>::type deep_copy_type;

            typedef typename fusion::result_of::at_c<typename std::result_of<get_labels(expr_type,fusion::nil_)>::type,0>::type label_first_type_ref;
            typedef typename fusion::result_of::at_c<typename std::result_of<get_labels(expr_type,fusion::nil_)>::type,1>::type label_second_type_ref;
            typedef typename std::remove_const<
                typename std::remove_reference<label_first_type_ref>::type>::type label_first_type;
            typedef typename std::remove_const<
                typename std::remove_reference<label_second_type_ref>::type>::type label_second_type;

            static_assert(mpl::not_equal_to<typename label_first_type::depth
                                           ,typename label_second_type::depth>::value,
                                         "label a depth equal to label b");

            typedef typename mpl::if_<
                                mpl::less<typename label_first_type::depth
                                         ,typename label_second_type::depth>
                                    ,label_first_type
                                    ,label_second_type>::type label_a_type;
            typedef typename mpl::if_<
                                mpl::greater<typename label_first_type::depth
                                            ,typename label_second_type::depth>
                                    ,label_first_type
                                    ,label_second_type>::type label_b_type;

            typedef typename label_a_type::particles_type particles_a_type;
            typedef typename label_b_type::particles_type particles_b_type;
            typedef typename particles_a_type::const_reference particle_a_reference;
            typedef typename particles_b_type::const_reference particle_b_reference;
            typedef typename particles_a_type::double_d double_d;
            typedef EvalCtx<fusion::map<
                fusion::pair<label_a_type,particle_a_reference>,
                fusion::pair<label_b_type,particle_b_reference>>,
                fusion::list<const double_d &>> bivariate_context_type;
            typedef typename proto::result_of::eval<expr_type, bivariate_context_type const>::type result;

            typedef typename std::remove_cv<typename std::remove_reference<
                        result
                        >::type>::type result_base_type;
             
        };


        /*
        template<typename Expr>
        struct SymbolicExpr<Expr,typename boost::enable_if<proto::matches<Expr, detail::univariate_expr> >::type> 
            : proto::extends<Expr, SymbolicExpr<Expr>, SymbolicDomain> {

            typedef typename std::result_of<detail::univariate_expr(Expr)>::type label_a_type_ref;
            typedef typename std::remove_reference<label_a_type_ref>::type label_a_type;
            typedef typename label_a_type::particles_type particles_a_type;

            explicit SymbolicExpr(Expr const &expr)
                : proto::extends<Expr, SymbolicExpr<Expr>, SymbolicDomain>(expr)
            {}

            template <typename unknown_tuple_type=std::tuple<>>
            typename proto::result_of::eval<Expr const, ParticleCtx<particles_a_type,unknown_tuple_type> const>::type
            eval( typename particles_a_type::const_reference particle, unknown_tuple_type unknown_tuple = std::tuple<>()) const {
                ParticleCtx<particles_a_type,unknown_tuple_type> const ctx(particle,unknown_tuple);
                return proto::eval(*this, ctx);
            }


        };

        template<typename Expr>
        struct SymbolicExpr<Expr,typename boost::enable_if<proto::matches<Expr, detail::bivariate_expr> >::type> 
            : proto::extends<Expr, SymbolicExpr<Expr>, SymbolicDomain> {

            typedef typename std::result_of<detail::bivariate_expr(Expr)>::type::first label_a_type_ref;
            typedef typename std::result_of<detail::bivariate_expr(Expr)>::type::second label_b_type_ref;
            typedef typename std::remove_reference<label_a_type_ref>::type label_a_type;
            typedef typename std::remove_reference<label_b_type_ref>::type label_b_type;
            typedef typename label_a_type::particles_type particles_a_type;
            typedef typename label_b_type::particles_type particles_b_type;
            typedef typename particles_a_type::double_d double_d;

            explicit SymbolicExpr(Expr const &expr)
                : proto::extends<Expr, SymbolicExpr<Expr>, SymbolicDomain>(expr)
            {}

            template <typename unknown_tuple_type=std::tuple<>>
            typename proto::result_of::eval<Expr const, TwoParticleCtx<particles_a_type,particles_b_type,unknown_tuple_type> const>::type
            eval( const double_d& dx, typename particles_a_type::const_reference particle1, typename particles_b_type::const_reference particle2,unknown_tuple_type unknown_tuple1 = std::tuple<>(), unknown_tuple_type unknown_tuple2 = std::tuple<>()) const {
                TwoParticleCtx<particles_a_type,particles_b_type, unknown_tuple_type> const ctx(dx,particle1,particle2,unknown_tuple1,unknown_tuple2);
                return proto::eval(*this, ctx);
            }

        }; 
        */



        template<typename Expr>
        struct GeometryExpr: proto::extends<Expr, GeometryExpr<Expr>, GeometryDomain> {
            explicit GeometryExpr(Expr const &expr)
                : proto::extends<Expr, GeometryExpr<Expr>, GeometryDomain>(expr)
            {}
        };


        template<typename Expr>
        struct LabelExpr: proto::extends<Expr, LabelExpr<Expr>, LabelDomain> {
            explicit LabelExpr(Expr const &expr)
                : proto::extends<Expr, LabelExpr<Expr>, LabelDomain>(expr)
            {}


            BOOST_PROTO_EXTENDS_USING_ASSIGN(LabelExpr)
        };


        template<unsigned int I, typename P>
        struct Label 
            : LabelExpr<typename proto::terminal<label<mpl::int_<I>,P > >::type> {

                typedef typename proto::terminal<label<mpl::int_<I>,P > >::type expr_type;
                typedef label<mpl::int_<I>,P > data_type;

                explicit Label(P& p)
                    : LabelExpr<expr_type>( expr_type::make(data_type(p)))
                {}


                //BOOST_PROTO_EXTENDS_USING_ASSIGN(Label)
            };



        #define SUBSCRIPT_TYPE \
        boost::proto::exprns_::basic_expr< \
            boost::proto::tagns_::tag::subscript \
            , boost::proto::argsns_::list2< \
                Aboria::detail::SymbolicExpr<boost::proto::exprns_::expr< \
                    boost::proto::tagns_::tag::terminal \
                    , boost::proto::argsns_::term<Aboria::detail::symbolic<VariableType> > \
                , 0l> >& \
            , Aboria::Label<0u, ParticlesType >& \
            > \
        , 2l \
        >


        // TODO: this seems a messy way to define a symbol subscripted by a label. might
        // be better to put a subscript operator in the Symbol class?
        // for: cleaner
        // against: messes up the grammer (already existing subscript expression). 
        // Could have Symbol *not* be an expression, but then can't assume labels 
        // within expressions anymore....
        template<typename VariableType, typename ParticlesType>
        struct SymbolicExpr<SUBSCRIPT_TYPE>
            : proto::extends<SUBSCRIPT_TYPE, SymbolicExpr<SUBSCRIPT_TYPE>, SymbolicDomain> {
                typedef SUBSCRIPT_TYPE Expr;

                typedef typename ParticlesType::position position;
                

                typedef typename std::remove_const< 
                    typename std::remove_reference<
                    typename proto::result_of::value< 
                    typename proto::result_of::child_c<Expr,1>::type 
                    >::type>::type>::type label_type;
                typedef typename proto::result_of::value< 
                    typename proto::result_of::child_c<Expr,0>::type 
                    >::type symbol_type;

                struct is_not_my_label
                    : proto::and_<
                        proto::terminal<label<_,_>>
                        , proto::if_< mpl::not_<
                                boost::is_same< proto::_value, label_type>>() >
                      >
                {};

                struct is_my_symbol
                    : proto::terminal<symbol_type>
                {};

                struct is_not_aliased
                    : proto::or_<
                        proto::and_<
                            proto::terminal<proto::_>
                            ,proto::not_< 
                                proto::and_<
                                    proto::if_<boost::is_same<symbol_type,
                                                       symbolic<position>>()>
                                    , proto::terminal<dx<_,_>>
                                >
                             >
                          >
                        , proto::and_<
                            proto::nary_expr< proto::_, proto::vararg<is_not_aliased>>
                            ,proto::not_<proto::subscript<is_my_symbol,is_not_my_label>>
                          >
                      >
                {};
                

        #undef SUBSCRIPT_TYPE

                explicit SymbolicExpr(Expr const &expr)
                    : proto::extends<Expr, SymbolicExpr<Expr>, SymbolicDomain>(expr),
                    msymbol(proto::value(proto::child_c<0>(expr))),
                    mlabel(proto::value(proto::child_c<1>(expr))) {}

                /*
                // Use the ParticleCtx to implement subscripting
                // of a Symbolic expression tree.
                template<typename ParticleType, typename unknown_tuple_type>
                typename proto::result_of::eval<Expr const, ParticleCtx<ParticleType,unknown_tuple_type> const>::type
                eval( typename ParticleType::const_reference particle, unknown_tuple_type unknown_tuple=unknown_tuple_type()) const {
                    ParticleCtx<ParticleType,unknown_tuple_type> const ctx(particle,unknown_tuple);
                    return proto::eval(*this, ctx);
                }
                template<typename ParticleType1, typename ParticleType, typename unknown_tuple_type>
                typename proto::result_of::eval<Expr const, TwoParticleCtx<ParticleType1,ParticleType2,unknown_tuple_type> const>::type
                eval( const typename ParticleType1::double_d& dx, typename ParticleType1::const_reference particle1,  typename ParticleType2::const_reference particle2, unknown_tuple_type unknown_tuple=unknown_tuple_type()) const {
                    TwoParticleCtx<ParticleType1,ParticleType2,unknown_tuple_type> const ctx(dx, particle1, particle2, unknown_tuple);
                    return proto::eval(*this, ctx);
                }
                */


                //mpl::for_each<mpl::range_c<int,1,label_type::number_of_particle_sets::value> > (this->name(proto::as_expr<SymbolicDomain>(expr,label))); \

                template< typename ExprRHS > 
                typename boost::enable_if<
                    mpl::not_<proto::matches<typename proto::result_of::as_expr<ExprRHS>::type, is_not_aliased>>
                    ,mpl::true_>::type
                alias_check(ExprRHS const & expr) const {
                    return mpl::true_();
                }

                template< typename ExprRHS > 
                typename boost::enable_if<
                    proto::matches<typename proto::result_of::as_expr<ExprRHS>::type, is_not_aliased>
                    ,mpl::false_>::type
                alias_check(ExprRHS const & expr) const {
                    return mpl::false_();
                }

                #define DEFINE_THE_OP(name,the_op) \
                template< typename ExprRHS > \
                const SymbolicExpr &operator the_op (ExprRHS const & expr) const { \
                    BOOST_MPL_ASSERT_NOT(( boost::is_same<VariableType,id > )); \
                    check_valid_assign_expr(proto::as_expr<SymbolicDomain>(expr)); \
                    this->name(proto::as_expr<SymbolicDomain>(expr), \
                            proto::matches<typename proto::result_of::as_expr<ExprRHS>::type, is_not_aliased>()); \
                    return *this; \
                } \

                DEFINE_THE_OP(assign,=)
                DEFINE_THE_OP(increment,+=)
                DEFINE_THE_OP(decrement,-=)
                DEFINE_THE_OP(divide,/=)
                DEFINE_THE_OP(multiply,*=)

                

                private:

                //label<mpl::int_<0>,ParticlesType> &mlabel;
                //symbolic<VariableType> &msymbol;
                label_type& mlabel;
                symbol_type& msymbol;

                void copy_from_buffer() const {
                    typedef typename VariableType::value_type value_type;
                    ParticlesType& particles = mlabel.get_particles();
                    std::vector<value_type>& buffer = msymbol.get_buffer(&particles);

                    const size_t n = particles.size();
                    #pragma omp parallel for
                    for (size_t i=0; i<n; i++) {
                        get<VariableType>(particles[i]) = buffer[i];	
                    }
                }

                void post() const {
                    ParticlesType& particles = mlabel.get_particles();

                    if (boost::is_same<VariableType,position>::value) {
                        particles.update_positions();
                    }

                    if (boost::is_same<VariableType,alive>::value) {
                        particles.delete_particles();
                    }
                }

                template< typename ExprRHS>
                typename boost::enable_if<detail::is_univariate<ExprRHS>,void >::type
                check_valid_assign_expr(ExprRHS const & expr) const {
                    const ParticlesType& particles = mlabel.get_particles();
                    const auto& particles_in_expr = fusion::at_c<0>(detail::get_labels()(expr,fusion::nil_())).get_particles();
                    CHECK(&particles_in_expr == &particles, "labels on LHS and RHS of assign expression do not refer to the same particles container");
                }

                template< typename ExprRHS>
                typename boost::enable_if<detail::is_const<ExprRHS>,void >::type
                check_valid_assign_expr(ExprRHS const & expr) const {
                }

                template< typename ExprRHS>
                typename boost::enable_if<detail::is_bivariate<ExprRHS>,void >::type
                check_valid_assign_expr(ExprRHS const & expr) const {
                    static_assert(!detail::is_bivariate<ExprRHS>::value,"asignment expression must be constant or univariate");
                }

                template< typename ExprRHS>
                const SymbolicExpr &assign(ExprRHS const & expr,
                        mpl::false_) const {
                    typedef typename VariableType::value_type value_type;

                    const ParticlesType& particles = mlabel.get_particles();
                    std::vector<value_type>& buffer = msymbol.get_buffer(&particles);
                    buffer.resize(particles.size());

                    const size_t n = particles.size();
                    #pragma omp parallel for
                    for (int i=0; i<n; i++) {
                        EvalCtx<fusion::map<fusion::pair<label_type,typename ParticlesType::const_reference>>> const ctx(fusion::make_map<label_type>(particles[i]));
                        buffer[i] = proto::eval(expr,ctx);
                    }

                    copy_from_buffer();
                    post();

                    return *this;
                }

                template< typename ExprRHS>
                const SymbolicExpr &assign(ExprRHS const & expr,
                        mpl::true_) const {
                    ParticlesType& particles = mlabel.get_particles();

                    const size_t n = particles.size();
                    #pragma omp parallel for
                    for (size_t i=0; i<n; i++) {
                        EvalCtx<fusion::map<fusion::pair<label_type,typename ParticlesType::const_reference>>> const ctx(fusion::make_map<label_type>(particles[i]));

                        get<VariableType>(particles)[i] = proto::eval(expr,ctx);
                    }

                    post();

                    return *this;
                }

                template< typename ExprRHS>
                const SymbolicExpr &increment(ExprRHS const & expr,
                        mpl::false_) const {
                    typedef typename VariableType::value_type value_type;

                    const ParticlesType& particles = mlabel.get_particles();
                    std::vector<value_type>& buffer = msymbol.get_buffer(&particles);
                    buffer.resize(particles.size());

                    const size_t n = particles.size();
                    #pragma omp parallel for
                    for (size_t i=0; i<n; i++) {
                        EvalCtx<fusion::map<fusion::pair<label_type,typename ParticlesType::const_reference>>> const ctx(fusion::make_map<label_type>(particles[i]));
                        buffer[i] = get<VariableType>(particles[i]) + proto::eval(expr,ctx);
                    }

                    copy_from_buffer();
                    post();

                    return *this;
                }

                template< typename ExprRHS>
                const SymbolicExpr &increment(ExprRHS const & expr,
                        mpl::true_) const {
                    ParticlesType& particles = mlabel.get_particles();

                    const size_t n = particles.size();
                    #pragma omp parallel for
                    for (size_t i=0; i<n; i++) {
                        EvalCtx<fusion::map<fusion::pair<label_type,typename ParticlesType::const_reference>>> const ctx(fusion::make_map<label_type>(particles[i]));

                        get<VariableType>(particles)[i] += proto::eval(expr,ctx);
                    }

                    post();

                    return *this;
                }

                template< typename ExprRHS>
                const SymbolicExpr &decrement(ExprRHS const & expr,
                        mpl::false_) const {
                    typedef typename VariableType::value_type value_type;

                    const ParticlesType& particles = mlabel.get_particles();
                    std::vector<value_type>& buffer = msymbol.get_buffer(&particles);
                    buffer.resize(particles.size());

                    const size_t n = particles.size();
                    #pragma omp parallel for
                    for (size_t i=0; i<n; i++) {
                        EvalCtx<fusion::map<fusion::pair<label_type,typename ParticlesType::const_reference>>> const ctx(fusion::make_map<label_type>(particles[i]));
                        buffer[i] = get<VariableType>(particles[i]) - proto::eval(expr,ctx);
                    }

                    copy_from_buffer();
                    post();

                    return *this;
                }

                template< typename ExprRHS>
                const SymbolicExpr &decrement(ExprRHS const & expr,
                        mpl::true_) const {
                    ParticlesType& particles = mlabel.get_particles();

                    const size_t n = particles.size();
                    #pragma omp parallel for
                    for (size_t i=0; i<n; i++) {
                        EvalCtx<fusion::map<fusion::pair<label_type,typename ParticlesType::const_reference>>> const ctx(fusion::make_map<label_type>(particles[i]));

                        get<VariableType>(particles)[i] -= proto::eval(expr,ctx);
                    }

                    post();

                    return *this;
                }

                template< typename ExprRHS>
                const SymbolicExpr &multiply(ExprRHS const & expr,
                        mpl::false_) const {
                    typedef typename VariableType::value_type value_type;

                    const ParticlesType& particles = mlabel.get_particles();
                    std::vector<value_type>& buffer = msymbol.get_buffer(&particles);
                    buffer.resize(particles.size());

                    const size_t n = particles.size();
                    #pragma omp parallel for
                    for (size_t i=0; i<n; i++) {
                        EvalCtx<fusion::map<fusion::pair<label_type,typename ParticlesType::const_reference>>> const ctx(fusion::make_map<label_type>(particles[i]));
                        buffer[i] = get<VariableType>(particles[i]) * proto::eval(expr,ctx);
                    }

                    copy_from_buffer();
                    post();

                    return *this;
                }

                template< typename ExprRHS>
                const SymbolicExpr &multiply(ExprRHS const & expr,
                        mpl::true_) const {
                    ParticlesType& particles = mlabel.get_particles();

                    const size_t n = particles.size();
                    #pragma omp parallel for
                    for (size_t i=0; i<n; i++) {
                        EvalCtx<fusion::map<fusion::pair<label_type,typename ParticlesType::const_reference>>> const ctx(fusion::make_map<label_type>(particles[i]));

                        get<VariableType>(particles)[i] *= proto::eval(expr,ctx);
                    }

                    post();

                    return *this;
                }

                template< typename ExprRHS>
                const SymbolicExpr &divide(ExprRHS const & expr,
                        mpl::false_) const {
                    typedef typename VariableType::value_type value_type;

                    const ParticlesType& particles = mlabel.get_particles();
                    std::vector<value_type>& buffer = msymbol.get_buffer(&particles);
                    buffer.resize(particles.size());

                    const size_t n = particles.size();
                    #pragma omp parallel for
                    for (size_t i=0; i<n; i++) {
                        EvalCtx<fusion::map<fusion::pair<label_type,typename ParticlesType::const_reference>>> const ctx(fusion::make_map<label_type>(particles[i]));
                        buffer[i] = get<VariableType>(particles[i]) * proto::eval(expr,ctx);
                    }

                    copy_from_buffer();
                    post();

                    return *this;
                }

                template< typename ExprRHS>
                const SymbolicExpr &divide(ExprRHS const & expr,
                        mpl::true_) const {
                    ParticlesType& particles = mlabel.get_particles();

                    const size_t n = particles.size();
                    #pragma omp parallel for
                    for (size_t i=0; i<n; i++) {
                        EvalCtx<fusion::map<fusion::pair<label_type,typename ParticlesType::const_reference>>> const ctx(fusion::make_map<label_type>(particles[i]));

                        get<VariableType>(particles)[i] *= proto::eval(expr,ctx);
                    }

                    post();

                    return *this;
                }
            };

    }

}

#endif /* SYMBOLIC_H_ */
