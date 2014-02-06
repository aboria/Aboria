/*
 * ptr.h
 *
 *  Created on: 30 Jan 2014
 *      Author: mrobins
 */

#ifndef PTR_H_
#define PTR_H_

#include <memory>
#include <boost/shared_ptr.hpp>

namespace Aboria {

//template<class T> using ptr = boost::shared_ptr<T>;
//#define ptr std::auto_ptr
template<class T>
class ptr: public boost::shared_ptr<T> {
public:
	template< class Y >
	explicit ptr( Y* ptr ):boost::shared_ptr<T>(ptr) {}
	ptr( const ptr& r ):boost::shared_ptr<T>(r) {}
};


}


#endif /* PTR_H_ */
