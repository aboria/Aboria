/*
 * ptr.h
 *
 *  Created on: 30 Jan 2014
 *      Author: mrobins
 */

#ifndef PTR_H_
#define PTR_H_

#include <memory>

namespace Aboria {

template<class T> using ptr = std::shared_ptr<T>;

}


#endif /* PTR_H_ */
