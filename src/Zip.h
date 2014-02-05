/*
* zip.h
*
* Copyright 2013 Martin Robinson
*
* This file is part of PDE_BD.
*
* PDE_BD is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* PDE_BD is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with PDE_BD. If not, see <http://www.gnu.org/licenses/>.
*
* Created on: 8 Mar 2013
* Author: robinsonm
*/

#ifndef ZIP_H_
#define ZIP_H_

#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>
#include <iterator>

template <typename... T>
auto zip(const T&... containers) -> boost::iterator_range<boost::zip_iterator<decltype(boost::make_tuple(std::begin(containers)...))>>
{
    auto zip_begin = boost::make_zip_iterator(boost::make_tuple(std::begin(containers)...));
    auto zip_end = boost::make_zip_iterator(boost::make_tuple(std::end(containers)...));
    return boost::make_iterator_range(zip_begin, zip_end);
}


#endif /* ZIP_H_ */
