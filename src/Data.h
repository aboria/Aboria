/*
 * Data.h
 *
 * Copyright 2012 Martin Robinson
 *
 * This file is part of RD_3D.
 *
 * RD_3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * RD_3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with RD_3D.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 25 Oct 2012
 *      Author: robinsonm
 */

#ifndef DATA_H_
#define DATA_H_


#include <boost/preprocessor.hpp>

#ifndef DATA_container
	#define DATA_container std::vector
#endif
#define DATA_reference(n)                          DATA_container<BOOST_PP_SEQ_ELEM(n,DATA_types)>::reference
#define DATA_param(z, n, notused)                  DATA_reference(n) BOOST_PP_SEQ_ELEM(n,DATA_names)
#define DATA_declarations(z, n, notused)           DATA_param(z, n, notused);
#define DATA_constructor_assignment(z, n, notused) BOOST_PP_SEQ_ELEM(n,DATA_names)(BOOST_PP_SEQ_ELEM(n,DATA_names))
#define DATA_assignment_operator(z, n, notused)    BOOST_PP_SEQ_ELEM(n,DATA_names) = rhs.BOOST_PP_SEQ_ELEM(n,DATA_names);
#define DATA_n                                     BOOST_PP_SEQ_SIZE(DATA_names)

#define DATA_entry
struct BOOST_PP_CAT(DATA_typename,_Entry) {
	BOOST_PP_CAT(DATA_typename,_Entry)(BOOST_PP_ENUM(DATA_n,DATA_param, ~)):
		BOOST_PP_ENUM(DATA_n, DATA_constructor_assignment, ~) {}
	BOOST_PP_CAT(DATA_typename,_Entry)& operator=(const BOOST_PP_CAT(DATA_typename,_Entry) &rhs) {
		if (this != &rhs) {
			BOOST_PP_REPEAT(DATA_n, DATA_assignment_operator, ~)
		}
		return *this;
	}
	BOOST_PP_REPEAT(DATA_n, DATA_declarations, ~)
};

#define DATA_push_back_params(z, n, notused)  const BOOST_PP_SEQ_ELEM(n,DATA_types)& BOOST_PP_CAT(_,BOOST_PP_SEQ_ELEM(n,DATA_names))
#define DATA_push_back_impl(z, n, notused)    BOOST_PP_SEQ_ELEM(n,DATA_names).push_back(BOOST_PP_CAT(_,BOOST_PP_SEQ_ELEM(n,DATA_names)));
#define DATA_pop_back_impl(z, n, notused)     BOOST_PP_SEQ_ELEM(n,DATA_names).pop_back();
#define DATA_clear_impl(z, n, notused)        BOOST_PP_SEQ_ELEM(n,DATA_names).clear();
#define DATA_data_index(z, n, notused)        BOOST_PP_SEQ_ELEM(n,DATA_names)[i]
#define DATA_define_containers(z, n, notused) std::vector<BOOST_PP_SEQ_ELEM(n,DATA_types)> BOOST_PP_SEQ_ELEM(n,DATA_names);

#define DATA_data
class DATA_typename {
public:
	DATA_typename() {}
	void push_back(BOOST_PP_ENUM(DATA_n,DATA_push_back_params, ~)) {
		BOOST_PP_REPEAT(DATA_n, DATA_push_back_impl, ~)
	}

	void pop_back() {
		BOOST_PP_REPEAT(DATA_n, DATA_pop_back_impl, ~)
	}

	void clear() {
		BOOST_PP_REPEAT(DATA_n, DATA_clear_impl, ~)
	}

	BOOST_PP_CAT(DATA_typename,_Entry) operator[](const int i) {
		return BOOST_PP_CAT(DATA_typename,_Entry)(BOOST_PP_ENUM(DATA_n,DATA_data_index, ~));
	}

	size_t size() const {
		return r.size();
	}

	BOOST_PP_REPEAT(DATA_n, DATA_define_containers, ~)
};

DATA_entry
DATA_data

#undef DATA_typename
#undef DATA_names
#undef DATA_types
#undef DATA_container
#undef DATA_reference
#undef DATA_param
#undef DATA_declarations
#undef DATA_constructor_assignment
#undef DATA_assignment_operator
#undef DATA_n
#undef DATA_push_back_params
#undef DATA_push_back_impl
#undef DATA_pop_back_impl
#undef DATA_clear_impl
#undef DATA_data_index
#undef DATA_define_containers
#undef DATA_entry
#undef DATA_data


#endif /* DATA_H_ */
