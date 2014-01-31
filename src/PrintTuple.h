/*
 * PrintTuple.h
 *
 *  Created on: 31 Jan 2014
 *      Author: mrobins
 */

#ifndef PRINTTUPLE_H_
#define PRINTTUPLE_H_

namespace Aboria {
// Define a type which holds an unsigned integer value
template<std::size_t> struct int_{};

template <class Tuple, size_t Pos>
std::ostream& print_tuple(std::ostream& out, const Tuple& t, int_<Pos> ) {
	out << std::get< std::tuple_size<Tuple>::value-Pos >(t) << ',';
	return print_tuple(out, t, int_<Pos-1>());
}

template <class Tuple>
std::ostream& print_tuple(std::ostream& out, const Tuple& t, int_<1> ) {
	return out << std::get<std::tuple_size<Tuple>::value-1>(t);
}

template <class... Args>
std::ostream& operator<<(std::ostream& out, const std::tuple<Args...>& t) {
	out << '(';
	print_tuple(out, t, int_<sizeof...(Args)>());
	return out << ')';
}
}

#endif /* PRINTTUPLE_H_ */
