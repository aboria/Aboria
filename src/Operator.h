/*
 * Operator.h
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
 *  Created on: 30 Oct 2012
 *      Author: robinsonm
 */

#ifndef OPERATOR_H_
#define OPERATOR_H_

#include <boost/timer/timer.hpp>
#include <memory>
#include <vector>

#include "Ptr.h"

namespace Aboria {

class Operator {
	struct OperatorConcept {
		virtual ~OperatorConcept() {}
		virtual void execute() = 0;
		virtual void reset() = 0;
		virtual void print(std::ostream& out) const = 0;
	};

	template<typename T>
	struct OperatorModel : OperatorConcept {
		OperatorModel( const T& t ) : object( t ) {}
		virtual ~OperatorModel() {}
		virtual void reset() {object.reset();}
		virtual void execute() {object.execute();}
		virtual void print(std::ostream& out) const {object.print(out);}
	private:
		T object;
	};

	std::shared_ptr<OperatorConcept> op;

public:
	template<typename T>
	Operator( const T& obj ) :
		op( new OperatorModel<T>( obj ) ) {
		reset();
	}

	void execute();
	void reset();
	void print(std::ostream& out) const;
	friend std::ostream& operator<<( std::ostream& out, const Operator& b ) {
		b.print(out);
		return out;
	}
private:
	void resume_timer();
	void stop_timer();

	boost::timer::cpu_timer timer;
	double total_time;
};


class OperatorList: public std::vector<ptr<Operator> > {
public:

	void print(std::ostream& out) const {
		out << "List of "<<this->size()<< " operators:"<< std::endl;
		for (auto i : *this) {
			out << "\t" << i <<std::endl;
		}
		out << "End list of "<<this->size()<< " operators";

	}
	void execute() {
		for (auto i : *this) {
			i->execute();
		}
	}
	void reset() {
		for (auto i : *this) {
			i->reset();
		}
	}
};

template <typename FunctionType>
class RepeatOperator {
public:
	RepeatOperator(ptr<Operator> op, FunctionType n_func): op(op),n_func(n_func) {}
	void print(std::ostream& out) const {
		out << "Repeat Operator (n="<<n_func()<<"): "<<op<<std::endl;
	}
	void execute() {
		const unsigned int n = n_func();
		for (int i = 0; i < n; ++i) {
			op->execute();
		}
	}
	void reset() {
		op->reset();
	}
private:
	ptr<Operator> op;
	FunctionType n_func;
};

template<FunctionType>
ptr<Operator> repeat(ptr<Operator> op, FunctionType f) {
	return ptr<Operator>(new Operator(
			RepeatOperator<FunctionType>(op,f)
	));
}


}

#endif /* OPERATOR_H_ */
