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

#include "Species.h"
//#include <boost/timer.hpp>
#include <boost/timer/timer.hpp>
#include <initializer_list>

namespace Tyche {
template<typename InputType, typename OutputType>
class BaseOperator {
	BaseOperator(InputType& input, OutputType& output):
		input(input),output(output) {}
	virtual ~BaseOperator() {};

	void execute();
	void reset();
	friend std::ostream& operator<<( std::ostream& out, const Operator& b ) {
		b.print(out);
		return out;
	}
	bool get_active() { return active;};
	void set_active(bool a) { active = a; };
	const InputType& get_input() {return input;}
	const OutputType& get_output() {return output;}

	protected:
		virtual void reset_impl();
		virtual void execute_impl();
		virtual void print_impl(std::ostream& out) const;

	private:
		void resume_timer();
		void stop_timer();

		boost::timer::cpu_timer timer;
		double total_time;
		bool active;
		InputType& input;
		OutputType& output;
	};
};

class Operator {
	struct OperatorConcept {
		virtual ~OperatorConcept() {}

		void execute();
		void reset();
		void print(std::ostream& out);
		friend std::ostream& operator<<( std::ostream& out, const Operator& b ) {
			b.print(out);
			return out;
		}
		bool get_active() { return active;};
		void set_active(bool a) { active = a; };


	protected:
		virtual void reset_impl() = 0;
		virtual void execute_impl() = 0;
		virtual void print_impl(std::ostream& out) const = 0;

	private:
		void resume_timer();
		void stop_timer();

		boost::timer::cpu_timer timer;
		double total_time;
		bool active;
	};

	template<typename T>
	struct OperatorModel : OperatorConcept {
		OperatorModel( const T& t ) : object( t ) {}
		virtual ~ObjectModel() {}
	protected:
		virtual void reset_impl() {object.reset();}
		virtual void execute_impl() {object.execute();}
		virtual void print_impl(std::ostream& out) {object.print();}
	private:
		T object;
	};

	boost::shared_ptr<OperatorConcept> op;

public:
	template<typename T>
	Operator( const T& obj ) :
		Operator( new OperatorModel<T>( obj ) ) {}
};


class OperatorList {
public:
	OperatorList() {}
	OperatorList(Operator& o) {
		list.push_back(&o);
	}
	OperatorList(std::initializer_list<Operator*> arg) {
		for (auto i: arg) {
			list.push_back(i);
			for (auto s: i->get_species()) {
				add_species(*s);
			}
		}

	}
	virtual ~OperatorList() {}

	static std::auto_ptr<Operator> New() {
		return std::auto_ptr<Operator>(new OperatorList());
	}

	void push_back(Operator* const i) {
		list.push_back(i);
		for (auto s: i->get_species()) {
			add_species(*s);
		}
	}
	void push_back(const OperatorList& i) {
		for (Operator* j: i.list) {
			list.push_back(j);
			for (auto s: j->get_species()) {
				add_species(*s);
			}
		}
	}

//	OperatorList& operator+=(const OperatorList &rhs) {
//		for (Operator* j: rhs.list) {
//			list.push_back(j);
//			for (auto s: j->get_species()) {
//				add_species(*s);
//			}
//		}
//		return *this;
//	}
//	OperatorList& operator+=(Operator &rhs) {
//		list.push_back(&rhs);
//		for (auto s: rhs.get_species()) {
//			add_species(*s);
//		}
//		return *this;
//	}
	std::vector<Operator*>& get_operators() {return list;}

protected:
	virtual void print(std::ostream& out) const {
		out << "List of "<<list.size()<< " operators:"<< std::endl;
		for (auto i : list) {
			out << "\t" << *i << " ("<<i->get_time_string()<<")"<<std::endl;
		}
		out << "End list of "<<list.size()<< " operators";

	}
	virtual void integrate(const double dt) {
		for (auto i : list) {
			i->operator ()(dt);
		}
	}
	virtual void reset_execute() {
		for (auto i : list) {
			i->reset();
		}
	}

	virtual void add_species_execute(Species &s) {
		for (auto i : list) {
			i->add_species(s);
		}
	}

	std::vector<Operator*> list;
};


//OperatorList operator+(Operator& arg1, Operator& arg2);
//OperatorList operator+(Operator& arg1, OperatorList& arg2);
//OperatorList operator+(OperatorList& arg1, Operator& arg2);
//OperatorList operator+(OperatorList& arg1, OperatorList& arg2);


}
#endif /* OPERATOR_H_ */
