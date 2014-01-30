/*
 * Operator.cpp
 *
 *  Created on: 30 Oct 2012
 *      Author: robinsonm
 */

#include "Operator.h"
#include <sstream>
#include <string>
#include "Log.h"

namespace Aboria {

template <typename T>
const std::string to_string(const T& data)
{
   std::ostringstream conv;
   conv << data;
   return conv.str();
}

void Operator::execute() {
	resume_timer();
	LOG(2, "Starting Operator: " << *this);

	op->execute();

	LOG(2, "Stopping Operator: " << *this);
	stop_timer();
}

void Operator::reset() {
	total_time = 0;
	op->reset();
}

void Operator::print(std::ostream& out) const {
	op->print(out);
	out << " (Time to execute: " + to_string(total_time) + " s)";
}

void Operator::resume_timer() {
	timer.start();
}

void Operator::stop_timer() {
	timer.stop();
	//global_timer.stop();
	const double time = (timer.elapsed().user + timer.elapsed().user)/double(1000000000);
	total_time += time;
}
}




