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


#ifndef DATAVECTOR_H_
#define DATAVECTOR_H_

#include "Particles.h"

namespace Aboria {

/// A thin wrapper class around Particles that acts like a vector
/// to the variable given by \p T
/// \param T the Variable type
/// \param ParticlesType the type of the Particles container. Note that
/// ParticlesType must have been templated on T
template<typename T, typename ParticlesType>
class DataVector {
public:
    /// A typedef of T, the Variable type
	typedef T variable_type;
    /// The value_type of the Variable
	typedef typename T::value_type value_type;

    /// pass in the particles continer for the data
	DataVector(ParticlesType& particles):
		particles(particles)
	{};

    /// \return a reference to the parent particles container
	ParticlesType &get_particles() {
		return particles;
	}

    /// \return a const reference to the parent particles container
	const ParticlesType &get_particles() const {
		return particles;
	}

    /// returns the size of the vector
	std::size_t size() const {
		return particles.size();
	}

    /// index the vector using index \p i. This returns
    /// the same as get<T>(particles[i])
    value_type operator[](const unsigned int i) const {
        return get<variable_type>(particles[i]);
    }
protected:
	ParticlesType &particles;
};


}
#endif /* DATAVECTOR_H_ */
