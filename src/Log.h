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


#ifndef LOG_H_
#define LOG_H_

#include <signal.h>
#include <iostream>
#include <assert.h>

#define ASSERT(condition, message) \
    assert(condition)

#define CHECK(condition, message) \
		if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            raise(SIGTRAP); \
        }

#define ERROR(message) \
            std::cerr << "Error at " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            raise(SIGTRAP);

#define ERROR_CUDA(message) \
            printf("%s\n",message); \
            raise(SIGTRAP);



//std::exit(EXIT_FAILURE);

#ifndef ABORIA_LOG_LEVEL
#	ifdef NDEBUG
#		define ABORIA_LOG_LEVEL 1
#	else
#		define ABORIA_LOG_LEVEL 2
#	endif
#endif

#define LOG(level, message) \
    if (level <= ABORIA_LOG_LEVEL) { \
    	std::cout << message << std::endl; \
    }

#define LOG_CUDA(level,message) \
    if (level <= ABORIA_LOG_LEVEL) { \
        printf("%s\n",message); \
    }



#endif /* LOG_H_ */
