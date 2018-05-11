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

#include <assert.h>
#include <iostream>
#include <signal.h>

#ifdef NDEBUG
#define ASSERT(condition, message)
#define ASSERT_CUDA(condition)
#else
#define ASSERT(condition, message)                                             \
  if (!(condition)) {                                                          \
    std::cerr << "Assertion `" #condition "` failed in " << __FILE__           \
              << " line " << __LINE__ << ": " << message << std::endl;         \
    raise(SIGTRAP);                                                            \
  }
#define ASSERT_CUDA(condition) assert(condition);

#endif

#define CHECK(condition, message)                                              \
  if (!(condition)) {                                                          \
    std::cerr << "Assertion `" #condition "` failed in " << __FILE__           \
              << " line " << __LINE__ << ": " << message << std::endl;         \
    raise(SIGTRAP);                                                            \
  }

#define CHECK_CUDA(condition, message)                                         \
  if (!(condition)) {                                                          \
    printf("Assertion %s failed in %s line %d message %s\n", #condition,       \
           __FILE__, __LINE__, message);                                       \
    assert(false);                                                             \
  }

#define ERROR(message)                                                         \
  std::cerr << "Error at " << __FILE__ << " line " << __LINE__ << ": "         \
            << message << std::endl;                                           \
  raise(SIGTRAP);

#define ERROR_CUDA(message)                                                    \
  printf("%s\n", message);                                                     \
  assert(false);

// std::exit(EXIT_FAILURE);

#ifndef ABORIA_LOG_LEVEL
#ifdef NDEBUG
#define ABORIA_LOG_LEVEL 1
#else
#define ABORIA_LOG_LEVEL 2
#endif
#endif

#define LOG(level, message)                                                    \
  if (level <= ABORIA_LOG_LEVEL) {                                             \
    char color[] = {0x1b, '[', '3', '8', ';', '5', ';', '7', 'm', 0};          \
    char reset[] = {0x1b, '[', '0', 'm', 0};                                   \
    switch (level) {                                                           \
    case 1:                                                                    \
      color[7] = '4';                                                          \
      break;                                                                   \
    case 2:                                                                    \
      color[7] = '2';                                                          \
      break;                                                                   \
    case 3:                                                                    \
      color[7] = '1';                                                          \
    default:                                                                   \
      color[5] = '5';                                                          \
    }                                                                          \
    std::cout << color << message << reset << std::endl;                       \
  }

#define LOG_CUDA(level, message)                                               \
  if (level <= ABORIA_LOG_LEVEL) {                                             \
    printf("%s\n", message);                                                   \
  }

#define LOG_CUDA1(level, message, variable)                                    \
  if (level <= ABORIA_LOG_LEVEL) {                                             \
    printf("%s %f\n", message, variable);                                                \
  }

#define LOG_CUDA2(level, message, variable1, variable2)                                   \
  if (level <= ABORIA_LOG_LEVEL) {                                             \
    printf("%s %f %f\n", message, variable, variable2);                        \
  }


// char color[] =  { 0x1b, '[', '1', ';', '3', '7', 'm', 0 };

#define LOG_BOLD(level, message)                                               \
  if (level <= ABORIA_LOG_LEVEL) {                                             \
    char bold[] = {0x1b, '[', '1', 'm', 0};                                    \
    char color[] = {0x1b, '[', '3', '8', ';', '5', ';', '7', 'm', 0};          \
    char reset[] = {0x1b, '[', '0', 'm', 0};                                   \
    switch (level) {                                                           \
    case 1:                                                                    \
      color[5] = '4';                                                          \
      break;                                                                   \
    case 2:                                                                    \
      color[5] = '2';                                                          \
      break;                                                                   \
    case 3:                                                                    \
      color[5] = '1';                                                          \
    default:                                                                   \
      color[5] = '5';                                                          \
    }                                                                          \
    std::cout << bold << color << message << reset << std::endl;               \
  }

#endif /* LOG_H_ */
