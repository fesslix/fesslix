/* Fesslix - Stochastic Analysis
 * Copyright (C) 2010-2026 Wolfgang Betz
 *
 * Fesslix is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Fesslix is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Fesslix.  If not, see <http://www.gnu.org/licenses/>. 
 */

#pragma once

#include "config.h"

#define RETURN_SUCCESS 0
#define RETURN_ERROR 1
#define RETURN_END 2

typedef double tdouble;
typedef tdouble* tdoublePtr;
typedef float tfloat;
typedef unsigned long long tulong;        // for very very long numbers
typedef unsigned int tnlong;        // for node numbering (and dof numbering)
typedef unsigned int tuint;
typedef const tuint* const_tuintPtr;

const tdouble ZERO = 0.;
const tdouble ONE = 1.;
const tdouble PI = 3.14159265358979323846264338327950;
const tdouble GAMMA = 0.5772156649015328606065120900824;   // Euler's constant
const tdouble Y_INFTY_1 = 30*ONE;
const tdouble Y_INFTY_2 = 200*ONE;

#if FLX_DEBUG_COUT
  #include <iostream>
  inline void debug(std::string dstr) { std::cout << dstr << std::endl;}
#endif


