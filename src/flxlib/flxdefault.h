/* Fesslix - Stochastic Analysis
 * Copyright (C) 2010-2025 Wolfgang Betz
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

  
#include "flxglobal.h"

const tuint DEFAULT_FPC_PREC = 6;
const tuint DEFAULT_FPC_TYPE = 2;
const tdouble DEFAULT_FPC_BVALU = 1e4;
const tdouble DEFAULT_FPC_BVALL = 1e-2;
const bool DEFAULT_FPC_DEL0 = true;
const bool DEFAULT_FPC_DELP = true;
const bool DEFAULT_FLX_PRGBAR = true;

const char DEFAULT_GAUSS_FILE [] = "{default}";
const tuint DEFAULT_GAUSS_MAXNUMB = 45;
const tuint DEFAULT_GAUSS_NUMB = 20;

const tuint DEFAULT_LEGENDRE_NUMB = 10;

const bool DEFAULT_LOG_INPUT=true;
const char DEFAULT_LOG_FILE [] = "fesslix.log";
const bool DEFAULT_LOG_OUTPUT = true;
const int DEFAULT_LOG_LEVEL = 4;
const bool DEFAULT_LOG_TRUNC = true;

const tuint DEFAULT_MT19937_INIT_CALLS = 1000;
const bool DEFAULT_MT19937_INIT_RAND = true;
const bool DEFAULT_MT19937_INIT_SEED = false;
const tuint DEFAULT_MT19937_INIT_SEEDVALUE = 0;

const tdouble DEFAULT_TOL = 1e-14;

