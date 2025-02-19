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




#include <pybind11/pybind11.h>

namespace py = pybind11;

double add(double a, double b) {
    return a + b;
}

PYBIND11_MODULE(fesslix, m) {
    m.def("add", &add, "A function that adds two numbers");
    m.attr("the_answer") = 42;
}


