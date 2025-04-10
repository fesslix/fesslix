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

#include "flxfunction.h"

PYBIND11_EXPORT FlxFunction* parse_py_para(const std::string& para_name, py::dict config, const bool required=true);
PYBIND11_EXPORT const bool parse_py_para_as_bool(const std::string& para_name, py::dict config, const bool required, const bool def_val=false);
PYBIND11_EXPORT const tdouble parse_py_para_as_float(const std::string& para_name, py::dict config, const bool required, const tdouble def_val=ZERO);
PYBIND11_EXPORT std::string parse_py_para_as_string(const std::string& para_name, py::dict config, const bool required, const std::string def_val="");
PYBIND11_EXPORT std::string parse_py_para_as_word(const std::string& para_name, py::dict config, const bool required, const bool lowercase, const bool emptyAllow=false,const bool numbegallow=false, const std::string def_val="");
PYBIND11_EXPORT void parse_py_para_as_word_lst(std::vector<std::string>& res, const std::string& para_name, py::dict config, const bool required, const bool lowercase, const bool emptyAllow=false,const bool numbegallow=false);

PYBIND11_EXPORT py::dict parse_py_obj_as_dict(py::object obj, std::string descr);
PYBIND11_EXPORT py::list parse_py_obj_as_list(py::object obj, std::string descr);
PYBIND11_EXPORT std::string parse_py_obj_as_string(py::object obj, std::string descr);

PYBIND11_EXPORT std::string parse_str_as_word(std::string strV, const bool lowercase, const bool emptyAllow=false,const bool numbegallow=false);
std::vector<std::string> parse_strseq_as_vec(const std::string& strseq, const char sep=',');

FLXLIB_EXPORT FlxFunction* parse_function(const std::string& funStr);
PYBIND11_EXPORT FLXLIB_EXPORT FlxFunction* parse_function(py::object pyobj, std::string descr="");


