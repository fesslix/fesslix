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

#include "flxfunction.h"


PYBIND11_EXPORT FlxFunction* parse_py_para(const std::string& para_name, py::dict config, const bool required=true, const tuint NumbOfPara=0);
PYBIND11_EXPORT const bool parse_py_para_as_bool(const std::string& para_name, py::dict config, const bool required, const bool def_val=false);
PYBIND11_EXPORT const tuint parse_py_para_as_tuint(const std::string& para_name, py::dict config, const bool required, const tuint def_val=0);
PYBIND11_EXPORT const tuint parse_py_para_as_tuintNo0(const std::string& para_name, py::dict config, const bool required, const tuint def_val=1);
PYBIND11_EXPORT const tulong parse_py_para_as_tulong(const std::string& para_name, py::dict config, const bool required, const tulong def_val=1);
PYBIND11_EXPORT const tdouble parse_py_para_as_float(const std::string& para_name, py::dict config, const bool required, const tdouble def_val=ZERO);
PYBIND11_EXPORT const tdouble parse_py_para_as_floatPosNo0(const std::string& para_name, py::dict config, const bool required, const tdouble def_val=ONE);
/**
* @brief returns a flxVec from a Python dict
*
* Unless def_val is returned, internally, a reference to the numpy array is returned.
*/
PYBIND11_EXPORT const flxVec parse_py_para_as_flxVec(const std::string& para_name, py::dict config, const bool required, const flxVec def_val=flxVec(0));
PYBIND11_EXPORT std::string parse_py_para_as_string(const std::string& para_name, py::dict config, const bool required, const std::string def_val="");
PYBIND11_EXPORT std::string parse_py_para_as_word(const std::string& para_name, py::dict config, const bool required, const bool lowercase, const bool emptyAllow=false,const bool numbegallow=false, const std::string def_val="");
PYBIND11_EXPORT void parse_py_para_as_word_lst(std::vector<std::string>& res, const std::string& para_name, py::dict config, const bool required, const bool lowercase, const bool emptyAllow=false,const bool numbegallow=false);

PYBIND11_EXPORT const tuint parse_py_obj_as_tuint(py::object obj, std::string descr);
PYBIND11_EXPORT py::dict parse_py_obj_as_dict(py::object obj, std::string descr);
PYBIND11_EXPORT py::list parse_py_obj_as_list(py::object obj, std::string descr);
/**
* @brief returns a flxVec from a Python object
*
* Internally, a reference to the numpy array is returned.
*/
PYBIND11_EXPORT flxVec parse_py_obj_as_flxVec(py::object obj, std::string descr);
PYBIND11_EXPORT std::string parse_py_obj_as_string(py::object obj, std::string descr);

PYBIND11_EXPORT std::string parse_str_as_word(std::string strV, const bool lowercase, const bool emptyAllow=false,const bool numbegallow=false);
std::vector<std::string> parse_strseq_as_vec(const std::string& strseq, const char sep=',');

FLXLIB_EXPORT FlxFunction* parse_function(const std::string& funStr);
PYBIND11_EXPORT FlxFunction* parse_function(py::object pyobj, std::string descr="", const tuint NumbOfPara=0);


class FlxObjBase;
PYBIND11_EXPORT FlxObjBase* parse_code(py::object pyobj, std::string descr="");

/**
* @brief returns a numpy array that is based on a pointer whose memory is managed externally
*
* @note for any template type, add a line below the definition!!!
*/
template <typename T>
PYBIND11_EXPORT py::array_t<T> py_wrap_array_no_ownership(T* ptr, size_t N);

PYBIND11_EXPORT py::array_t<tdouble> convert_flxVec_to_pyArray(flxVec& vec);

